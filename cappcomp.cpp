/*
 *	cappdyn.cpp
 *	test with delayed neutron
 *	and thermal hydraulic feedback
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <chrono>
#include <petscksp.h>
#include "CoreDesc.h"
#include "GeomCAPP.h"
#include "XScapp.h"
#include "NeutFVM.h"
#include "Delay.h"

#include "Cond3D.h"
#include "SegCond.h"

#include "MPupdate.h"

int main(int argc, char **argv)
{
  /* initial protocol to call petsc */
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,(char*)0,NULL);
  if (ierr) {
    printf("*** problem during MPI initialization.\n");
    return ierr;
  }
  auto start = std::chrono::high_resolution_clock::now();

  int mpid;
  MPI_Comm_rank (MPI_COMM_WORLD, &mpid);

  int prlev = 0;

  CoreDesc desc;
  double Totpow = 15.0e+6;      // (Watt)
  double Psys = 3.0e+6;         // (Pa)
  double T0 = 273.15;
  double Tin =  300.0 +T0;      // 573.15 K
  double Tout = 630.0 +T0;	// 913.15 K
  desc.setMMR(Totpow,Psys,Tin,Tout);
  // if(mpid == 0) desc.print();
  double Mdot = desc.Mdot; // 8.75;  // (kg/s)
  GeomCAPP *geom = new GeomCAPP();
  int mat[]  = {0, 1, 0,1, 1,1,1, 0,0,0,0, -1,0,0,0,0};
  double pitch = 0.322;  // (meter)
  double side = pitch/sqrt(3.0);
  geom->setPlane(side, 6, mat);
  int chan[] = {0, 1, 0,1, 1,1,1, 0,0,0,0, -1,0,0,0,0};
  geom->setChannel(chan);
  /*  set control rod number  */
  int crod[] = {1, 0, 2,0, 0,0,0, 0,0,0,0, -1,0,0,0,0};
  geom->setCrod(crod);

  /* units in meter  */
  double dz[] = {0.350,0.700,0.01736,0.66528,0.03472,0.66528,0.03472,0.66528,
        0.03472,0.66528,0.03472,0.66528,0.03472,0.66528,0.01736,0.70,0.35};
  int zseg[]  = {   2,   4,    1,     4,    1,     4,    1,     4,
            1,     4,    1,     4,    1,     4,    1,   4,   2};
  int ispowz[] = {0,0,0, 1,0,1,0,1, 0,1,0,1,0,1,0,0,0};
  int nz = sizeof(dz)/sizeof(double);
  int nz2 = sizeof(ispowz)/sizeof(int);
  assert(nz==nz2);
  geom->extrude(nz,zseg,dz,ispowz);
  /* set material index for fuel assemblies */
  /* 0: graphite,  1: gra. w. hole,  2: fuel block  */
  int matz[] = {1,1,1, 2,1,2,1,2,1, 2,1,2,1,2, 1,1,1};
  int nz3 = sizeof(matz)/sizeof(int);
  assert(nz3==nz);
  geom->setHexMat(1,matz);
  geom->setHexMat(3,matz);
  geom->setHexMat(4,matz);
  geom->setHexMat(5,matz);
  geom->setHexMat(6,matz);

  int ntri = geom->getNtri();

  /*	prepare Neutronics  */
  XScapp *xs = new XScapp();
  xs->setGeom(geom);
  int efpd = 50;
  double Tmref = 800;
  double Tfref = 800;
  xs->open("MMR1g.db",efpd,Tmref,Tfref);
  if(mpid==0) printf("# Totpow=%.1le efpd=%d\n",Totpow,efpd);

  NeutFVM *neut = new NeutFVM();
  neut->setGeom(geom);
  neut->setXS(xs,Tmref,Tfref);
  neut->prepare(1);     // set no. of neutron groups

  /*	prepare TH  */
  Cond3D *cond = new Cond3D();
  cond->setGeom(geom);
  cond->setFluence(xs->getFlu());
  double six = 6;

  THpro *thpro = new THpro();
  double fkweit[] = {1.0,0.9,0.8};      // multipliers for conductivity
  double fcweit[] = {1.0,0.9,0.9};      // multipliers for heat capacity
  thpro->setFactors(3,fkweit,fcweit);
  double Dcool=1.40e-02; double wetden=desc.wet_den;
  double ncool= desc.nchan;       // number of coolant holes
  thpro->setChannel(Dcool,wetden,ncool);
  thpro->setMdot(Mdot,Psys);

  cond->setTHpro(thpro);
  /*  we give 1/6 core flow rate, power  */
  cond->setMdot(Mdot/six);
  double corevol = geom->getPowVol();
  cond->setTin(Tin);
  cond->setup();

  double Mdotcp = Mdot * 5195.0;

  int npows = geom->getNpows();
  if(mpid==0) printf(" npows=%d\n",npows);
  /*	segment	*/
  SegCond *seg = new SegCond();
  seg->setGeom();
  seg->setup(npows);
  seg->setPfactor(4.5); // multiplier to convert core power density to compact

  /*  prepare Thermal property  */

  MPupdate *mpupd = new MPupdate();
  /* system pressure  */
  mpupd->set(3.0e+6);

  /*  triso dimension (in m) */
  double dKer = 500.0e-6;
  double tlay[] = {100.0,40.0,35.0,40.0};
  double PF = 0.4;
  mpupd->setTriso(dKer,4,tlay,PF);
  /*  prepare vector storages  */
  Vec powdenvec;
  ierr = VecCreate(PETSC_COMM_WORLD,&powdenvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(powdenvec,PETSC_DECIDE,npows);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(powdenvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(powdenvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  Vec Tfuelvec;
  ierr = VecDuplicate(powdenvec,&Tfuelvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  /*  prepare flu and fima in vector form  */
  Vec fluvec,fimavec;
  VecDuplicate(powdenvec,&fluvec);
  VecDuplicate(powdenvec,&fimavec);
  VecCopy(cond->contractPdouble(xs->getFlu()),fluvec);
  VecCopy(cond->contractPdouble(xs->getFima()),fimavec);
  mpupd->init(npows);

  /*	solve eigenvalue problem  */
  double keff;
  neut->setEigShift(0.7);
  double sigrod[] = {0.0046,0.0046};
  neut->setRodXS(2, sigrod);
  double rodpos[] = {0.0,3.0};
  neut->setRodPos(2,rodpos);
  neut->setupAB();
  Vec powden,hcoefvec,Tbulkvec;
  Vec Tcomp,DTker;
  double errth = 999.9;
  double errseg = 999.9;
  double *TK,*TF;
  double tolTemp = 0.01;
  if(mpid==0) printf(" Total vol=%.2le power=%.3le (W)  Mdot=%.2le Tin=%.1lf dT=%.2lf K\n",
        corevol,Totpow,Mdot/six,Tin,Totpow/Mdotcp);

  for(int iter=0; iter < 10; iter++) {
    keff = neut->solveEigen();
    if(mpid==0) printf("== iter=%d keff=%.5lf",iter,keff);
    neut->makePower(Totpow/six);
    powden = neut->getPowden();
    cond->steady(powden);
    double Tpeak = cond->getTpeak();
    double TPavg = cond->getTPavg();
    double Texit = cond->getTexit();
    if(mpid==0) printf("  %.2lf %.2lf %.2lf",Tpeak,TPavg,Texit);

    hcoefvec = cond->getPhcoef();
    Tbulkvec = cond->getPTbulk();
    errth = cond->updateTemp();
    seg->setHcoef(hcoefvec);
    seg->steady(Tbulkvec,powden);
    double TSpeak = seg->getTpeak();
    double TSpin = seg->getTpin();
    if(mpid==0) printf("  %.2lf %.2lf",TSpeak,TSpin);

    Tcomp = seg->getTpinVec();

    errseg = mpupd->updateProp(seg,NULL,fluvec,fimavec);
    ierr = VecCopy(Tcomp,Tfuelvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);

    TK = cond->getTK();
    TF = cond->expandPvec(Tfuelvec);
    neut->setTemp(TK,TF);
    if(mpid==0) printf( "  errth=%.2le errseg=%.2lf", errth,errseg);
    if(mpid==0) printf("\n");

     if((errth < tolTemp) && (errseg < tolTemp)) break;
  }
  if(mpid==0) printf("\n== keff=%.5lf errth=%.1le errseg=%.1le\n",
	keff,errth,errseg);

  /*	prepare transient  */
  double dlam[] = {1.33373E-02, 3.27241E-02, 1.20789E-01,
                   3.03011E-01, 8.50566E-01, 2.85639E+00};
  double dbeta[]= {2.33250E-04, 1.20779E-03, 1.15433E-03,
                   2.59444E-03, 1.07174E-03, 4.48447E-04};
  double gvel[] = {1.042574e+4};
  int nprec = sizeof(dbeta)/sizeof(double);
  double beta = 0;
  for(int d=0; d < nprec; d++) beta += dbeta[d];
  neut->start(Totpow/six,keff,beta,gvel);
  Vec powvec = neut->makePower(Totpow/six);

  Vec fissvec = neut->getFiss();
  Delay *delay = new Delay();
  delay->setParam(npows,nprec,dbeta,dlam);
  delay->init(fissvec);
  Vec dnsrc = delay->getNeutron();
  neut->setDsrc(dnsrc);
  neut->setLTE(1.0e-4,1.0e-10);

  /*   TH  */
  cond->start();
  seg->start();

  double endtime = 1000.0;
  double dtmin = 1.0e-6;
  double dtmax = 50.0;
  double dt = 1.0e-6;
  double time = 0;
  /*  control rod movement  */
  rodpos[0] = 0.0;
  rodpos[1] = 3.1;
  neut->setRodPos(2,rodpos);	// rod in 
  printf("# time   dt    nLTE totpow peakden  ");
  printf(" cLTE Tpeak TPavg Texit");
  printf(" sLTE TSpeak TPin");
  printf(" (errth errseg)\n");
  for(int nt=0; (time < endtime) && (nt < 1000); nt++) {
    double nLTE = neut->step(dt);
    neut->march();
    time += dt;
   if(mpid==0)  printf("%.3le %.1le",time,dt);
    double peakden = neut->getPeakPowDen();
    double totpow = neut->getTotPow();
    if(mpid==0) printf("  %.1le %.4le %.2le",nLTE,totpow,peakden);
    fissvec = neut->getFiss();
    delay->step(dt,fissvec);
    delay->march();
    dnsrc = delay->getNeutron();
    neut->setDsrc(dnsrc);

    powden = neut->getPowden();

    double cLTE = cond->step(dt,powden);
    cond->march();
    double Tpeak = cond->getTpeak();
    double TPavg = cond->getTPavg();
    double Texit = cond->getTexit();
    if(mpid==0) printf("  %.1le %.2lf %.2lf %.2lf",cLTE,Tpeak,TPavg,Texit);

    hcoefvec = cond->getPhcoef();
    Tbulkvec = cond->getPTbulk();
    errth = cond->updateTemp();

    seg->setHcoef(hcoefvec);
    double sLTE = seg->step(dt,Tbulkvec,powden);
    seg->march();
    double TSpeak = seg->getTpeak();
    double TSpin = seg->getTpin();
    if(mpid==0) printf("  %.1le %.2lf %.2lf",sLTE,TSpeak,TSpin);

    Tcomp = seg->getTpinVec();

    errseg = mpupd->updateProp(seg,NULL,fluvec,fimavec);

    if(mpid==0) printf(" (%.1le %.1le)", errth,errseg);

    ierr = VecCopy(Tcomp,Tfuelvec);     CHKERRABORT(PETSC_COMM_SELF,ierr);

    TK = cond->getTK();
    TF = cond->expandPvec(Tfuelvec);
    neut->setTemp(TK,TF);


    if(mpid==0) printf("\n");
    /*  control time step  */
    double rLTE = nLTE;
    if(rLTE < cLTE) rLTE = cLTE;
    if(rLTE < sLTE) rLTE = sLTE;

    double hstar = dt*pow(fabs(rLTE),-0.333333);
    if(hstar < 0.7*dt) dt = hstar;
    if(hstar > 1.5*dt) dt = hstar;
    if(dt < dtmin) dt = dtmin;
    if(dt > dtmax) dt = dtmax;
    double remtime = endtime-time;
    if(dt > remtime) dt = remtime;
  }

  delete geom;
  delete xs;
  delete neut;
  delete cond;
  delete thpro;
  delete seg;
  delete mpupd;
  auto stop = std::chrono::high_resolution_clock::now();
  double xtime =std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
  if(mpid==0) printf("== execution %.3lf sec\n",xtime/1000);

};
