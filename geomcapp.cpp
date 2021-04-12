/*
 * 	test driver for capp geometry
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "GeomCAPP.h"
#include "PlotGeom.h"

int main()
{
  int prlev = 0;
  GeomCAPP *geom = new GeomCAPP();
  geom->setPrlev(0);
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
  int nzp = geom->getNz();
  printf("== ntri=%d  nzp=%d\n",ntri,nzp);
  /*  check FVM nodes */
  int *chtri = geom->Chtri();
  int *ispp = geom->Ispp();
  int *jspowz = geom->Ispowz();
  double totvol = 0.0;
  double powvol = 0.0;
  int npowvol = 0;
  for(int t=0; t < ntri; t++) {
    int c = chtri[t];
    for(int z=0; z < nzp; z++) {
      int tz = t*nzp+z;
            /* 5 faces for a tri-prism  */
      int tb[5];
      double area[5],dist[5];
      double vol = geom->getFVM(t,z,tb,area,dist);
      if(prlev) printf("t%d z%d v=%.2le:",t,z,vol);
      for(int e=0; e < 5; e++) {
        int etb = tb[e];
        if(etb >=0) {
          double distL = geom->getDistL(etb,tz);
          if(prlev) printf(" %.3le(%d)",distL,etb);
        }
      }
      if(prlev) printf("\n");
      if(c >= 0) {
        if(ispp[c] && jspowz[z]) {
          powvol += vol;
          npowvol += 1;
        }
      }
      totvol += vol;
    }
  }
  double chanarea = geom->getChanArea();
  printf(" Total vol=%.3le  powvol=%.3le %d vols %.3le ChanArea=%.2le\n",
	totvol,powvol, npowvol,powvol/npowvol, chanarea);
  // geom->printMat3D();
// #define PLOT
#ifdef PLOT
  PlotGeom *plot = new PlotGeom();
  plot->plot(geom,"capp.png");
  delete plot;
#endif
  delete geom;
};
