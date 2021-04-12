/*
 * 	Geom.cpp
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "Geom.h"

Geom::Geom()
{
  prlev = 0;
  tri = NULL;
  px = NULL;
  py = NULL;
  dz = NULL;
  zpos = NULL;
  coner = NULL;
  top = NULL;
  nchan = 0;
  tchan = NULL;
  chtri = NULL;
  ispp = NULL;
  ispowz = NULL;
  zpow = NULL;
  /*  boundary triangles */
  nbdry = 0;
  btri = NULL;
  bedge = NULL;
  /*  material */
  mat2d = NULL;
  mat3d = NULL;
  /*  contrl rod */
  crod2d = NULL;
  /*  zseg  */
  zseg = NULL;
  /*  FVM  */
  nbtri = NULL;
  FVM_area = NULL;
  FVM_edge = NULL;
  FVM_dist = NULL;
};	// Geom::Geom

Geom::~Geom()
{
  delete [] tri;
  delete [] px;
  delete [] py;
  delete [] dz;
  delete [] zpos;
  delete [] coner;
  delete [] top;
  delete [] tchan;
  delete [] chtri;  // memory ??
  chtri = NULL;
  delete [] ispp;
  delete [] ispowz;
  delete [] zpow;
  delete [] btri;
  delete [] bedge;
  delete [] mat2d;
  delete [] mat3d;	// memory ??
  mat3d = NULL;
  delete [] crod2d;
  delete [] zseg;
  /*  FVM  */
  delete [] nbtri;
  delete [] FVM_area;
  delete [] FVM_edge;
  delete [] FVM_dist;
};	// Geom::~Geom

void Geom::setPrlev(int prl)
{
  prlev = prl;
};	// Geom::setPrlev

int Geom::extrude(int nzinp_i, int zseg_i[], double dz_i[], int ispowz_i[])
{
  /* we expand zone for multiplicity */
  nzinp = nzinp_i;
  zseg = new int[nzinp];
  int nz_exp = 0;
  for(int zi=0; zi < nzinp; zi++) {
    assert(zseg_i[zi] > 0);
    zseg[zi] = zseg_i[zi];
    nz_exp += zseg[zi];
  }
  double *dz_e = new double[nz_exp];
  int *ispowz_e = new int[nz_exp];
  int z = 0;
  for(int zi=0; zi < nzinp; zi++) {
    for(int i=0; i < zseg[zi]; i++) {
      dz_e[z] = dz_i[zi]/zseg[zi];
      ispowz_e[z] = ispowz_i[zi];
      z++;
    }
  }
  assert(z == nz_exp);
  if(prlev) {
    printf(" z %d input expanded to %d\n",nzinp,nz_exp);
    printf("dz:");
    for(int i=0; i < nz_exp; i++) printf(" %.2lf",dz_e[i]);
    printf("\n");
    printf("ispp:");
    for(int i=0; i < nz_exp; i++) printf(" %d",ispowz_e[i]);
    printf("\n");
  }
  int nptsx = extrude(nz_exp,dz_e,ispowz_e);
  delete [] dz_e;
  delete [] ispowz_e;
  return nptsx;
};	// Geom::extrude

int Geom::extrude(int nz_i, double dz_i[], int ispowz_i[])
{
  /* ispowz_i : axial power production flag.	*/
  /*	if NULL assume all 1 (powprod)		*/
  nz = nz_i;
  dz = new double[nz];
  zpos = new double[nz+1];
  for(int k=0; k < nz; k++) {
    dz[k] = dz_i[k];
    zpos[k+1] = zpos[k] + dz[k];
  }
  /* check if p is a coner point */
  coner = new int[npts2d];
  for(int p=0; p < npts2d; p++) coner[p] = 0;
  for(int t=0; t < ntri; t++) {
    for(int c=0; c < 3; c++) {
      int p = tri[t*6+c];
      assert(p < npts2d);
      coner[p] = 1;
    }
  };
  /*  assign 3D point number  */
  top = new int[npts2d];
  int k = 0;
  for(int p=0; p < npts2d; p++) {
    top[p] = k;
    if(coner[p]) k += 2*nz+1;
    else k += nz+1;
  }
  npts3d = k;
  if(prlev) printf("npts3d=%d\n",npts3d);
  /*  possible power production axe */
  ispowz = new int[nz];
  if(ispowz_i == NULL) for(int k=0; k < nz; k++) ispowz[k] = 1;
  else for(int k=0; k < nz; k++) ispowz[k] = ispowz_i[k];
  prepareTH();
  /*  material  */
  /*  we assume axial material zone is uniform !!  */
  mat3d = new int[ntri*nz];
  for(int t=0; t < ntri; t++) {
    for(int z=0; z < nz; z++) mat3d[t*nz+z] = mat2d[t];
  }
  return npts3d;
};	// Geom::extrude

void Geom::prepareTH()
{
  /*  count number of power production element  */
  int kf = 0;
  for(int c=0; c < nchan; c++) {
    if(ispp[c]) {
      for(int z=0; z < nz; z++) {
        if(ispowz[z]) kf++;
      }
    }
  }
  npows = kf;
  /* get number of power generating axial sections */
  nzpow = 0;
  for(int z=0; z < nz; z++) if(ispowz[z]) nzpow++;
  /*  set mid point elevation  */
  zpow = new double[nzpow];
  double z0 = 0.0;
  int k = 0;
  for(int z=0; z < nz; z++) {
    double z1 = z0+dz[z];
    if(ispowz[z]) {
      zpow[k] = (z0+z1)/2;
      k++;
    }
    z0 = z1;
  }
  assert(k == nzpow);
  if(prlev) {
    printf("zpow:");
    for(k=0; k < nzpow; k++) printf(" %.2lf",zpow[k]);
    printf("\n");
    if(prlev > 5) exit(0);
  }
  
  /*   channel no. to tri  */
  assert(chtri == NULL);
  chtri = new int[ntri];
  for(int t=0; t < ntri; t++) chtri[t] = -1;
  for(int c=0; c < nchan; c++) {
    int t = tchan[c];
    assert(t < ntri);
    chtri[t] = c;
    assert((t >= 0) && ( t < ntri));
  }
  if(prlev > 1) {
    printf("chtri: ntri=%d nchan=%d npows=%d:\n",ntri,nchan,npows);
    for(int c=0; c < nchan; c++) printf(" %d",tchan[c]);
    printf("\n");
  }
};	//  Geom::prepareTH

double Geom::getNodes(int t, int z, double x2d[], double y2d[], int node[])
{
  int p0 = tri[6*t];
  int p1 = tri[6*t+1];
  int p2 = tri[6*t+2];
  int p3 = tri[6*t+3];
  int p4 = tri[6*t+4];
  int p5 = tri[6*t+5];
  /*  get position  */
  x2d[0] = px[p0];	y2d[0] = py[p0];
  x2d[1] = px[p1];      y2d[1] = py[p1];
  x2d[2] = px[p2];      y2d[2] = py[p2];
  /* get 3d point numbers  */
  node[0] = top[p0] + 2*z;  node[6] = node[0]+1;  node[3] = node[0]+2;
  node[1] = top[p1] + 2*z;  node[7] = node[1]+1;  node[4] = node[1]+2;
  node[2] = top[p2] + 2*z;  node[8] = node[2]+1;  node[5] = node[2]+2;
  node[9]  = top[p3] + z;    	node[12] = node[9]+1;
  node[10] = top[p4] + z;     	node[13] = node[10]+1;
  node[11] = top[p5] + z;	node[14] = node[11]+1;
  return dz[z];
};	// Geom::getNodes


void Geom::printMat3D()
{
  printf("#z:");
  double zh=0;
  for(int z=0; z < nz; z++) {
    zh = zh+dz[z];
    printf(" %.1lf",zh);
  }
  printf("\n");
  double areasum = 0;
  double toparea = 0;
  double fvol = 0;
  for(int t=0; t < ntri; t++) {
    double area = Area(t);
    areasum += area;
    if(chtri[t] >= 0) toparea += area;
    printf("t%2d:",t);
    for(int z=0; z < nz; z++) {
      char ch = ' ';
      if(ispowz[z]) ch = '*';
      if(chtri[t] < 0) ch = ' ';
      printf(" %d%c",mat3d[t*nz+z],ch);
      if(ch == '*') fvol += area*dz[z];
    }
    printf("\n");
  }
  printf("  areasum=%.3le toparea=%.3le fvol=%.3le\n",areasum,toparea,fvol);
};	// Geom::printMat3D

double Geom::Area(int t)
{
  int p0 = tri[6*t];	int p1 = tri[6*t+1];	int p2 = tri[6*t+2];
  double x0 = px[p0];	double x1 = px[p1];	double x2 = px[p2];
  double y0 = py[p0];   double y1 = py[p1];     double y2 = py[p2];
  double area = (x0*y1-x1*y0 + x1*y2-x2*y1 + x2*y0-x0*y2)/2;
  return area;
};	// Geom::Area

void Geom::findNeighbor()
{
  nbtri = new int[6*ntri];
  for(int ij=0; ij < 6*ntri; ij++) nbtri[ij] = -1;
  /*  for each triangle check  */
  for(int t1=0; t1 < ntri; t1++) {
    /* check edge 0 */
    int fnd = 0;
    int p1 = tri[6*t1];  int p2 = tri[6*t1+1];
    for(int t2=t1+1; !fnd && (t2 < ntri); t2++) {
      int p3 = tri[6*t2+1]; 
      if(p3 == p1) {
        int p4 = tri[6*t2+0];
        if(p4 == p2) {
          assert(nbtri[6*t1] < 0);
          nbtri[6*t1] = t2;
          nbtri[6*t1+1] = 0;
          fnd = 1;
          break;
        }
      }
      p3 = tri[6*t2+0];
      if(p3 == p1) {
        int p4 = tri[6*t2+2];
        if(p4 == p2) {
          assert(nbtri[6*t1] < 0);
          nbtri[6*t1] = t2;
          nbtri[6*t1+1] = 1;
          fnd = 1;
          break;
        }
      }
      p3 = tri[6*t2+2];
      if(p3 == p1) {
        int p4 = tri[6*t2+1];
        if(p4 == p2) {
          assert(nbtri[6*t1] < 0);
          nbtri[6*t1] = t2;
          nbtri[6*t1+1] = 2;
          fnd = 1;
          break;
        }
      }
    }	// for t2
    if(fnd) {	// symmetrize
      int ts = nbtri[6*t1];
      int js = nbtri[6*t1+1];
      assert(nbtri[6*ts+2*js] < 0);
      nbtri[6*ts+2*js] = t1;
      nbtri[6*ts+2*js+1] = 0;
    }
    /* check edge 1 */
    fnd = 0;
    p1 = tri[6*t1+1];  p2 = tri[6*t1+2];
    for(int t2=t1+1; !fnd && (t2 < ntri); t2++) {
      int p3 = tri[6*t2+1];
      if(p3 == p1) {
        int p4 = tri[6*t2+0];
        if(p4 == p2) {
          assert(nbtri[6*t1+2] < 0);
          nbtri[6*t1+2] = t2;
          nbtri[6*t1+3] = 0;
          fnd = 1;
          break;
        }
      }
      p3 = tri[6*t2+0];
      if(p3 == p1) {
        int p4 = tri[6*t2+2];
        if(p4 == p2) {
          assert(nbtri[6*t1+2] < 0);
          nbtri[6*t1+2] = t2;
          nbtri[6*t1+3] = 2;
          fnd = 1;
          break;
        }
      }
      p3 = tri[6*t2+2];
      if(p3 == p1) {
        int p4 = tri[6*t2+1];
        if(p4 == p2) {
          assert(nbtri[6*t1+2] < 0);
          nbtri[6*t1+2] = t2;
          nbtri[6*t1+3] = 1;
          fnd = 1;
          break;
        }
      }
    }   // for t2
    if(fnd) {   // symmetrize
      int ts = nbtri[6*t1+2];
      int js = nbtri[6*t1+3];
      assert(nbtri[6*ts+2*js] < 0);
      nbtri[6*ts+2*js] = t1;
      nbtri[6*ts+2*js+1] = 1;
    }
    /*  check edge 2  */
    fnd = 0;
    p1 = tri[6*t1+2];  p2 = tri[6*t1];
    for(int t2=t1+1; !fnd && (t2 < ntri); t2++) {
      int p3 = tri[6*t2+2];
      if(p3 == p1) {
        int p4 = tri[6*t2+1];
        if(p4 == p2) {
          nbtri[6*t1+4] = t2;
          nbtri[6*t1+5] = 1;
          fnd = 1;
          break;
        }
      }
      p3 = tri[6*t2+1];
      if(p3 == p1) {
        int p4 = tri[6*t2];
        if(p4 == p2) {
          nbtri[6*t1+4] = t2;
          nbtri[6*t1+5] = 0;
          fnd = 1;
          break;
        }
      }
      p3 = tri[6*t2+0];
      if(p3 == p1) {
        int p4 = tri[6*t2+2];
        if(p4 == p2) {
          nbtri[6*t1+4] = t2;
          nbtri[6*t1+5] = 2;
          fnd = 1;
          break;
        }
      }
    }   // for t2
    if(fnd) {   // symmetrize
      int ts = nbtri[6*t1+4];
      int js = nbtri[6*t1+5];
      assert(nbtri[6*ts+2*js] < 0);
      nbtri[6*ts+2*js] = t1;
      nbtri[6*ts+2*js+1] = 2;
    }
  }	// for t1
  /*  notify Dirichlet Boundary as -2  */
  for(int b=0; b < nbdry; b++) {
    int tb = btri[b];
    int eb = bedge[b];
    nbtri[6*tb+2*eb] = -2;
  };
  if(prlev > 1) {
    for(int t=0; t < ntri; t++) {
      printf("nbtri%2d:",t);
      for(int i=0; i < 6; i++) printf(" %d",nbtri[6*t+i]);
      printf("\n");
    }
  }
  int bad = 0;
  for(int t=0; t < ntri; t++) {
    /* check conformance  */
    for(int j=0; j < 3; j++) {
      int t2 = nbtri[6*t+2*j];
      int e2 = nbtri[6*t+2*j+1];
      if(t2 >= 0) {
        int t1 = nbtri[6*t2+2*e2];
        int e1 = nbtri[6*t2+2*e2+1];
        if((t1 != t) || (e1 != j)) bad++;
        assert(t1 == t);
        assert(e1 == j);
      }
    }
  }
  if(bad) prlev = 9;
  if(prlev) {
    for(int t=0; t < ntri; t++) {
      printf("nbtri%2d:",t);
      for(int i=0; i < 6; i++) printf(" %d",nbtri[6*t+i]);
      printf("\n");
    }
    if(prlev > 8) {
      printf("*** %s:%d Confomation error.\n",__FILE__,__LINE__);
      exit(0);
    }
  }
};	// Geom::findNeighbor

void Geom::makeFVM()
{
  FVM_area = new double[ntri];
  FVM_edge = new double[3*ntri];
  FVM_dist = new double[3*ntri];
  for(int t=0; t < ntri; t++) {
    int p1 = tri[6*t];  int p2 = tri[6*t+1];  int p3 = tri[6*t+2];
    double x1 = px[p1];	double x2 = px[p2];  double x3 = px[p3];
    double y1 = py[p1]; double y2 = py[p2];  double y3 = py[p3];
    double area = (x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;
    double x0 = (x1+x2+x3)/3;  double y0 = (y1+y2+y3)/3;
    double m12 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    double d12 =    -((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/m12;
    double m23 = sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));
    double d23 =    -((x3-x2)*(y2-y0)-(x2-x0)*(y3-y2))/m23;
    double m31 = sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3));
    double d31 =    -((x1-x3)*(y3-y0)-(x3-x0)*(y1-y3))/m31;
    if(prlev) {
      printf("t%d: area=%.3lf  m=%.2lf %.2lf %.2lf  d=%.2lf %.2lf %.2lf\n",
	t,area,m12,m23,m31, d12,d23,d31);
    }
    FVM_area[t] = area;
    FVM_edge[3*t] = m12;  FVM_edge[3*t+1] = m23;  FVM_edge[3*t+2] = m31;
    FVM_dist[3*t] = d12;  FVM_dist[3*t+1] = d23;  FVM_dist[3*t+2] = d31;
  }
};	// Geom::makeFVM

double Geom::getFVM(int t, int z, int tb[], double area[], double dist[])
{
  /* returns FVM quantities  */
  double vol = FVM_area[t]*dz[z];
  /*  upside  */
  if(z <= 0) tb[0] = -2;	// enforce Zero BC
  else tb[0] = t*nz+z-1;
  area[0] = FVM_area[t]; dist[0] = dz[z]/2;
  /*  downside */
  if(z >=  nz-1) tb[1] = -2;	// enforce Zero BC
  else tb[1] = t*nz+z+1;
  area[1] = FVM_area[t]; dist[1] = dz[z]/2;
  /*  3 side  */
  for(int e=0; e < 3; e++) {
    int nbt = nbtri[6*t+2*e];
    if(nbt < 0) tb[2+e] = nbt;
    else tb[2+e] = nbt*nz+z;
    area[2+e] = FVM_edge[3*t+e]*dz[z];
    dist[2+e] = FVM_dist[3*t+e];
  }
  return vol;
};	// Geom::getFVM

double Geom::getVol(int t, int z)
{
  return FVM_area[t]*dz[z];
};	// Geom::getVol

double Geom::getDistL(int etb, int tz)
{
  /*  compute edge distance from etb to tz  */
  int t0 = etb/nz;
  int z0 = etb%nz;
  int t1 = tz/nz;
  int z1 = tz%nz;
  if(t0 == t1) {
    assert(abs(z0-z1) == 1);
    return dz[z0]/2;
  }
  else if(z0 == z1) {
    for(int e=0; e < 3; e++) if(nbtri[6*t0+2*e]==t1) return FVM_dist[3*t0+e];
    printf("\n*** %s:%d look %d and from %d\n",__FILE__,__LINE__,t0,t1);
    printf("t%d: points",t1);
    for(int p=0; p < 3; p++) printf(" %d",tri[6*t1+p]);
    printf("\n");
    printf("  nbh:");
    for(int e=0; e < 3; e++) printf(" %d",nbtri[6*t1+2*e]);
    printf("\n");
    printf("t%d:",t0);
    for(int p=0; p < 3; p++) printf(" %d",tri[6*t0+p]);
    printf("\n");
    printf("  nbh:");
    for(int e=0; e < 3; e++) printf(" %d",nbtri[6*t0+2*e]);
    printf("\n");
    exit(0);
  }
  else {
    printf("   nz=%d\n",nz);
    printf("*** %s:%d impossible etb=%d tz=%d\n",__FILE__,__LINE__,etb,tz);
    printf("  t0=%d z0=%d  t1=%d z1=%d\n",t0,z0, t1,z1);
    exit(0);
  }
};	// Geom::getDistL

double Geom::getChanArea()
{
  /*  returns total channel area  */
  double totarea = 0;
  for(int c=0; c < nchan; c++) {
    int t = tchan[c];
    totarea += FVM_area[t];
  }
  return totarea;
}	// Geom::getChanArea

double Geom::getPowVol()
{
  /*  returns total volume of power producting elements */
  double pvol = 0.0;
  for(int c=0; c < nchan; c++) {
    if(ispp[c]) {
      double area = FVM_area[tchan[c]];
      for(int z=0; z < nz; z++) {
        if(ispowz[z]) pvol += area*dz[z];
      }	// for z
    }	 // if ispp
  }	// for c
  return pvol;
};	//  Geom::getPowVol

