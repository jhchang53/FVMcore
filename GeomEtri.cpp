/*
 * 	GeomEtri.cpp
 *	generate equilateral triangles
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "GeomEtri.h"

GeomEtri::GeomEtri()
{
  prlev = 0;
  rt3 = sqrt(3.0);
  isbdry = 1;	// default to zero-incoming BC
};

GeomEtri::~GeomEtri()
{
};

void GeomEtri::setBdry(int isbdry_i)
{
  /*  set BD condition */
  isbdry = isbdry_i;
};

int GeomEtri::generate(double side_i, int nlev_i)
{
  side = side_i;
  nlev = nlev_i;
  if(prlev) printf("= Etri side=%.1lf nlev=%d\n",side,nlev);
  genPoints();
  genTris();
  if(isbdry) setBdry();
  setMat2D();
  setChan(nlev);
  setCrod();
  /*  for FVM  */
  findNeighbor();
  makeFVM();
  return ntri;
};	// GeomEtri::generate

void GeomEtri::genPoints()
{
  /*  generate points  */
  npts2d = (2*nlev+1)*(nlev+1);
  if(prlev) printf("nlev=%d nprs2d=%d\n",nlev,npts2d);
  px = new double[npts2d];
  py = new double[npts2d];
  int p =0;
  for(int n=0; n < 2*nlev+1; n++) {
    double x = n*side;	double y = 0.0;
    for(int i=0; i < n+1; i++) {
      px[p] = x;	py[p] = y;
      x -= 0.5*side;	y += 0.5*rt3*side;
      p++;
    }
  }
  if(prlev) {
    printf(" p=%d\n",p);
    for(int i=0; i < p; i++) {
      printf("P%d: %.2lf %.2le\n",i,px[i],py[i]);
    }
  }
  assert(p == npts2d);
};	// GeomEtri::genPoints

void GeomEtri::genTris()
{
  /*  generate triangles */
  ntri = nlev*nlev;
  tri = new int[6*ntri];
  int t = 0;
  for(int n=0; n < nlev; n++) {
    if(prlev > 1) printf("n=%d t=%d\n",n,t);
    int p0 = 2*n*n + n;
    int p1 = 2*n*n + 3*n + 1;
    int p2 = 2*n*n + 5*n + 3;
    if(prlev > 1) printf(" n=%d p0=%d p1=%d p2=%d\n",n,p0,p1,p2);
    for(int j=0; j <= 2*n ; j++) {
      if(j%2==0) {
        /* upside triangle  */
        tri[6*t]  =p0+j;  tri[6*t+1]=p2+j;   tri[6*t+2]=p2+j+2;
        tri[6*t+3]=p1+j;  tri[6*t+4]=p2+j+1; tri[6*t+5]=p1+j+1;
      }
      else {
        /*  downside triangle */
        tri[6*t]  =p0+j-1; tri[6*t+1]=p2+j+1; tri[6*t+2]=p0+j+1;
        tri[6*t+3]=p1+j;   tri[6*t+4]=p1+j+1; tri[6*t+5]=p0+j;
      }
      t++;
    }
  }
  if(prlev) {
    printf("ntri=%d t=%d\n",ntri,t);
    for(int k=0;  k < ntri; k++) {
      printf("t%2d:",k);
      for(int i=0; i < 6; i++) printf(" %3d",tri[6*k+i]);
      printf("\n");
    }
  }
};	// GeomEtri::genTris

void GeomEtri::setBdry()
{
  /*  set boundary triangle and edge  */
  nbdry = nlev;
  btri = new int[nbdry];
  bedge = new int[nbdry];
  for(int b=0; b < nbdry; b++) {
    int t = (nlev-1)*(nlev-1) + 2*b;
    btri[b] = t;
    bedge[b] = 1;
  }
  if(prlev) {
    printf("Btri:");
    for(int b=0; b < nbdry; b++) printf(" %2d(%d)",btri[b],bedge[b]);
    printf("\n");
  }
};	//  GeomEtri::setBdry

void GeomEtri::setChan(int nchin)
{
  /*  set channel assignment for triangle in nchin */
  /*  nchan : no. of channels  */
  /*  tchan : triangle no. of channel */
  /*  ispp : power production flag */
  assert(nchin <= nlev);
  nchan = nchin*nchin;
  tchan = new int[nchan];
  ispp = new int[nchan];
  int c = 0;
  int t = 0;
  for(int n=0; n < nchin*nchin; n++) {
    tchan[c] = t;  ispp[c] = 1;
    t++;   c++;
  }
};	// GeomEtri::setChan

void GeomEtri::setMat2D()
{
  mat2d = new int[ntri];
  for(int t=0; t < ntri; t++) mat2d[t] = 0;
};	// GeomEtri::setMat2D

void GeomEtri::setCrod()
{
  /*  control rod at all triangle  */
  /*  set control rod number (start from 1)  */
  crod2d = new int[ntri];
  for(int t=0; t < ntri; t++) crod2d[t] = 1;
};	// GeomEtri::setCrod

