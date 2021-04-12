/*
 * 	GeomCAPP.cpp
 * 	CAPP type geometry
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "GeomCAPP.h"

GeomCAPP::GeomCAPP()
{
  prlev = 0;
  rt3 = sqrt(3.0);
  hex = NULL;
  map = NULL;
  hbdrytype = NULL;
  htri = NULL;
};

GeomCAPP::~GeomCAPP()
{
  delete [] hex;
  delete [] map;
  delete [] hbdrytype;
  delete [] htri;
};

void GeomCAPP::setPrlev(int prl)
{
  prlev = prl;
};	// GeomCAPP::setPrlev

void GeomCAPP::setPlane(double side, int nrow_i, int mat[])
{
  nrow = nrow_i;
  /* assign hex number  */
  hex = new int[nrow*nrow];
  for(int i=0; i < nrow*nrow; i++) hex[i] = -1;
  int j = 0;
  int jm = 0;
  for(int r=0; r < nrow; r++) {
    for(int c=0; c <= r; c++) {
      hex[j] = jm;
      int matj = mat[j];
      if(matj < 0) hex[j] = -1;
      else jm++;
      j++;
    }
  }
  nentry = j;
  if(prlev > 1) {
    printf(" nentry=%d:",nentry);
    for(int i=0; i < nentry; i++) printf(" %d",hex[i]);
    printf("\n");
  }
  genMap();
  genHbdry();
  /*  count number of triagles and points to be generated */
  int nt = 1;	// central one is a trinigle
  int np = 6;
  for(int nr=1; nr < nrow; nr++) {
    for(int nc=0; nc < nr; nc++) {
      int rc = nr*nrow+nc;
      if(hbdrytype[rc] > 0) {	// 0 is used obly for firstr triangle
        nt += 6;
        /* 19 points for full hexagon  */
        switch(hbdrytype[rc]) {
         case 1:  np += 16;  break;
         case 2:  np += 14;  break;
         case 3:  np += 12;  break;
         case 4:  np += 14;  break;
         default:
           printf("*** %s:%d unknown hex boundary type at %d = %d\n",
		__FILE__,__LINE__,rc,hbdrytype[rc]);
           exit(0);
        }	// switch
      }	// if hbdrytype
    }	// for nc
  }	// for nc
  ntri = nt;
  npts2d = np;
  if(prlev) printf(" %d triangles and %d point are expected.\n",nt,np);
  makePlane(side);
  setBdry();
  setMat2D(mat);

  findNeighbor();
  makeFVM();
};	// GeomCAPP::setPlane

void GeomCAPP::genMap()
{
  /*  use rectangular map */
  map = new int[nrow*nrow];
  for(int ij=0; ij < nrow*nrow; ij++) map[ij] = -1;
  /*  first row is a special triangle  */
  map[0] = hex[0];
  int nh = 1;
  for(int nr=1; nr < nrow; nr++) {
    for(int nc=0; nc < nr; nc++) {
      map[nr*nrow+nc] = hex[nh];
      nh++;
    }
  }
  if(prlev > 1) {
    printf("map:");
    for(int nr=0; nr < nrow; nr++) {
      for(int nc=0; nc < nrow; nc++) printf(" %d",map[nr*nrow+nc]);
      printf("\n");
    }
  }
  nhex = nh;
  if(prlev) printf(" we have %d non empty hexagons\n",nhex);
};	// GeomCAPP::genMap

void GeomCAPP::genHbdry()
{
  /*  generate boundary type  */
  hbdrytype =  new int[nrow*nrow];
  for(int ij=0; ij < nrow*nrow; ij++) hbdrytype[ij] = -1;
  /*  first row is a special triangle  */
  hbdrytype[0] = 0;
  for(int nr=1; nr < nrow; nr++) {
    for(int nc=0; nc < nr; nc++) {
      int rc = nr*nrow+nc;
      if(map[rc] >= 0) {
        if(nc==0) {
          hbdrytype[rc] = 1;
        }
        else if(nc == nr-1) {
          hbdrytype[rc] = 2;
        }
        else {
          if(map[rc-1] < 0) hbdrytype[rc] = 4;
          else  hbdrytype[rc] = 3;
        }
      }
    }
  }	// for nr
  if(prlev > 1) {
    printf("hbdrytype:");
    for(int nr=0; nr < nrow; nr++) {
      for(int nc=0; nc < nrow; nc++) printf(" %d",hbdrytype[nr*nrow+nc]);
      printf("\n");
    }
  }
};	// GeomCAPP::genHbdry

void GeomCAPP::makePlane(double s)
{
  double h = rt3/2*s;
  tri = new int[6*ntri];
  px = new double[npts2d];
  py = new double[npts2d];
  /*  set to nan */
  for(int ij=0; ij < npts2d; ij++) {
    px[ij] = 0.0/0.0;  py[ij] = 0.0/0.0;
  }
  htri = new int[6*nhex];	// triangle number of a hexagon
  for(int ij=0; ij < 6*nhex; ij++) htri[ij] = -1;
  /*  special case : first triangle  */
  assert(hbdrytype[0] == 0);
  int p = 0;
  px[0] =  0.0;	py[0] = 0.0;	p++;
  px[1] = -s/4;	py[1] = -h/2;	p++;
  px[2] =  s/4;	py[2] = -h/2;	p++;
  px[3] = -s/2;	py[3] = -h;	p++;
  px[4] =  0.0;	py[4] = -h;	p++;
  px[5] =  s/2;	py[5] = -h;	p++;
  tri[0] = 3; tri[1] = 5; tri[3] = 0; tri[3] = 4; tri[4] = 2; tri[5] = 1;
  htri[0] = 0;
  int t = 1;
  int nh = 1;
  int hexup,tup,hexleft,tleft,hexdn,tdn;
  int ptvar[4];
  ptvar[0] = nh;  ptvar[1] = p;  ptvar[2] = t;
  for(int nr=1; nr < nrow; nr++) {
    double cy = -2*h*nr-h;	// hexagon center
    for(int nc=0; nc < nr; nc++) {
      double cx = 1.5*s*nc;	
      cy += h;
      int rc = nr*nrow+nc;
      if(map[rc] > 0) {
        ptvar[0] = nh;
        switch(hbdrytype[rc]) {
        case 1: 
          hexup = map[(nr-1)*nrow+nc];
          tup = htri[6*hexup];
          makeType1(ptvar, tup, cx,cy,s,h); 
          break;	//  case 1
         case 2:
          hexleft = map[(nr-1)*nrow+nc-1];
          tleft = htri[6*hexleft+1];
          hexdn = map[nr*nrow+nc-1];
          tdn = htri[6*hexdn+2];
          makeType2(ptvar, tleft,tdn, cx,cy,s,h);
          break;	// case 2
         case 3:
          hexup = map[(nr-1)*nrow+nc];
          tup = htri[6*hexup];
          hexleft = map[(nr-1)*nrow+nc-1];
          tleft = htri[6*hexleft+1];
          hexdn = map[nr*nrow+nc-1];
          tdn = htri[6*hexdn+2];
          makeType3(ptvar, tup,tleft,tdn, cx,cy,s,h);
          break;	// case 3
         case 4:
          hexup = map[(nr-1)*nrow+nc];
          tup = htri[6*hexup];
          hexleft = map[(nr-1)*nrow+nc-1];
          tleft = htri[6*hexleft+1];
          makeType4(ptvar, tup,tleft, cx,cy,s,h);
          break;        // case 4
         default:
           printf("*** %s:%d hbdrytype %d is not known at nr=%d, nc=%d\n",
		__FILE__,__LINE__, hbdrytype[rc], nr,nc);
           exit(0);
        }	// for switch
        nh++;
      }	// if non empty hex
    }	// for nc
  }	// for nr
  assert(ptvar[1] ==  npts2d);
  assert(ptvar[2] == ntri);
  if(prlev > 1) {
    for(int t=0; t < ntri; t++) {
      printf("t%d:",t);
      for(int i=0; i < 6; i++) {
        int p =  tri[6*t+i];
        printf(" %d",p);
      }
      printf("\n");
    }
    printf("== npts2d=%d\n",npts2d);
    for(int t=0; t < ntri; t++) {
      printf("t%d:",t);
      for(int i=0; i < 6; i++) {
        int p =  tri[6*t+i];
        printf(" %d(%.2lf,%.2lf)",p,px[p],py[p]);
      }
      printf("\n");
    }
  }
};	// GeomCAPP::makePlane

void GeomCAPP::makeType1(int ptvar[], int tup,
	double cx, double cy, double s, double h)
{
  int nh = ptvar[0];
  int p = ptvar[1];
  int t = ptvar[2];
  /*				*/
  /*       x   x   x    	*/
  /*     0   1   2   3		*/
  /*   4   5   6   7   8	*/
  /*     9  10  11  12		*/
  /*      13  14  15		*/
  /*				*/
  int p0 = p; // 16 points will be defined
  px[p] = cx-3*s/4;       py[p] = cy+h/2;         p++;
  px[p] = cx  -s/4;       py[p] = cy+h/2;         p++;
  px[p] = cx  +s/4;       py[p] = cy+h/2;         p++;
  px[p] = cx+3*s/4;       py[p] = cy+h/2;         p++;
  px[p] = cx  -s;         py[p] = cy;             p++;
  px[p] = cx  -s/2;       py[p] = cy;             p++;
  px[p] = cx;             py[p] = cy;             p++;
  px[p] = cx  +s/2;       py[p] = cy;             p++;
  px[p] = cx  +s;         py[p] = cy;             p++;
  px[p] = cx-3*s/4;       py[p] = cy-h/2;         p++;
  px[p] = cx  -s/4;       py[p] = cy-h/2;         p++;
  px[p] = cx  +s/4;       py[p] = cy-h/2;         p++;
  px[p] = cx+3*s/4;       py[p] = cy-h/2;         p++;
  px[p] = cx  -s/2;       py[p] = cy-h;           p++;
  px[p] = cx;             py[p] = cy-h;           p++;
  px[p] = cx  +s/2;       py[p] = cy-h;           p++;
  assert(p-p0 == 16);
  /*  triangle 0 */
  int t0 = t;
  tri[6*t]  = p0+13;  tri[6*t+1]= p0+15;  tri[6*t+2]= p0+6;
  tri[6*t+3]= p0+14;  tri[6*t+4]= p0+11;  tri[6*t+5]= p0+10;
  htri[6*nh] = t;         t++;
  tri[6*t]  = p0+15;  tri[6*t+1]= p0+8;  tri[6*t+2]= p0+6;
  tri[6*t+3]= p0+12;  tri[6*t+4]= p0+7;  tri[6*t+5]= p0+11;
  htri[6*nh+1] = t;        t++;
  tri[6*t]  = p0+8;   tri[6*t+1]= tri[6*tup+1]; tri[6*t+2] = p0+6;
  tri[6*t+3]= p0+3;   tri[6*t+4]= p0+2;  tri[6*t+5]= p0+7;
  htri[6*nh+2] = t;     t++;
  tri[6*t] = tri[6*tup+1];  tri[6*t+1]=tri[6*tup]; tri[6*t+2]= p0+6;
  tri[6*t+3]= tri[6*tup+3]; tri[6*t+4]=p0+1; tri[6*t+5]=p0+2;
  htri[6*nh+3] = t;     t++;
  tri[6*t] = tri[6*tup];  tri[6*t+1]= p0+4;  tri[6*t+2]= p0+6;
  tri[6*t+3]= p0;    tri[6*t+4]= p0+5;  tri[6*t+5]= p0+1;
  htri[6*nh+4] = t;	t++;
  tri[6*t]  = p0+4;  tri[6*t+1]= p0+13;  tri[6*t+2]= p0+6;
  tri[t*6+3]= p0+9;  tri[6*t+4]= p0+10;  tri[6*t+5]= p0+5;
  htri[6*nh+5] = t;	t++;
  assert(t-t0 == 6);
  if(prlev > 1) {
    printf("htri(%d):",nh);
    for(int i=0; i < 6; i++) printf(" %d",htri[6*nh+i]);
    printf("\n");
  }
  ptvar[1] = p;
  ptvar[2] = t;
};	// makeType1

void GeomCAPP::makeType2(int ptvar[], int tleft, int tdn,
        double cx, double cy, double s, double h)
{
  int nh = ptvar[0];
  int p = ptvar[1];
  int t = ptvar[2];
  /*                           	*/
  /*       x   0   1        	*/
  /*     x   2   3   4     	*/
  /*   x   5   6   7   8    	*/
  /*     x   9  10  11        	*/
  /*       x  12  13        	*/
  /*				*/
  int p0 = p; // 14 points will be defined
  px[p] = cx;		py[p] = cy+h;		p++;
  px[p] = cx  +s/2;	py[p] = cy+h;		p++;
  px[p] = cx  -s/4;     py[p] = cy+h/2;         p++;
  px[p] = cx  +s/4;     py[p] = cy+h/2;         p++;
  px[p] = cx+3*s/4;     py[p] = cy+h/2;         p++;
  px[p] = cx  -s/2;     py[p] = cy;             p++;
  px[p] = cx;           py[p] = cy;             p++;
  px[p] = cx  +s/2;     py[p] = cy;             p++;
  px[p] = cx  +s;       py[p] = cy;             p++;
  px[p] = cx  -s/4;     py[p] = cy-h/2;         p++;
  px[p] = cx  +s/4;     py[p] = cy-h/2;         p++;
  px[p] = cx+3*s/4;     py[p] = cy-h/2;         p++;
  px[p] = cx;           py[p] = cy-h;           p++;
  px[p] = cx  +s/2;     py[p] = cy-h;           p++;
  assert(p-p0 == 14);
  /*  triangle 0 */
  int t0 = t;
  tri[6*t]  = tri[6*tdn];  tri[6*t+1]= p0+13; tri[6*t+2] = p0+6;
  tri[6*t+3]= p0+12;  tri[6*t+4]= p0+10;  tri[6*t+5]= p0+9;
  htri[6*nh] = t;         t++;
  tri[6*t]  = p0+13;  tri[6*t+1]= p0+8;  tri[6*t+2]= p0+6;
  tri[6*t+3]= p0+11;  tri[6*t+4]= p0+7;  tri[6*t+5]= p0+10;
  htri[6*nh+1] = t;        t++;
  tri[6*t]  = p0+8;   tri[6*t+1]= p0+1;  tri[6*t+2]= p0+6;
  tri[6*t+3]= p0+4;   tri[6*t+4]= p0+3;  tri[6*t+5]= p0+7;
  htri[6*nh+2] = t;     t++;
  tri[6*t]  = p0+1;   tri[6*t+1]=tri[6*tleft+1];     tri[6*t+2]= p0+6;
  tri[6*t+3]= p0;     tri[6*t+4]=p0+2;   tri[6*t+5]=p0+3;
  htri[6*nh+3] = t;     t++;
  tri[6*t] = tri[6*tleft+1];  tri[6*t+1]= tri[6*tleft];  tri[6*t+2]= p0+6;
  tri[6*t+3]= tri[6*tleft+3]; tri[6*t+4]= p0+5;  tri[6*t+5]= p0+2;
  htri[6*nh+4] = t;     t++;
  tri[6*t]  = tri[6*tleft];  tri[6*t+1]= tri[6*tdn];  tri[6*t+2]= p0+6;
  tri[t*6+3]= tri[6*tdn+3];  tri[6*t+4]= p0+9; tri[6*t+5]= p0+5;
  htri[6*nh+5] = t;     t++;
  assert(t-t0 == 6);
  if(prlev > 1) {
    printf("htri(%d):",nh);
    for(int i=0; i < 6; i++) printf(" %d",htri[6*nh+i]);
    printf("\n");
  }
  ptvar[1] = p;
  ptvar[2] = t;
};      // makeType2

void GeomCAPP::makeType3(int ptvar[], int tup, int tleft, int tdn,
        double cx, double cy, double s, double h)
{
  int nh = ptvar[0];
  int p = ptvar[1];
  int t = ptvar[2];
  /*                            */
  /*       x   x   x            */
  /*     x   0   1   2          */
  /*   x   3   4   5   6        */
  /*     x   7   8   9          */
  /*       x  10  11            */
  /*                            */
  int p0 = p; // 12 points will be defined
  px[p] = cx  -s/4;     py[p] = cy+h/2;         p++;
  px[p] = cx  +s/4;     py[p] = cy+h/2;         p++;
  px[p] = cx+3*s/4;     py[p] = cy+h/2;         p++;
  px[p] = cx  -s/2;     py[p] = cy;             p++;
  px[p] = cx;           py[p] = cy;             p++;
  px[p] = cx  +s/2;     py[p] = cy;             p++;
  px[p] = cx  +s;       py[p] = cy;             p++;
  px[p] = cx  -s/4;     py[p] = cy-h/2;         p++;
  px[p] = cx  +s/4;     py[p] = cy-h/2;         p++;
  px[p] = cx+3*s/4;     py[p] = cy-h/2;         p++;
  px[p] = cx;           py[p] = cy-h;           p++;
  px[p] = cx  +s/2;     py[p] = cy-h;           p++;
  assert(p-p0 == 12);
  /*  triangle 0 */
  int t0 = t;
  tri[6*t]  = tri[6*tdn];  tri[6*t+1]= p0+11; tri[6*t+2] = p0+4;
  tri[6*t+3]= p0+10;  tri[6*t+4]= p0+8;   tri[6*t+5]= p0+7;
  htri[6*nh] = t;         t++;
  tri[6*t]  = p0+11;  tri[6*t+1]= p0+6;  tri[6*t+2]= p0+4;
  tri[6*t+3]= p0+9;   tri[6*t+4]= p0+5;  tri[6*t+5]= p0+8;
  htri[6*nh+1] = t;        t++;
  tri[6*t]  = p0+6;   tri[6*t+1]= tri[6*tup+1];  tri[6*t+2]= p0+4;
  tri[6*t+3]= p0+2;   tri[6*t+4]= p0+1;  tri[6*t+5]= p0+5;
  htri[6*nh+2] = t;     t++;
  tri[6*t]  = tri[6*tup+1]; tri[6*t+1]=tri[6*tleft+1];  tri[6*t+2]= p0+4;
  tri[6*t+3]= tri[6*tup+3]; tri[6*t+4]=p0;   tri[6*t+5]=p0+1;
  htri[6*nh+3] = t;     t++;
  tri[6*t] = tri[6*tleft+1];  tri[6*t+1]= tri[6*tleft];  tri[6*t+2]= p0+4;
  tri[6*t+3]= tri[6*tleft+3]; tri[6*t+4]= p0+3;  tri[6*t+5]= p0;
  htri[6*nh+4] = t;     t++;
  tri[6*t]  = tri[6*tleft];  tri[6*t+1]= tri[6*tdn];  tri[6*t+2]= p0+4;
  tri[6*t+3]= tri[6*tdn+3];  tri[6*t+4]= p0+7;  tri[6*t+5]= p0+3;
  htri[6*nh+5] = t;     t++;
  assert(t-t0 == 6);
  if(prlev > 1) {
    printf("htri(%d):",nh);
    for(int i=0; i < 6; i++) printf(" %d",htri[6*nh+i]);
    printf("\n");
  }
  ptvar[1] = p;
  ptvar[2] = t;
};      // makeType3

void GeomCAPP::makeType4(int ptvar[], int tup, int tleft,
        double cx, double cy, double s, double h)
{
  int nh = ptvar[0];
  int p = ptvar[1];
  int t = ptvar[2];
  /*                            */
  /*       x   x   x            */
  /*     x   0   1   2          */
  /*   x   3   4   5   6        */
  /*     7   8   9  10        	*/
  /*      11  12  13            */
  /*                            */
  int p0 = p; // 14 points will be defined
  px[p] = cx  -s/4;     py[p] = cy+h/2;         p++;
  px[p] = cx  +s/4;     py[p] = cy+h/2;         p++;
  px[p] = cx+3*s/4;     py[p] = cy+h/2;         p++;
  px[p] = cx  -s/2;     py[p] = cy;             p++;
  px[p] = cx;           py[p] = cy;             p++;
  px[p] = cx  +s/2;     py[p] = cy;             p++;
  px[p] = cx  +s;       py[p] = cy;             p++;
  px[p] = cx-3*s/4;	py[p] = cy-h/2;		p++;
  px[p] = cx  -s/4;     py[p] = cy-h/2;         p++;
  px[p] = cx  +s/4;     py[p] = cy-h/2;         p++;
  px[p] = cx+3*s/4;     py[p] = cy-h/2;         p++;
  px[p] = cx  -s/2;	py[p] = cy-h;		p++;
  px[p] = cx;           py[p] = cy-h;           p++;
  px[p] = cx  +s/2;     py[p] = cy-h;           p++;
  assert(p-p0 == 14);
  /*  triangle 0 */
  int t0 = t;
  tri[6*t]  = p0+11;  tri[6*t+1]= p0+13; tri[6*t+2] = p0+4;
  tri[6*t+3]= p0+12;  tri[6*t+4]= p0+9;   tri[6*t+5]= p0+8;
  htri[6*nh] = t;         t++;
  tri[6*t]  = p0+13;  tri[6*t+1]= p0+6;  tri[6*t+2]= p0+4;
  tri[6*t+3]= p0+10;  tri[6*t+4]= p0+5;  tri[6*t+5]= p0+9;
  htri[6*nh+1] = t;        t++;
  tri[6*t]  = p0+6;   tri[6*t+1]= tri[6*tup+1];  tri[6*t+2]= p0+4;
  tri[6*t+3]= p0+2;   tri[6*t+4]= p0+1;  tri[6*t+5]= p0+5;
  htri[6*nh+2] = t;     t++;
  tri[6*t]  = tri[6*tup+1]; tri[6*t+1]=tri[6*tleft+1];  tri[6*t+2]= p0+4;
  tri[6*t+3]= tri[6*tup+3]; tri[6*t+4]=p0;   tri[6*t+5]=p0+1;
  htri[6*nh+3] = t;     t++;
  tri[6*t] = tri[6*tleft+1];  tri[6*t+1]= tri[6*tleft];  tri[6*t+2]= p0+4;
  tri[6*t+3]= tri[6*tleft+3]; tri[6*t+4]= p0+3;  tri[6*t+5]= p0;
  htri[6*nh+4] = t;     t++;
  tri[6*t]  = tri[6*tleft];  tri[6*t+1]= p0+11;  tri[6*t+2]= p0+4;
  tri[6*t+3]= p0+7;  tri[6*t+4]= p0+8;  tri[6*t+5]= p0+3;
  htri[6*nh+5] = t;     t++;
  assert(t-t0 == 6);
  if(prlev > 1) {
    printf("htri(%d):",nh);
    for(int i=0; i < 6; i++) printf(" %d",htri[6*nh+i]);
    printf("\n");
  }
  ptvar[1] = p;
  ptvar[2] = t;
};      // GeomCAPP::makeType4

void GeomCAPP::setBdry()
{
  /*  we assign  boundary using hbdrytype array  */
  int b = 0;
  for(int nc=0; nc < nrow; nc++) {
    int bn = 0;
    for(int nr=nc; nr < nrow; nr++) {
      int rc = nr*nrow+nc;
      if(map[rc] > 0) {
        switch(hbdrytype[rc]) {
         case 1: bn = 1;  break;
         case 2: case 3: bn = 2; break;
         case 4: bn = 3;  break;
         default:
           printf("*** %s:%d unknown hbdrytype %d\n",__FILE__,__LINE__,
		hbdrytype[rc]);
           exit(0);
        }
      }
    }	// for nr
    b += bn;
  }	// for nc
  nbdry = b;
  if(prlev) printf(" %d boundaries found.\n",nbdry);
  /*  now find boundary hex number */
  btri = new int[nbdry];
  bedge = new int[nbdry];
  b = 0;
  for(int nc=0; nc < nrow-1; nc++) {
    /*  we know that last column is not assigned  */
    int lastrc = -1;
    for(int nr=nc; nr < nrow; nr++) {
      int rc = nr*nrow+nc;
      if(map[rc] > 0) {
        lastrc = rc;
      }
    }   // for nr
    if(lastrc < 0) printf("lastrc=%d at nc=%d \n",lastrc,nc);
    assert(lastrc > 0);
    if(prlev > 2) printf("hex=%d type=%d\n",map[lastrc],hbdrytype[lastrc]);
    int hex = map[lastrc];
    switch(hbdrytype[lastrc]) {
     case 1:
       btri[b] = htri[6*hex]; 	bedge[b] = 0;	b++;
       break;
     case 2: case 3:
       btri[b] = htri[6*hex];	bedge[b] = 0;	b++;
       btri[b] = htri[6*hex+1]; bedge[b] = 0;	b++;
       break;
     case 4:
       btri[b] = htri[6*hex+5];	bedge[b] = 0;	b++;
       btri[b] = htri[6*hex];	bedge[b] = 0;	b++;
       btri[b] = htri[6*hex+1];	bedge[b] = 0;	b++;
       break;
     default:
	assert(1==0);
    }
  }	// for nc
  if(prlev > 1) {
    printf("btri:");
    for(b=0; b < nbdry; b++) printf(" %d",btri[b]);
    printf("\n");
  }
};	// GeomCAPP::setBdry

void GeomCAPP::setChannel(int ischan[])
{
  int c = 0;
  for(int n=0; n < nhex; n++) {
    for(int i=0; i < 6; i++) {
      if((ischan[n] > 0) && (htri[6*n+1] >= 0)) c++;
    }
  }
  if(prlev) printf(" nchan=%d\n",c);
  nchan = c;
  tchan = new int[nchan];
  ispp = new int[nchan];
  c = 0;
  int nh = 1;
  if(ischan[0] > 0) {
    tchan[c] = htri[0];
    c++;
  }
  /*  ischan follows input rule while htri follows non empty */
  int n = 1;
  for(int nr=1; nr < nrow; nr++) {
    for(int nc=0; nc < nr; nc++) {
      int rc = nr*nrow+nc;
      if(map[rc] > 0) {
        if(ischan[n] > 0) {
          for(int i=0; i < 6; i++) {
            if(htri[6*nh+i] >= 0) {
              tchan[c] = htri[6*nh+i];
              ispp[c] = 1;
              c++;
            }
          }	// for i
        }	// if ischan
        nh++;
      }
      n++;
    }
  }
  if(prlev > 1) {
    printf("tchan:");
    for(int i=0; i < c; i++) printf(" %d",tchan[i]);
    printf("\n");
  }
  assert(c == nchan);
};	// GeomCAPP::setChannel

void GeomCAPP::setCrod(int crod[])
{
  /*  set control rod number (start from 1) similar way to mat2d  */
  crod2d = new int[ntri];
  for(int t=0; t < ntri; t++) crod2d[t] = 0;
  crod2d[0] = crod[0];
  int n= 1;
  int nh = 1;
  for(int nr=1; nr < nrow; nr++) {
    for(int nc=0; nc < nr; nc++) {
      int rc = nr*nrow+nc;
      if(map[rc] > 0) {
        for(int i=0; i < 6; i++) {
          if(htri[6*nh+i] >= 0) {
            int t = htri[6*nh+i];
            crod2d[t] = crod[n];
          }     // for i
        }       // if ischan
        nh++;
      }
      n++;
    }
  }
  if(prlev > 1) {
    printf("crod2d:");
    for(int i=0; i < ntri; i++) printf(" %d",crod2d[i]);
    printf("\n");
    for(int i=0; i < ntri; i++) assert(crod2d[i] >= 0);
  }
};	// GeomCAPP::setCrod
  

void GeomCAPP::setMat2D(int mat[])
{
  /* set materials  */
  mat2d = new int[ntri];
  for(int t=0; t < ntri; t++) mat2d[t] = -1;
  mat2d[0] = mat[0];
  int n= 1;
  int nh = 1;
  for(int nr=1; nr < nrow; nr++) {
    for(int nc=0; nc < nr; nc++) {
      int rc = nr*nrow+nc;
      if(map[rc] > 0) {
        for(int i=0; i < 6; i++) {
          if(htri[6*nh+i] >= 0) {
            int t = htri[6*nh+i];
            mat2d[t] = mat[n];
          }     // for i
        }       // if ischan
        nh++;
      }
      n++;
    }
  }
  if(prlev > 1) {
    printf("mat2d:");
    for(int i=0; i < ntri; i++) printf(" %d",mat2d[i]);
    printf("\n");
    for(int i=0; i < ntri; i++) assert(mat2d[i] >= 0);
  }
};	//  GeomCAPP::setMat2D

void GeomCAPP::setHexMat(int h, int matz[])
{
  /*  set mat number */
  if((h < 0) || (h >= nhex)) {
    printf("*** %s:%d hexagon number %d is our of range. nhex=%d\n",
	__FILE__,__LINE__, h, nhex);
    exit(0);
  }
  /*  find triangle numbers  */
  for(int j=0; j < 6; j++) {
    int ht = htri[6*h+j];
    if(ht < 0) continue;
    int z = 0;
    for(int k=0; k < nzinp; k++) {
      for(int ki=0; ki < zseg[k]; ki++) {
        mat3d[ht*nz+z] = matz[k];
        z++;
      }
    }
    assert(z == nz);
  }
};	// GeomCAPP::setHexMat
#ifdef XX
void GeomCAPP::findNeighbor()
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

  for(int t=0; t < ntri; t++) {
    printf("nbtri%2d:",t);
    for(int i=0; i < 6; i++) printf(" %d",nbtri[6*t+i]);
    printf("\n");
  }
  int bad = 0;
  for(int t=0; t < ntri; t++) {
    printf("= check t%d:",t);
    /* check conformance  */
    for(int i=0; i < 3; i++) printf(" p%d",tri[6*t+i]);
    printf(" : ");
    for(int j=0; j < 3; j++) {
      int t2 = nbtri[6*t+2*j];
      int e2 = nbtri[6*t+2*j+1];
      printf(" t2=%d e2=%d ||",t2,e2);
      if(t2 >= 0) {
        for(int i=0; i < 3; i++) printf(" q%d",tri[6*t2+i]);
        printf(" |");
        int t1 = nbtri[6*t2+2*e2];
        int e1 = nbtri[6*t2+2*e2+1];
        printf("  t1=%d t=%d  e1=%d j=%d",t1,t, e1,j);
        if((t1 != t) || (e1 != j)) bad++;
        assert(t1 == t);
        assert(e1 == j);
      }
    }
    printf("\n");
  }
  if(bad) prlev = 9;
  if(prlev > 1) {
    for(int t=0; t < ntri; t++) {
      printf("nbtri%2d:",t);
      for(int i=0; i < 6; i++) printf(" %d",nbtri[6*t+i]);
      printf("\n");
    }
    if(prlev > 8) exit(0);
  }
};	// GeomCAPP::findNeighbor

void GeomCAPP::makeFVM()
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


    printf("t%d: area=%.3lf  m=%.2lf %.2lf %.2lf  d=%.2lf %.2lf %.2lf\n",
	t,area,m12,m23,m31, d12,d23,d31);
    FVM_area[t] = area;
    FVM_edge[3*t] = m12;  FVM_edge[3*t+1] = m23;  FVM_edge[3*t+2] = m31;
    FVM_dist[3*t] = d12;  FVM_dist[3*t+1] = d23;  FVM_dist[3*t+2] = d31;
  }
};
#endif
