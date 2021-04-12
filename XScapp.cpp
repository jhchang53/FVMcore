/*
 * 	XScapp.cpp
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "XScapp.h"

XScapp::XScapp()
{
  prlev = 0;
  db = NULL;
  tzs = NULL;
  zindbeg = NULL;
  sigbase = NULL;
  Tm_coef = NULL;
  Tf_coef = NULL;
  flu = NULL;
  fima = NULL;
};	// XScapp::XScapp

XScapp::~XScapp()
{
  sqlite3_close(db);
  delete [] tzs;
  delete [] zindbeg;
  delete [] sigbase;
  delete [] Tm_coef;
  delete [] Tf_coef;
  delete [] flu;
  delete [] fima;
};	// XScapp::~XScapp

void XScapp::setPrlev(int prl)
{
  prlev = prl;
};	// XScapp::setPrlev

void XScapp::setGeom(GeomCAPP *geom_i)
{
  geom = geom_i;
  ntri = geom->getNtri();
  nz = geom->getNz();
  if(prlev) printf("ntri*nz=%d\n",ntri*nz);
  ndata = ntri*nz;
  /*  map translate  */
  nrow = geom->getNrow();
  map = geom->Map();
  if(prlev) {
    for(int m=0; m < nrow; m++) {
      printf("map%d:",m);
      for(int i=0; i < nrow; i++) printf(" %d",map[m*nrow+i]);
      printf("\n");
    }
  }
  /*  nz translate  */
  nzinp = geom->getNzinp();
  int *zseg = geom->Zseg();
  zindbeg = new int[nzinp+1];
  zindbeg[0] = 0;
  for(int zs=0; zs < nzinp; zs++) zindbeg[zs+1] = zindbeg[zs] + zseg[zs];
  if(prlev) {
    printf("zindbeg:");
    for(int zs=0; zs < nzinp+1; zs++) printf(" %d",zindbeg[zs]);
    printf("\n");
  }
};	// XScapp::setGeom

void XScapp::open(const char *dbfn, int efpd_i, double Tm_ref_i, double Tf_ref_i)
{
  efpd = efpd_i;
  Tm_ref = Tm_ref_i;
  Tf_ref = Tf_ref_i;
  rt_ref = sqrt(Tf_ref);
  /*  check if db file exists really. It is required    */
  /*  because sqlite3_open may create an empty file     */
  struct stat stbuf;
  stat(dbfn,&stbuf);
  if(!S_ISREG(stbuf.st_mode) || (stbuf.st_size < 1000)) {
    printf("*** DB file %s does not exist.\n",dbfn);
    exit(0);
  }
  /*  open data base  */
  int rc = sqlite3_open(dbfn,&db);
  if(rc) {
    printf("*** Cannot open DB %s\n",sqlite3_errmsg(db));
    exit(0);
  }
  /*  get number of idx  */
  nids = getNumIds();
  setCappIndex();
  /*  we are dealing with 1G CS : sigtr, sigab, signf */
  nwords = 3;
  /*  get base cross sections  */
  getXSbase(Tm_ref,Tf_ref);
  /*  get Tm coef  */
  getTcoefs(Tm_ref,Tf_ref,1);
  getTcoefs(Tm_ref,Tf_ref,2);
  /*  get flu and fima  */
  getFluence();
};	// XScapp::open

int XScapp::getNumIds()
{
  /* count number of ids in SegIndex  */
  /* TABLE SegIndex (id INT PRIMARY KEY, i INT, j INT, k INT, h INT, z INT)  */
  char zSql[100];
  sprintf(zSql,"SELECT count(*) FROM SegIndex");
  sqlite3_stmt *stmt;
  int rc = sqlite3_prepare_v2(db,zSql,-1,&stmt,NULL);
  if(rc != SQLITE_OK) {
    printf("*** %s:%d rc=%d\n",__FILE__,__LINE__,rc);
    printf("*** %s\n",zSql);
    exit(0);
  }
  /* rc can be either SQLITE_ROW (100)  or SQLITE_DONE (101) */
  rc = sqlite3_step(stmt);
  if(rc != SQLITE_ROW) {
    printf("*** %s:%d rc=%d\n",__FILE__,__LINE__,rc);
    exit(0);
  }
  int nids = sqlite3_column_int(stmt,0);
  if(prlev) printf(" nids=%d\n",nids);
  sqlite3_finalize(stmt);
  return nids;
};      //  XScapp::getIndex

void XScapp::setCappIndex()
{
  /* TABLE SegIndex (id INT PRIMARY KEY, i INT, j INT, k INT, h INT, z INT)  */

  char zSql[100];
  sprintf(zSql,"SELECT id,i,j,k,h,z FROM SegIndex\n");
  sqlite3_stmt *stmt;
  int rc = sqlite3_prepare_v2(db,zSql,-1,&stmt,NULL);
  if(rc != SQLITE_OK) {
    printf("* %s:%d rc=%d\n",__FILE__,__LINE__,rc);
    exit(0);
  }
  /* rc can be either SQLITE_ROW (100)  or SQLITE_DONE (101) */
  tzs = new int[nids];
  for(int n=0; n < nids; n++) tzs[n] = -1;
  for(int n=0; n < nids; n++) {
    rc = sqlite3_step(stmt);
    if(rc != SQLITE_ROW) {
      printf("*** %s:%d rc=%d\n",__FILE__,__LINE__,rc);
      exit(0);
    }
    /*  translate i,j,k,h,z to tz  */
    int id = sqlite3_column_int(stmt,0);
    int ii = sqlite3_column_int(stmt,1);
    int jj = sqlite3_column_int(stmt,2);
    int kk = sqlite3_column_int(stmt,3);
    int hh = sqlite3_column_int(stmt,4);
    int zz = sqlite3_column_int(stmt,5);
    /* n-th entry has id */
    tzs[id] = indplane(ii,jj,hh)*nz + indaxial(kk,zz);
  }
  sqlite3_finalize(stmt);
  /*  check if all nodes are assigned  */
  for(int n=0; n < nids; n++) assert(tzs[n] >= 0);
};      //  XScapp::setCappIndex

int XScapp::indplane(int i, int j, int h)
{
  /*  translate plane to triangle number  */
  /*  special case to i=0 */
  if(i==0) {
    assert(j==0);
    assert(h == 0);
    return 0;
  }
  int hexno = map[i*nrow+j];
  int t = (hexno-1)*6 + h+1;
  if(prlev > 1) printf(" i=%d j=%d h=%d t=%d\n",i,j,h,t);
  if((t < 0) || (t >= ntri)) {
    printf("*** problem in i=%d j=%d h=%d, t=%d\n",i,j,h,t);
    exit(0);
  }
  return t;
};	// XScapp::indplane

int XScapp::indaxial(int k, int zz)
{
  /*  translate axial to axial number  */
  /*  we use zseg for axial translation  */
  int indz = zindbeg[k] + zz;
  return indz;
};	// XScapp::indaxial

/* TABLE XSdata (id INT, efpd INT, kind INT, Tmod FLOAT, Tfuel FLOAT, xsblob BLOB) */

void XScapp::getXSbase(double Tm, double Tf)
{
  /*  retrieve type 0 base XS with given Tm and Tf */
  sigbase = new double[nids*nwords];
  double *sigblob = new double[nwords];
  size_t sig_bytes = nwords*sizeof(double);
  char zSql[100];
  sprintf(zSql,"SELECT id,xsblob FROM XSdata WHERE efpd=%d AND kind=%d AND Tmod=%.0lf AND Tfuel=%.0lf",
        efpd,0, Tm,Tf);      // kind = 1 for base
  sqlite3_stmt *stmt;
  int rc = sqlite3_prepare_v2(db,zSql,-1,&stmt,NULL);
  if(rc != SQLITE_OK) {
    printf("*** %s:%d rc=%d\n",__FILE__,__LINE__,rc);
    printf("*** %s\n",zSql);
    exit(0);
  }
  for(int n=0; n < nids; n++) {
    rc = sqlite3_step(stmt);
    if(rc != SQLITE_ROW) {
      printf("*** %s:%d rc=%d\n",__FILE__,__LINE__,rc);
      exit(0);
    }
    int id = sqlite3_column_int(stmt,0);
    /*  xs is packed in BLOB  */
    const void* blob = sqlite3_column_blob(stmt,1);
    size_t blob_bytes = sqlite3_column_bytes(stmt,1);
    assert(blob_bytes == sig_bytes);
    memcpy(sigblob,blob,blob_bytes);
    /*  save to sig array  */
    int tz = tzs[id];
    for(int k=0; k < nwords; k++) sigbase[tz*nwords+k] = sigblob[k];
  }
  sqlite3_finalize(stmt);
  delete [] sigblob;
  if(prlev > 1) {
    for(int n=0; n < 20; n++) {
      printf("%d:",n);
      for(int k=0; k < nwords; k++) printf(" %.5le",sigbase[n*nwords+k]);
      printf("\n");
    }
  }
};      // XScapp::getXSbase

void XScapp::getTcoefs(double Tm, double Tf, int kind)
{
  /*  retrieve type 1 or 2 base XS with given Tm and Tf */
  double *coef;
  if(kind  == 1) {
    Tm_coef = new double[nids*nwords*2];        // 2 words for temp. variation
    coef = Tm_coef;
  }
  else if(kind == 2) {
    Tf_coef = new double[nids*nwords*2];
    coef = Tf_coef;
  }
  else {
    assert(1==0);
  }
  double *coefblob = new double[nwords*2];
  size_t coef_bytes = 2*nwords*sizeof(double);
  char zSql[100];
  sprintf(zSql,"SELECT id,xsblob FROM XSdata WHERE efpd=%d AND kind=%d AND Tmod=%.0lf AND Tfuel=%.0lf",
        efpd,kind, Tm,Tf);   // kind = 1/2 for Tmod/Tfuel coefficients
  sqlite3_stmt *stmt;
  int rc = sqlite3_prepare_v2(db,zSql,-1,&stmt,NULL);
  if(rc != SQLITE_OK) {
    printf("*** %s:%d rc=%d\n",__FILE__,__LINE__,rc);
    printf("*** %s\n",zSql);
    exit(0);
  }
  for(int n=0; n < nids; n++) {
    rc = sqlite3_step(stmt);
    if(rc != SQLITE_ROW) {
      printf("*** %s:%d at n=%d rc=%d\n",__FILE__,__LINE__,n,rc);
      exit(0);
    }
    int id = sqlite3_column_int(stmt,0);
    int tz = tzs[id];
    /*  xs is packed in BLOB  */
    const void* blob = sqlite3_column_blob(stmt,1);
    size_t blob_bytes = sqlite3_column_bytes(stmt,1);
    assert(blob_bytes == coef_bytes);
    memcpy(coefblob,blob,blob_bytes);
    for(int k=0; k < 2*nwords; k++) coef[2*tz*nwords+k] = coefblob[k];
  }
  sqlite3_finalize(stmt);
  delete [] coefblob;
  if(prlev > 1) {
    for(int n=0; n < 20; n++) {
      printf("%d:",n);
      for(int k=0; k < 2*nwords; k++) printf(" %.5le",coef[n*2*nwords+k]);
      printf("\n");
    }
  }
};      // XScapp::getTcoefs

void XScapp::getSig(int tz, double Tm, double Tf, double sig[])
{
  /* we know that sig is arranged in triangluar prism coord. tz  */
  /*  returns xs 3-ple  */
  double xm = Tm - Tm_ref;
  double xf = sqrt(Tf) - rt_ref;        // rt_ref = sqrt(Tf_ref)
  for(int k=0; k < 3; k++) {
    int nwk = tz*nwords+k;
    sig[k] = sigbase[nwk] \
        + xm*(Tm_coef[2*nwk] + xm*Tm_coef[2*nwk+1]) \
        + xf*(Tf_coef[2*nwk] + xf*Tf_coef[2*nwk+1]);
  }
  if(tz==-120) {
    printf("tz=%d:",tz);
    for(int k=0; k < 3; k++) {
      int nwk = tz*nwords+k;
      printf(" %.5le :",sigbase[nwk]);
      double x3 = 300.0-Tm_ref;
      double sig3 = sigbase[nwk] + x3*(Tm_coef[2*nwk]+x3*Tm_coef[2*nwk+1]);
      double x19 = 1900.0-Tm_ref;
      double sig19 = sigbase[nwk] + x19*(Tm_coef[2*nwk]+x19*Tm_coef[2*nwk+1]);
      printf(" %.5le %.5le\n",sig3,sig19);
    }
  }
};	// XScapp::getSig

double XScapp::getSigtr(int tz, double Tm, double Tf)
{
  double xm = Tm - Tm_ref;
  double xf = sqrt(Tf) - rt_ref;        // rt_ref = sqrt(Tf_ref)
  int k=0;	// sigtr is the first
  int nwk = tz*nwords+k;
  double str = sigbase[nwk] \
        + xm*(Tm_coef[2*nwk] + xm*(Tm_coef[2*nwk+1])) \
        + xf*(Tf_coef[2*nwk] + xf*(Tf_coef[2*nwk+1]));
  return str;
};	//  XScapp::getSigtr

double XScapp::getSigab(int tz, double Tm, double Tf)
{
  double xm = Tm - Tm_ref;
  double xf = sqrt(Tf) - rt_ref;        // rt_ref = sqrt(Tf_ref)
  int k=1;      // sigab is the second
  int nwk = tz*nwords+k;
  double sab = sigbase[nwk] \
        + xm*(Tm_coef[2*nwk] + xm*(Tm_coef[2*nwk+1])) \
        + xf*(Tf_coef[2*nwk] + xf*(Tf_coef[2*nwk+1]));
  return sab;
};	// XScapp::getSigab

double XScapp::getSignf(int tz, double Tm, double Tf)
{
  double xm = Tm - Tm_ref;
  double xf = sqrt(Tf) - rt_ref;        // rt_ref = sqrt(Tf_ref)
  int k=2;      // sigab is the second
  int nwk = tz*nwords+k;
  double snf = sigbase[nwk] \
        + xm*(Tm_coef[2*nwk] + xm*(Tm_coef[2*nwk+1])) \
        + xf*(Tf_coef[2*nwk] + xf*(Tf_coef[2*nwk+1]));
  return snf;
};	// XScapp::getSignf

/*  TABLE Fluence (id INT, efpd INT, flu FLOAT, fima FLOAT, UNIQUE(id,efpd)) */
void XScapp::getFluence()
{
  /*  retrieve fluence and fima  */
  flu = new double[nids];
  fima = new double[nids];
  char zSql[100];
  sprintf(zSql,"SELECT id,flu,fima FROM Fluence WHERE efpd=%d",efpd);
  sqlite3_stmt *stmt;
  int rc = sqlite3_prepare_v2(db,zSql,-1,&stmt,NULL);
  if(rc != SQLITE_OK) {
    printf("*** %s:%d rc=%d\n",__FILE__,__LINE__,rc);
    printf("    %s\n",zSql);
    printf("*** database may not opened.\n");
    exit(0);
  }
  for(int n=0; n < nids; n++) {
    rc = sqlite3_step(stmt);
    if(rc != SQLITE_ROW) {
      printf("*** %s:%d rc=%d\n",__FILE__,__LINE__,rc);
      exit(0);
    }
    int id = sqlite3_column_int(stmt,0);
    int tz = tzs[id];
    flu[tz] = sqlite3_column_double(stmt,1);
    fima[tz] = sqlite3_column_double(stmt,2);
  }
  sqlite3_finalize(stmt);
};	// XScapp::getFluence

double XScapp::getFlu(int tz)
{
  return flu[tz];
};	// XScapp::getFlu

double XScapp::getFima(int tz)
{
  return fima[tz];
};      // XScapp::getFima

