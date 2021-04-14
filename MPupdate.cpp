/*
 * 	MPupdate.cpp
 *	update SegCdond and TriCond property
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "MPupdate.h"

MPupdate::MPupdate()
{
  prlev = 0;
  npows = 0;
};

MPupdate::~MPupdate()
{
};

void MPupdate::init(int npows_i)
{
  npows = npows_i;
  ierr = VecCreate(PETSC_COMM_WORLD,&Tcompvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(Tcompvec,PETSC_DECIDE,npows);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(Tcompvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(Tcompvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSet(Tcompvec,0.0);	CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecDuplicate(Tcompvec,&Terrvec); CHKERRABORT(PETSC_COMM_SELF,ierr);

  for(int m=0; m < 3; m++) {
    VecDuplicate(Tcompvec,&Tempvec[m]);
    VecDuplicate(Tcompvec,&condSvec[m]);
    VecDuplicate(Tcompvec,&rhocpSvec[m]);
  }
  for(int m=0; m < 6; m++) {
    VecDuplicate(Tcompvec,&condTvec[m]);
    VecDuplicate(Tcompvec,&rhocpTvec[m]);
  }

};	// MPupdate::init

double MPupdate::updateProp(SegCond *seg, TriCond *triso, Vec fluvec, Vec fimavec)
{
  assert(npows > 0);
  Vec Tcomp = seg->getTpinVec();
  Vec Tgra  = seg->getTgraVec();
  /*  update thermal properties */
  ierr = VecCopy(Tgra,Tempvec[0]);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecCopy(Tcomp,Tempvec[2]);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecWAXPY(Tempvec[1],1.0,Tgra,Tcomp);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScale(Tempvec[1],0.5);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  Vec radvec = seg->getRadVec();
  makeSegProp(fluvec,fimavec,3,Tempvec,radvec,condSvec,rhocpSvec);
  seg->setCond(3,condSvec,rhocpSvec);
  ierr = VecCopy(Tcomp,Terrvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAXPY(Terrvec,-1.0,Tcompvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("Tcomp:\n");
    VecView(Tcomp,PETSC_VIEWER_STDOUT_WORLD);
    printf("Tcompvec:\n");
    VecView(Tcompvec,PETSC_VIEWER_STDOUT_WORLD);
    printf("Terrvec:\n");
    VecView(Terrvec,PETSC_VIEWER_STDOUT_WORLD);
  }
  double errsum;
  VecNorm(Terrvec,NORM_1,&errsum);
  double errseg = errsum/npows;
  ierr = VecCopy(Tcomp,Tcompvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  makeTriProp(fluvec,fimavec,Tcomp, 6,condTvec,rhocpTvec);
  if(triso) triso->setCond(6,condTvec,rhocpTvec);
  return errseg;
};	// MPupdate::updateProp
