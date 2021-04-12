/*
 * 	MPupdate.h
 *	update SegCdond and TriCond property
 *     Usage:
  MPupdate();
  ~MPupdate();
  void init(int npows);
  double updateProp(SegCond *seg, TriCond *tri, Vec fluvec, Vec fimavec);

 */
#include "MPwork.h"
#include "SegCond.h"
#include "TriCond.h"

class MPupdate : public MPwork {
public:
  MPupdate();
  ~MPupdate();
  void init(int npows);
  double updateProp(SegCond *seg, TriCond *tri, Vec fluvec, Vec fimavec);
private:
  int prlev;
  int npows;
  PetscErrorCode ierr;
  Vec Tempvec[3],condSvec[3],rhocpSvec[3];
  Vec condTvec[6],rhocpTvec[6];
  Vec Tcompvec,Terrvec;
};
