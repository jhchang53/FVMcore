/*
 * 	MPupdate.h - derived class of MPwork
 *	update SegCdond and TriCond property
 *     Usage:
 * 	MPupdate(); - construtor
 *	~MPupdate(); - destructor
 * 	void init(int npows); - initialize for size npows vectors
 * 	  npows : size of vectors
 *	double updateProp(SegCond *seg, TriCond *tri, Vec fluvec, Vec fimavec);
 *	 - update conductivity/heat capacity for SegCond/TriCond model
 *	  returns difference in temperature
 *	  seg : SegCond class
 *	  tri : TriCond class
 *	  output:
 *	  fluvec : fluence vector
 *	  fimavec : % fima vector
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
