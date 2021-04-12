/*
 * 	GeomEtri.h
 * 	large equilateral triangles
 */
#include "Geom.h"

class GeomEtri : public Geom
{
public:
  GeomEtri();
  ~GeomEtri();
  void setBdry(int isbdry);
  int generate(double side, int nlev);
  void setChan(int nchin);
private:
  void genPoints();
  void genTris();
  void setBdry();
  void setMat2D();
  void setCrod();

  int prlev;
  int nlev,isbdry;
  double side;
  double rt3;
};

