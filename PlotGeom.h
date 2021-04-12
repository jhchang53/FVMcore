/*
 * 	PlotGeom.h
 */
extern "C" {
#include <cairo/cairo-svg.h>
};

#include "Geom.h"

class PlotGeom {
public:
  PlotGeom();
  ~PlotGeom();
  void plot(Geom *geom, const char pngfn[]);
private:
  void getBounds();
  void plotPoints();
  void plotMesh();
  void plotBdry();
  void plotChan();
  void plotTriNum();
  void plotCrod();
  void plotSave(const char pngfn[]);

  double TWOPI;
  Geom *geom;
  int npts,ntri;
  double *px,*py;
  int *tri;
  /*  cairo plot  */
  cairo_surface_t* surface;
  cairo_t* cr;
  double scale,offx,offy;

};
