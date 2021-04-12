/*
 * 	PlotGeom.cpp
 */
#include <stdio.h>
#include <math.h>
#include "PlotGeom.h"

PlotGeom::PlotGeom()
{
  TWOPI = 4.0*acos(0.0);
  px = NULL;
  py = NULL;
  cr = NULL;
  surface = NULL;
};

PlotGeom::~PlotGeom()
{
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
};

void PlotGeom::plot(Geom *geom_i, const char pngfn[])
{
  geom = geom_i;
  npts = geom->getNpts2d();
  px = geom->Px();
  py = geom->Py();
  getBounds();
  plotPoints();
  ntri = geom->getNtri();
  tri = geom->Tri();
  plotMesh();
  plotBdry();
  plotChan();
  plotTriNum();
  if(geom->Crod2D()) plotCrod();
  plotSave(pngfn);
};	//  PlotGeom::plot

void PlotGeom::getBounds()
{
  double xmin = 10.0;   double xmax = -10.0;
  double ymin = 10.0;   double ymax = -10.0;
  for(int p=0; p < npts; p++) { 
    if(px[p] < xmin) xmin = px[p];
    if(px[p] > xmax) xmax = px[p];
    if(py[p] < ymin) ymin = py[p];
    if(py[p] > ymax) ymax = py[p];
  }
  /*  make picture frame  */
  int width = 2000;
  int heit = 2000;
  int offset = 100;
  double windx = xmax-xmin;
  double windy = ymax-ymin;
  if(windx > windy) {
    scale = (width-2*offset)/windx;
    heit = scale*windy + 2*offset;
  }
  else {
    scale = (heit-2*offset)/windy;
    width = scale*windx + 2*offset;
  }
  printf("width=%d heit=%d\n",width,heit);
  printf(" ymin=%.3lf ymax=%.3lf xmin=%.3lf xmax=%.3lf\n",
        ymin,ymax,xmin,xmax);
  offx = offset-scale*xmin;
  offy = heit - (offset-scale*ymin);

  surface =  cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
        width,heit);

  cr = cairo_create(surface);
  /*    white back ground */
  cairo_save(cr);
  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_paint(cr);
  cairo_restore(cr);
};	// PlotGeom::getBounds

void PlotGeom::plotPoints()
{
  cairo_set_source_rgba(cr, 1,0.2,0.2,0.6);
  
  for(int p=0; p < npts; p++) {
    cairo_arc(cr, offx+scale*px[p],offy-scale*py[p], 10.0, 0.0,TWOPI);
    cairo_fill(cr);
    cairo_stroke(cr);
  }
};	// PlotGeom::plotPoints

void PlotGeom::plotMesh()
{
  printf("offx=%.2le offy=%.2le scale=%.2le\n",offx,offy,scale);
  cairo_set_line_width(cr,2.0);
  cairo_set_source_rgb(cr, 0.0,0.0,0.0);
  /* plot order 2 triangles  */
  for(int t2=0; t2 < ntri; t2++) {
    double x0 = px[tri[6*t2]];    double y0 = py[tri[6*t2]];
    double x1 = px[tri[6*t2+1]];  double y1 = py[tri[6*t2+1]];
    double x2 = px[tri[6*t2+2]];  double y2 = py[tri[6*t2+2]];
    double x3 = px[tri[6*t2+3]];  double y3 = py[tri[6*t2+3]];
    double x4 = px[tri[6*t2+4]];  double y4 = py[tri[6*t2+4]];
    double x5 = px[tri[6*t2+5]];  double y5 = py[tri[6*t2+5]];

    cairo_move_to(cr,offx+scale*x0,offy-scale*y0);
    cairo_line_to(cr,offx+scale*x3,offy-scale*y3);
    cairo_line_to(cr,offx+scale*x1,offy-scale*y1);
    cairo_line_to(cr,offx+scale*x4,offy-scale*y4);
    cairo_line_to(cr,offx+scale*x2,offy-scale*y2);
    cairo_line_to(cr,offx+scale*x5,offy-scale*y5);
    cairo_line_to(cr,offx+scale*x0,offy-scale*y0);
    cairo_stroke(cr);
  }
  printf("== plotMsh with ntri=%d\n",ntri);
};      //  PlotGeom::plotMesh

void PlotGeom::plotBdry()
{
  /* plot boundary lines */
  int nbdry = geom->getNbdry();
  printf("= ndbry=%d\n",nbdry);
  int *btri = geom->Btri();
  int *bedge = geom->Bedge();
  cairo_set_source_rgb(cr, 1.0,0.0,0.0);
  cairo_set_line_width(cr,5.0);
  for(int b=0; b < nbdry; b++) {
    int t = btri[b];
    int e = bedge[b];
    int p0 = tri[6*t+e];
    int p1 = tri[6*t+(e+1)%3];
    printf(" b=%d p0=%d p1=%d\n",b,p0,p1);
    cairo_move_to(cr,offx+scale*px[p0],offy-scale*py[p0]);
    cairo_line_to(cr,offx+scale*px[p1],offy-scale*py[p1]);
    cairo_stroke(cr);
  }
};	// PlotGeom::plotBdry

void PlotGeom::plotChan()
{
  int nchan = geom->getNchan();
  printf("= nchan=%d\n",nchan);
  int *tchan = geom->Tchan();
  int *ispp = geom->Ispp();
  cairo_set_line_width(cr,10.0);
  for(int c=0; c < nchan; c++) {
    int t = tchan[c];
    int p0 = tri[6*t];
    int p1 = tri[6*t+1];
    int p2 = tri[6*t+2];
    double pxc = (px[p0]+px[p1]+px[p2])/3;
    double pyc = (py[p0]+py[p1]+py[p2])/3;
    printf("  c%d (%.2lf,%.2lf)\n",c,pxc,pyc);
    if(ispp[c]) cairo_set_source_rgb(cr, 1.0,0.0,0.0);
    else        cairo_set_source_rgb(cr, 0.0,0.0,1.0);
    cairo_arc(cr, offx+scale*pxc,offy-scale*pyc, 40.0, 0.0,TWOPI);
    cairo_stroke(cr);
  }
};	// PlotGeom::plotChan

void PlotGeom::plotTriNum()
{
  /*  write tri number at the center */
  cairo_set_source_rgb(cr, 0.0,0.0,0.0);
// #define SIMPLE
#ifdef SIMPLE
  cairo_select_font_face(cr,"Purisa",CAIRO_FONT_SLANT_NORMAL,
        CAIRO_FONT_WEIGHT_BOLD);
#else
  cairo_font_face_t *cff = cairo_toy_font_face_create
	("Purisa", CAIRO_FONT_SLANT_NORMAL,
        CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_face(cr,cff);
#endif

  cairo_set_font_size(cr,16);
  char tstr[10];
  for(int t=0; t < ntri; t++) {
    sprintf(tstr,"%d",t);
    int p0 = tri[6*t];
    int p1 = tri[6*t+1];
    int p2 = tri[6*t+2];
    double pxc = (px[p0]+px[p1]+px[p2])/3;
    double pyc = (py[p0]+py[p1]+py[p2])/3;

    cairo_move_to(cr,offx+scale*pxc,offy-scale*pyc);
    cairo_show_text(cr,tstr);
    cairo_stroke(cr);
  }
  printf("== plotMsh with ntri=%d\n",ntri);
#ifndef SIMPLE
  cairo_font_face_destroy(cff);
#endif
};	// PlotGeom::plotTriNum

void PlotGeom::plotCrod()
{
  /*  write Crod number  */
  cairo_font_face_t *cff = cairo_toy_font_face_create
        ("Purisa", CAIRO_FONT_SLANT_NORMAL,
        CAIRO_FONT_WEIGHT_BOLD);
  cairo_set_font_face(cr,cff);

  cairo_set_font_size(cr,16);
  cairo_set_source_rgb(cr, 0.0,0.0,0.0);
  char tstr[10];
  int *crod2d = geom->Crod2D();
  for(int t=0; t < ntri; t++) {
    int rodnum = crod2d[t];
    if(rodnum == 0) continue;
    sprintf(tstr,"Rod %d",rodnum);
    int p0 = tri[6*t];
    int p1 = tri[6*t+1];
    int p2 = tri[6*t+2];
    double pxc = (px[p0]+px[p1]+px[p2])/3;
    double pyc = (py[p0]+py[p1]+py[p2])/3;

    cairo_move_to(cr,offx+scale*pxc -20,offy-scale*pyc +20);
    cairo_show_text(cr,tstr);
    cairo_stroke(cr);
  }
  printf("== plotCrod with ntri=%d\n",ntri);
  cairo_font_face_destroy(cff);
};	// PlotGeom::plotCrod

void PlotGeom::plotSave(const char *pngfn)
{
  cairo_surface_write_to_png(surface,pngfn);
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
  printf("== plot file %s written.\n",pngfn);
  cr = NULL;
  surface= NULL;
};      // PlotGeom::plotSave
