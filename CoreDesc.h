/*
 * 	CoreDesc.h
 * 	core description
 */
class CoreDesc {
public:
  CoreDesc();
  ~CoreDesc();
  void setPrlev(int prlev);
  void setMMR(double Totpow, double Psys, double Tin, double Tout);
  void print();
  /*  global model  */
  double Totpow,Psys,Tin,Mdot;
  double side;		// side lenght of hexagon (m)
  double fblk;  	// height of fuel zone (m_
  double fgra;  	// graphite only zone of fuel block
  double Dcool;		// coolant hole diameter (m)
  double powden;	// core power density (W/m3)
  double nchan;	// number of coolant channels in full core
  double wet_den;	// wetted perimeter density
  double ffuel;		// fraction of fuel for conductivity calc.
  double fgraf;		// fraction of graf for conductivity calc.
  /* meso model */
  double rF,rgap,rD;	// radii (m)
  double f_qcomp,qcomp;	// compact power density (W/m3)
  /*  super model  */
  double rC,rA,rB;	// radii (m)
  double ffrac;		// fraction of fuel in mixture region
  double f_qmix,qmix;	// mixture power density
  /*  triso model  */
  int ncoat;		// no. of coating layers
  double dKer;		// kernel diameter (um)
  double tcoat[10];	// coating thickness (um)
  double PF;		// packing fraction
  double f_pker,pker;	// kernel power (W)
private:
  void makeTH3D();
  void makeMeso();
  void makeSuper();
  void makeTriso();

  int prlev;
  double rt3,PI;
  int block_nfuel,block_nlayer;
  int no_pellets,no_holes;
  double block_pitch,block_height;
  double block_Dcool;
  double pellet_stack_height;
  double pellet_OD,pellet_height,block_Dfuel;
  double pellet_over;
  double cell_pitch;
  double hole_radius;
  double block_area,cell_area,compact_area;
  double triso_Dkernel;
  int triso_ncoat;
  double triso_OD,triso_tbuf, triso_tipyc, triso_tsic, triso_topyc;
  double packing_fraction;

  double triso_compact;
};
