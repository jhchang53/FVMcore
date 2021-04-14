PETSC_DIR=/home/jhchang/petsc-3.14.2
PETSC_ARCH=arch-linux-c-debug

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

#  lapacke is used at MPwork to compute effective compact conductivity
LAINC = -I/usr/include/lapacke
LALIB = -llapacke -llapack -lblas

LSQL=-L/home/jhchang/sqlite -lsqlite3 -lpthread -ldl

CXXFLAGS = -g

#  Neutronics only
#  test drivers
cappsimple: cappsimple.o CoreDesc.o Geom.o GeomCAPP.o XS.o XScapp.o \
  NeutFVM.o Eigen.o EigenSIpower.o Solver.o
	-$(CLINKER) -o cappsimple -g $^ ${PETSC_KSP_LIB} $(LSQL)
cappworth: cappworth.o CoreDesc.o Geom.o GeomCAPP.o XS.o XScapp.o \
  NeutFVM.o Eigen.o EigenSIpower.o Solver.o
	-$(CLINKER) -o cappworth -g $^ ${PETSC_KSP_LIB} $(LSQL)
capptcoef: capptcoef.o CoreDesc.o Geom.o GeomCAPP.o XS.o XScapp.o \
  NeutFVM.o Eigen.o EigenSIpower.o Solver.o
	-$(CLINKER) -o capptcoef -g $^ ${PETSC_KSP_LIB} $(LSQL)
cappkin: cappkin.o CoreDesc.o Geom.o GeomCAPP.o XS.o XScapp.o \
  NeutFVM.o Eigen.o EigenSIpower.o Delay.o Solver.o
	-$(CLINKER) -o cappkin -g $^ ${PETSC_KSP_LIB} $(LSQL)
#
#  full scale test
#
cappcomp: cappcomp.o CoreDesc.o Geom.o GeomCAPP.o XS.o XScapp.o \
  NeutFVM.o Eigen.o EigenSIpower.o Delay.o Solver.o \
  Cond3D.o SegCond.o TriCond.o THpro.o MatPro.o MPwork.o MPupdate.o
	-$(CLINKER) -o cappcomp -g $^ ${PETSC_KSP_LIB} $(LSQL)
cappdyn: cappdyn.o CoreDesc.o Geom.o GeomCAPP.o XS.o XScapp.o \
  NeutFVM.o Eigen.o EigenSIpower.o Delay.o Solver.o \
  Cond3D.o SegCond.o TriCond.o  THpro.o MatPro.o MPwork.o MPupdate.o
	-$(CLINKER) -o cappdyn -g $^ ${PETSC_KSP_LIB} $(LSQL)
#  to check gap radiation effect
capprad: capprad.o CoreDesc.o Geom.o GeomCAPP.o XS.o XScapp.o \
  NeutFVM.o Eigen.o EigenSIpower.o Delay.o Solver.o \
  Cond3D.o SegCond.o TriCond.o THpro.o MatPro.o MPwork.o MPupdate.o
	-$(CLINKER) -o capprad -g $^ ${PETSC_KSP_LIB} $(LSQL)

#
#  Thermal hydraulic
#
cond3d: cond3d.o Cond3D.o Geom.o GeomCAPP.o Solver.o THpro.o MatPro.o
	-$(CLINKER) -o cond3d -g $^ ${PETSC_KSP_LIB}
segcond: segcond.o SegCond.o Solver.o
	g++ -o segcond -g $^ $(PETSC_KSP_LIB)
tricond: tricond.o TriCond.o Solver.o
	g++ -o tricond -g $^ $(PETSC_KSP_LIB)
#  multiple case version
segmult: segmult.o SegCond.o Solver.o MPwork.o MatPro.o
	-$(CLINKER) -o segmult -g $^ $(PETSC_KSP_LIB)
trimult: trimult.o TriCond.o Solver.o MPwork.o MatPro.o
	-$(CLINKER) -o trimult -g $^ $(PETSC_KSP_LIB)
# TH properties
thproflow: thproflow.o THpro.o MatPro.o
	g++ -o thproflow -g $^ -lm
thprocond: thprocond.o THpro.o MatPro.o
	g++ -o thprocond -g $^ -lm
mpwork: mpwork.o MPwork.o MatPro.o
	g++ -o mpwork -g $^ -lm $(LALIB)
# simple test
etricond: etricond.o Cond3D.o Geom.o GeomEtri.o Solver.o THpro.o MatPro.o
	-$(CLINKER) -o etricond -g $^ ${PETSC_KSP_LIB}
etriseg: etriseg.o Cond3D.o Geom.o GeomEtri.o Solver.o THpro.o MatPro.o \
 SegCond.o
	-$(CLINKER) -o etriseg -g $^ ${PETSC_KSP_LIB}
# MMR test
condall:  condall.o Cond3D.o Geom.o GeomCAPP.o Solver.o THpro.o MatPro.o \
  SegCond.o TriCond.o
	-$(CLINKER) -o condall -g $^ ${PETSC_KSP_LIB}
# check geometry
geomcapp: geomcapp.o Geom.o GeomCAPP.o PlotGeom.o
	g++ -o geomcapp -g $^ -lm -lcairo
#  core description
coredesc: coredesc.o CoreDesc.o
	g++ -o coredesc -g $^ 
# special handling to use LApack
MPwork.o: MPwork.cpp
	mpicxx -c $^ $(LAINC) $(PETSC_CC_INCLUDES)
tgz:
	tar -zcvf cond.tgz README Makefile *.h *.cpp MMR1g.db
clear:
	rm -f *.o

