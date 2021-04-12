/*
 * 	coredesc.cpp
 * 	generate paramters
 */
#include <stdio.h>
#include "CoreDesc.h"

int main()
{
  CoreDesc desc;
  double Totpow = 15.0e+6;	// (Watt)
  double Psys = 3.0e+6;		// (Pa)
  double Tin =  300.0;		// (degC)
  double Tout = 630.0;
  // double Mdot = 8.691;		// (kg/sec)
  desc.setMMR(Totpow, Psys, Tin, Tout);
  desc.print();
};
