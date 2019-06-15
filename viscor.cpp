/************************************************************************
viscor - for NetFlowV1
Computation of segment viscosity based on Pries et al. Circ Res. 75,
904-915, 1994. Diameter corrected for differences in MCV
TWS 2012
************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float viscor(float d, float h)
{
	extern float mcvcorr,optw,vplas;
	extern float *cpar,*viscpar;

	float dcorr,c,eta45,etarel,hdfac,hdref=0.45;
	dcorr = d*mcvcorr;
	c = (cpar[1] + exp(cpar[2]*dcorr))*
		(-1. + 1./(1. + pow(10.F,cpar[3])*pow(dcorr,cpar[4])))
		+ 1./(1. + pow(10.F,cpar[3])*pow(dcorr,cpar[4]));
	eta45 = viscpar[1]*exp(viscpar[2]*dcorr)
		+ viscpar[3] + viscpar[4]*exp(viscpar[5]*pow(dcorr,viscpar[6]));
	hdfac = (pow(1.F - h,c) - 1.)/(pow(1.F - hdref,c) - 1.);
	etarel = (1. + (eta45 - 1.)*hdfac*SQR(dcorr/(dcorr-optw)))*SQR(dcorr/(dcorr-optw));
	return etarel*vplas;
}
