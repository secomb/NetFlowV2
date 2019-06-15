/************************************************************************
setuparrays1 - for AngioAdapt07.  TWS January 08
Set up arrays with dimensions nnod and nseg
Revised TWS2010
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays1(int nseg, int nnod)
{
	extern int nsp,nodsegm;
	extern int *nk,*nodrank,*nodout,*nodtyp;
	extern int *ista, *iend;
	extern int **nodnod,**nodseg;
	extern float *lseg,*qold,*hdold,*cond,*qinput,*segvar,*nodvar,*segpress;
	extern float *qq, *tau, *length_weight, *histogramdisplay, *histogramweight;
	extern double *nodpress;

	ista = ivector(1,nseg); 
	iend = ivector(1,nseg); 
	nk = ivector(1,nnod); 
	nodout = ivector(1,nnod); 
	nodrank = ivector(1,nnod); 
	nodtyp = ivector(1,nnod); 

	nodnod = imatrix(1,nodsegm,1,nnod); 
	nodseg = imatrix(1,nodsegm,1,nnod); 

	cond = vector(1, nseg);
	qq = vector(1, nseg);
	tau = vector(1, nseg);
	segpress = vector(1, nseg);
	hdold = vector(1,nseg);
	lseg = vector(1,nseg); 
	qold = vector(1,nseg);
	segvar = vector(1,nseg);
	nodvar = vector(1, nnod);
	length_weight = vector(1, nnod);
	histogramdisplay = vector(1, LMAX(nseg, nnod)); 	// added for histogram display
	histogramweight = vector(1, LMAX(nseg, nnod)); 	// added for histogram display

	nodpress = dvector(1,nnod); 
}
