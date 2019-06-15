/************************************************************************
flow - for NetFlowV1
TWS 2012
Flows in nl/min, pressures in mmHg, diameters in um, viscosity in cP
For nonlinear iterations, if convergence is slow, hematocrit is increasingly underrelaxed.
This eventually forces convergence even when system is unstable due to rheological effects
(see work of R.L. Carr et al.).
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

float viscor(float d, float h);
//void analyzenet(int allsegs);
void solve(double *nodpress);
void dishem_generalized(void);
void dishem(void);
void putrank(void);

void flow()
{
	extern int nseg,nnod,varyviscosity,phaseseparation,nitmax1;
	extern int *segtyp,*ista,*iend;
	extern float constvisc,consthd,facfp,qtol,hdtol;
	extern float *diam,*q,*hd,*cond,*lseg,*qold,*hdold;
	extern double *nodpress;

	int iseg,inod,niter,errsegq,errseghd;
	float visc,maxqerr,maxhderr,qchange,hdchange,relax;

	int dishemversion = 1;	//0 for dishem, 1 for dishem_generalized

	for(iseg=1; iseg<=nseg; iseg++) q[iseg] = 0.;
	for(inod=1; inod<=nnod; inod++) nodpress[inod] = 50.;
	relax = 1.;
	for(niter=1; niter<=nitmax1; niter++){
		if(niter%5 == 0) relax = 0.8*relax;
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
			qold[iseg] = q[iseg];
			hdold[iseg] = hd[iseg];
			if(varyviscosity == 1) visc = viscor(diam[iseg],hd[iseg]);
			else visc = constvisc;
			cond[iseg] = facfp*pow(diam[iseg],4)/lseg[iseg]/visc;
		}
		solve(nodpress);
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5)
			q[iseg] = (nodpress[ista[iseg]] - nodpress[iend[iseg]])*cond[iseg];
//calculate segment hematocrits
		if(phaseseparation == 1){
			putrank();
  			if(dishemversion == 1) dishem_generalized();
			else dishem();
		}
		else for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5) hd[iseg] = consthd;
//compare hd and q with previous values
		maxqerr = 0.;
		maxhderr = 0.;
		errsegq = 0;
		errseghd = 0;
		for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
			qchange = q[iseg] - qold[iseg];
			hdchange = hd[iseg] - hdold[iseg];
			hd[iseg] = hdold[iseg] + relax*hdchange;
			if(fabs(qchange) >= maxqerr){
				maxqerr = fabs(qchange);
				errsegq = iseg;
			}
			if(fabs(hdchange) >= maxhderr){
				maxhderr = fabs(hdchange);
				errseghd = iseg;
			}
		}
		if(maxqerr < qtol && maxhderr < hdtol) goto converged;
	}
	printf("*** Warning: Nonlinear iteration not converged\n");
	printf("*** Flow error = %f at segment %i, h'crit error = %f at segment %i\n",maxqerr,errsegq,maxhderr,errseghd);
	converged:;
	printf("Flow: %i iterations\n",niter);

//recalculate segment hematocrits - otherwise errors can occur if relax < 1.  TWS  April 2011
    if(phaseseparation == 1  && relax < 0.9999){
        putrank();
  		if(dishemversion == 1) dishem_generalized();
		else dishem();
	}
}
