/************************************************************************
solve - for NetFlowV1
iterative solution for pressures in network
TWS 2012

Note: No nondimensionalization
Lengths, diameters in microns, times in s
Flows in nanoliters/min
pressures in mmHg
Viscosity in cp

Need double precision pressures for networks with many segments
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void solve(double *nodpress)
{
	extern int nseg,nnod,nnodbc,nitmax,nodsegm;
	extern int *bcnod,*bctyp,*nodtyp,*segtyp;
	extern int **nodseg,**nodnod;
	extern float tol,omega,constvisc,consthd;
	extern float *bcprfl,*cond,*hd,*nodvar,*segvar;

	int iseg,inod,inodbc,i,j,errnode,niter;
	float maxerr;
	double **wk,press1,pcondsum,condsum;
	wk = dmatrix(1,nodsegm,1,nnod);

//set up coefficients
	for(inod=1; inod<=nnod; inod++){
		if(nodtyp[inod] == 1) wk[1][inod] = 0.;
		if(nodtyp[inod] > 1){
			condsum = 0.;
			for(i=1; i<=nodtyp[inod]; i++){
				iseg = abs(nodseg[i][inod]);
				condsum += cond[iseg];
				wk[i][inod] = cond[iseg];
			}
			for(i=1; i<=nodtyp[inod]; i++) wk[i][inod] = wk[i][inod]/condsum;
		}
	}
//pressure and flow boundary nodes.  Temporarily set nodtyp of pressure nodes to -1.
	for(inodbc=1; inodbc<=nnodbc; inodbc++){
		inod = bcnod[inodbc];
		if(bctyp[inodbc] == 0){
			nodpress[inod] = bcprfl[inodbc];
			nodtyp[inod] = -1;
		}
		else wk[1][inod] = bcprfl[inodbc]/cond[nodseg[1][inod]];
	}
//iterative solution for pressures
	for(niter=1; niter<=nitmax; niter++){
		maxerr = 0.;
		for(inod=1; inod<=nnod; inod++){
			if(nodtyp[inod] == 1) press1 = omega*(nodpress[nodnod[1][inod]] + wk[1][inod] - nodpress[inod]);
			if(nodtyp[inod] >= 2){
				pcondsum = 0.;
				for(i=1; i<=nodtyp[inod]; i++) pcondsum += wk[i][inod]*nodpress[nodnod[i][inod]];
  				press1 = omega*(pcondsum - nodpress[inod]);
			}
 			if(nodtyp[inod] >= 1){
				nodpress[inod] += press1;
				if(fabs(press1) >= maxerr){
					maxerr = fabs(press1);
					errnode = inod;
				}
			}
		}
		if(maxerr < tol) goto converged;
	}
	printf("*** Warning: linear iteration not converged, maxerr = %g at node %i\n", maxerr,errnode);
	converged:;	
	for(inod=1; inod<=nnod; inod++)	if(nodtyp[inod] == -1) nodtyp[inod] = 1;//reset nodtyp of pressure nodes
	free_dmatrix(wk,1,nodsegm,1,nnod);
}
