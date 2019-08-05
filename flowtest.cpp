/************************************************************************
flowtest - for FlowEstimateV1
Use linear solver to identify segments without flow and set segtyp to 1
Set conductances to constant values for this calculation
Since not all boundary conditions are known, random pressures are assigned to all boundary nodes
TWS June 2019
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"

void flowtest()
{
	extern int nseg, nnod;
	extern int *segtyp, *nodtyp, *ista, *iend, *nodname, *segname;
	extern int **nodnod, **segnodname;
	extern float *q, *diam, *hd, **cnode;
	extern double *nodpress;

	int iseg, iseg1, inod, inod1, inod2, i, errnode, niter, nnflow, nitmax = 10000, nnod0;
	float maxerr, tol = 1.e-12, omega = 1.95;
	double press1, pcondsum;

	//basic version of analyzenet
	for (iseg = 1; iseg <= nseg; iseg++) {		//Search for nodes corresponding to this segment
		for (inod = 1; inod <= nnod; inod++) if (nodname[inod] == segnodname[1][iseg]) {
			ista[iseg] = inod;
			goto foundit1;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[1][iseg]);
	foundit1:;
		for (inod = 1; inod <= nnod; inod++) if (nodname[inod] == segnodname[2][iseg]) {
			iend[iseg] = inod;
			goto foundit2;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[2][iseg]);
	foundit2:;
	}
	for (inod = 1; inod <= nnod; inod++) nodtyp[inod] = 0;
	for (iseg = 1; iseg <= nseg; iseg++) {
		inod1 = ista[iseg];
		inod2 = iend[iseg];
		nodtyp[inod1]++;
		nodtyp[inod2]++;
		nodnod[nodtyp[inod1]][inod1] = inod2;
		nodnod[nodtyp[inod2]][inod2] = inod1;
	}
	//basic version of solve
	for (inod = 1; inod <= nnod; inod++) {
		if (nodtyp[inod] == 1) nodpress[inod] = rand()*100. / RAND_MAX;	//put random pressures on all type 1 nodes
		else nodpress[inod] = 50.;
	}
	for (niter = 1; niter <= nitmax; niter++) {					//iterative solution for pressures
		maxerr = 0.;
		for (inod = 1; inod <= nnod; inod++) if (nodtyp[inod] > 1) {//don't process type 1 nodes
			pcondsum = 0.;
			for (i = 1; i <= nodtyp[inod]; i++) pcondsum += nodpress[nodnod[i][inod]] / (float)nodtyp[inod];
			press1 = omega * (pcondsum - nodpress[inod]);
			nodpress[inod] += press1;
			if (fabs(press1) >= maxerr) {
				maxerr = fabs(press1);
				errnode = inod;
			}
		}
		if (maxerr < tol) goto converged;
	}
	printf("*** Warning: Flowtest - linear iteration not converged, maxerr = %g at node %i\n", maxerr, errnode);
converged:;
	printf("Flowtest: %i iterations, maxerr = %e\n", niter, maxerr);

	nnflow = 0;
	for (iseg = 1; iseg <= nseg; iseg++) {
		if (fabs(nodpress[ista[iseg]] - nodpress[iend[iseg]]) < 1.e-6) {	//remove segment
			nnflow++;
			if (nnflow == 1) printf("Zero-flow segments set to type 0\n");
			printf("%i ", iseg);
			if (nnflow % 20 == 0) printf("\n");
			segtyp[iseg] = 0;
		}
	}
	if (nnflow > 0) printf("\n%i segments set to type 0\n", nnflow);
}