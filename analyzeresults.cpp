/************************************************************************
analyzeresults for NetFlowV2
TWS, August 2018, May 2019.
*************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nrutil.h"

void histogram(float *var, float *weight, int n, const char filename[]);
float viscor(float d, float h);

void analyzeresults()
{
	extern int nnod, nseg, varyviscosity, *ista, *iend;
	extern int *nodname, *nodtyp, *segname, *segtyp, **nodseg;
	extern float constvisc, pi1;
	extern float *histogramdisplay, *histogramweight, *q, *qq, *diam, *tau, *segpress, *hd;
	extern float *length_weight, *lseg, *cond;
	extern double *nodpress;

	int i, iseg, inod, numberzeroflows = 0;
	int minpressnod, maxpressnod, maxflowseg, minflowseg, maxshearseg, minshearseg;
	float meanflow, meannodpress, meansegpress, meanshear, nodpressdeviation, segpressdeviation, sheardeviation;
	float maxflow, minflow, maxpress, minpress, maxshear, minshear, totallength, totallength1;
	FILE *ofp;

	totallength = 0.;
	for (iseg = 1; iseg <= nseg; iseg++) {
		totallength += lseg[iseg];
		segpress[iseg] = (nodpress[ista[iseg]] + nodpress[iend[iseg]]) / 2;
		qq[iseg] = fabs(q[iseg]);
		tau[iseg] = (nodpress[ista[iseg]] - nodpress[iend[iseg]])*1333.*diam[iseg] / lseg[iseg] / 4.;
		if (qq[iseg] < 1.e-4) {
			segtyp[iseg] = 0;		//label zero-flow segments as type 0
			numberzeroflows++;
		}
	}
	printf("Zero flow segments (final): %4i\n", numberzeroflows);

	//node weighting factors: 1/2 the sum of the lengths of the segments connected to the node
	for (inod = 1; inod <= nnod; inod++) {
		length_weight[inod] = 0;
		for (iseg = 1; iseg <= nodtyp[inod]; iseg++) length_weight[inod] += 0.5*lseg[nodseg[iseg][inod]];
	}

	for (inod = 1; inod <= nnod; inod++) {
		histogramdisplay[inod] = nodpress[inod];
		histogramweight[inod] = length_weight[inod];
	}
	histogram(histogramdisplay, histogramweight, nnod, "Current/histogram-pressures.out");

	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) if(qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(qq[iseg]);
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "Current/histogram-logflows.out");

	i = 0;	//velocities in mm/s
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) if (qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(4000. / 60. * qq[iseg] / pi1 / SQR(diam[iseg]));
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "Current/histogram-logvelocities.out");

	i = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) if (qq[iseg] > 0.) {
		i++;
		histogramdisplay[i] = log10(fabs(tau[iseg]));
		histogramweight[i] = lseg[iseg];
	}
	histogram(histogramdisplay, histogramweight, i, "Current/histogram-logstress.out");

	meannodpress = 0.;
	minpress = 1.e6;
	maxpress = -1.e6;
	for (inod = 1; inod <= nnod; inod++) {		// calculate average and max/min pressures
		meannodpress += nodpress[inod] * length_weight[inod];
		if (nodpress[inod] < minpress) {
			minpress = nodpress[inod];
			minpressnod = inod;
		}
		if (nodpress[inod] > maxpress) {
			maxpress = nodpress[inod];
			maxpressnod = inod;
		}
	}
	meannodpress = meannodpress / totallength;
	nodpressdeviation = 0.;
	for (inod = 1; inod <= nnod; inod++) nodpressdeviation += length_weight[inod] * SQR(nodpress[inod] - meannodpress);
	nodpressdeviation = sqrt(nodpressdeviation / totallength);

	totallength1 = 0.;
	meanflow = 0.;
	meansegpress = 0.;
	meanshear = 0.;
	maxflow = 0.;
	minflow = 1.e6;
	maxshear = 0.;
	minshear = 1.e6;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {		//mean calculations
		totallength1 += lseg[iseg];
		meanflow += qq[iseg] * lseg[iseg];
		meansegpress += segpress[iseg] * lseg[iseg];
		meanshear += fabs(tau[iseg]) * lseg[iseg];
		if (fabs(maxflow) < qq[iseg]) {
			maxflow = qq[iseg];
			maxflowseg = iseg;
		}
		if (fabs(minflow) > qq[iseg]) {
			minflow = qq[iseg];
			minflowseg = iseg;
		}
		if (fabs(maxshear) < fabs(tau[iseg])) {
			maxshear = fabs(tau[iseg]);
			maxshearseg = iseg;
		}
		if (fabs(minshear) > fabs(tau[iseg])) {
			minshear = fabs(tau[iseg]);
			minshearseg = iseg;
		}
	}
	meansegpress = meansegpress / totallength1;
	meanflow = meanflow / totallength1;
	meanshear = meanshear / totallength1;
	segpressdeviation = 0.;
	sheardeviation = 0.;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {					// standard deviation calculations
		segpressdeviation += lseg[iseg] * SQR(segpress[iseg] - meansegpress);
		sheardeviation += lseg[iseg] * SQR(fabs(tau[iseg]) - meanshear);
	}
	segpressdeviation = sqrt(segpressdeviation / totallength1);
	sheardeviation = sqrt(sheardeviation / totallength1);

	ofp = fopen("Current/Results.out", "w");
	fprintf(ofp, "All quantities weighted by segment length\n");
	fprintf(ofp, "Segment flow mean = %6f\n", meanflow);
	fprintf(ofp, "Nodal pressure mean +- s.d.: %6f +- %6f\n", meannodpress, nodpressdeviation);
	fprintf(ofp, "Segment pressure mean +- s.d.: %6f +- %6f\n", meansegpress, segpressdeviation);
	fprintf(ofp, "Shear stress mean +- s.d.: %6f +- %6f\n", meanshear, sheardeviation);
	fprintf(ofp, "Maximum pressure: %6f at node %4i\n", maxpress, nodname[maxpressnod]);
	fprintf(ofp, "Minimum pressure: %6f at node %4i\n", minpress, nodname[minpressnod]);
	fprintf(ofp, "Maximum flow: %6f at segment %4i\n", maxflow, segname[maxflowseg]);
	fprintf(ofp, "Minimum flow: %6f at segment %4i\n", minflow, segname[minflowseg]);
	fprintf(ofp, "Maximum shear stress: %6f at segment %4i\n", maxshear, segname[maxshearseg]);
	fprintf(ofp, "Minimum shear stress: %6f at segment %4i\n", minshear, segname[minshearseg]);
	fprintf(ofp, "segment name flow    length    conductance pressure shear stress dishem visc\n");
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		fprintf(ofp, "%4i %4i %6f %6f %10.2e %6f %6f %6f ", iseg, segname[iseg], q[iseg], lseg[iseg], cond[iseg], segpress[iseg], tau[iseg], hd[iseg]);
		if (varyviscosity == 1) fprintf(ofp, "%6f\n", viscor(diam[iseg], hd[iseg]));
		else fprintf(ofp, "%6f\n", constvisc);
	}
	fclose(ofp);
}