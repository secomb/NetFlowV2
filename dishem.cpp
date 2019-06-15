/************************************************************************
dishem - for AngioAdapt07
TWS October 07
Hematocrit partition according to Pries and Secomb, Handbook of Physiology 2008
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void dishem()
{
	extern int nnodbc,nseg,nnodfl,nodsegm;
	extern int *bcnod,*nodtyp,*nodrank,*nodout,*segtyp;
	extern int **nodseg;
	extern float *bifpar,*hd,*q,*bchd,*diam;
 
	int iseg,bnod,in,i,isegk,inod,nodt,nin,nout;
	int *segs;
	float hq,diaquot,hdd,x0,a,b,rbcrat,qikdash,flowsum;
	float *flow;

	segs = ivector(1,nodsegm);
	flow = vector(1,nodsegm);

	isegk = 0;
//boundary nodes
	for(bnod=1; bnod<=nnodbc; bnod++){
		inod = bcnod[bnod];
		if(nodout[inod] == 1){
			iseg = nodseg[1][inod];
			hd[iseg] = bchd[bnod];
			isegk++;
		}
	}
//interior nodes
	for(in=1; in<=nnodfl; in++){
		inod = nodrank[in];
		nodt = nodtyp[inod];
		nout = nodout[inod];
		nin = nodt - nout;
		if(nodt >= 2) for(i=1; i<=nodt; i++){
			segs[i] = nodseg[i][inod];
			flow[i] = fabs(q[segs[i]]);
		}
//2-segment nodes
		if(nodt == 2){
			if(nout == 2) hd[segs[2]] = bchd[1];	//nout = 2 could happen during iteration
			hd[segs[1]] = hd[segs[2]];
		}
//3-segment nodes
		if(nodt == 3){
			if(nout == 1) hd[segs[1]] = (flow[2]*hd[segs[2]] + flow[3]*hd[segs[3]])/flow[1];	//convergent
			if(nout == 2){	//divergent
				hdd = (1. - hd[segs[3]])/diam[segs[3]];
				diaquot = SQR(diam[segs[1]]/diam[segs[2]]);
				a = bifpar[3]*(diaquot - 1.)/(diaquot + 1.)*hdd;	
				b = 1. + bifpar[2]*hdd;
				x0 =  bifpar[1]*hdd;
				qikdash = (flow[1]/flow[3] - x0)/(1. - 2.*x0);
				if(qikdash <= 0.){
					hd[segs[1]] = 0;
					hd[segs[2]] = hd[segs[3]]*flow[3]/flow[2];
				}
				else if(qikdash >= 1.){
					hd[segs[2]] = 0;
					hd[segs[1]] = hd[segs[3]]*flow[3]/flow[1];
				}
				else{
					rbcrat = 1./(1. + exp(-a-b*log(qikdash/(1.-qikdash))));
					hd[segs[1]] = rbcrat*hd[segs[3]]*flow[3]/flow[1];
					hd[segs[2]] = (1.-rbcrat)*hd[segs[3]]*flow[3]/flow[2];
				}
			}
			if(nout == 3) for(i=1; i<=3; i++) hd[segs[i]] = bchd[1];
		}
//For nodes with more than three segments, no phase separation
		if(nodt > 3){
			if(nout == nodt) for(i=1; i<=nodt; i++) hd[segs[i]] = bchd[1];
			else{
				flowsum = 0.;
				hq = 0;
				for(i=nout+1; i<=nodt; i++){
					flowsum = flowsum + fabs(q[segs[i]]);
					hq = hq + hd[segs[i]]*fabs(q[segs[i]]);
				}
				if(nout >= 1) for(i=1; i<=nout; i++) hd[segs[i]] = hq/flowsum;
			}
		}
		if(nodt >= 2) isegk += nout;
	}
	if(isegk != nseg)
		printf("*** Error in dishem, %i of %i segments processed\n",isegk,nseg);
	free_ivector(segs,1,nodsegm);
	free_vector(flow,1,nodsegm);
}