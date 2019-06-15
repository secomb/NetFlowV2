/************************************************************************
dishem - for AngioAdapt07
TWS October 07
Generalized version, for arbitrary node type.
For nout = 3 or more, treats as nout-1 sequential bifurcations 
JA,TWS - January 2012
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void dishem_generalized()
{
	extern int nnodbc,nseg,nnodfl,nsegfl,nodsegm;
	extern int *bcnod,*nodtyp,*nodrank,*nodout,*segtyp;
	extern int **nodseg;
	extern float *bifpar,*hd,*q,*bchd,*diam;
 
	int iseg,bnod,in,i,ibif,isegk,inod,nodt,nout,seg1,seg2;
	float hq,diaquot,hdd,x0,a,b,rbcrat,qikdash,flowsum,diammax,flow1;
	
	isegk = 0;
//boundary nodes
	for(bnod=1; bnod<=nnodbc; bnod++){
		inod = bcnod[bnod];
		if(nodout[inod] >= 1){
			iseg = nodseg[1][inod];
			hd[iseg] = bchd[bnod];
			isegk++;
		}
	}
//interior nodes
	for(in=1; in<=nnodfl; in++){
		inod = nodrank[in];
		nodt = nodtyp[inod];
		if(nodt >= 2){
			nout = nodout[inod];
			if(nout == nodt) for(i=1; i<=nodt; i++)	hd[nodseg[i][inod]] = bchd[1];//This should not normally happen
			else{
				flowsum = 0.;
				hq = 0;
				diammax = 0.;
				for(i=nout+1; i<=nodt; i++){//inflow nodes
					iseg = nodseg[i][inod];
					flowsum += fabs(q[iseg]);
					hq += hd[iseg]*fabs(q[iseg]);
					if(diam[iseg] > diammax) diammax = diam[iseg];	//maximum of the inflows
				}
				if(nout == 1) hd[nodseg[1][inod]] = hq/flowsum;	//single outflow
				else for(ibif=1; ibif<=nout-1; ibif++){			//outflow nodes
					seg1 = nodseg[ibif][inod];
					seg2 = nodseg[ibif+1][inod];
					hdd = (1. - hq/flowsum)/diammax;
					diaquot = SQR(diam[seg1]/diam[seg2]);
					a = bifpar[3]*(diaquot - 1.)/(diaquot + 1.)*hdd;	
					b = 1. + bifpar[2]*hdd;
					x0 = bifpar[1]*hdd;
					flow1 = fabs(q[seg1]);
					qikdash = (flow1/flowsum - x0)/(1. - 2.*x0);
					if(qikdash <= 0.) rbcrat = 0.;
					else if(qikdash >= 1.) rbcrat = 1.;
					else rbcrat = 1./(1. + exp(-a-b*log(qikdash/(1.-qikdash))));
					hd[seg1] = rbcrat*hq/flow1;					
					if(ibif < nout-1){
						flowsum -= flow1;		//total flow in remaining outflows
						hq = (1. - rbcrat)*hq;	//total red cell flux in remaining outflows
						diammax = diam[seg2];	//inflow diameter to next bifurcation
					}
					else hd[seg2] = (1.-rbcrat)*hq/fabs(q[seg2]);	//last remaining outflow
				}
			}
			isegk += nout;
		}
	}
	if(isegk != nsegfl)
		printf("*** Error in dishem, %i of %i segments processed\n",isegk,nsegfl);
}
