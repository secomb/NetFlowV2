/************************************************************************
putrank - for NetFlowV1.  TWS 2012
Generate list of nodes in order of flow direction
Considers only type 4 and 5 segments
If nodrank[i] < nodrank[j], node j is not upstream of node i
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void putrank(void)
{
	extern int nseg,nnod,nnodfl,nsegfl;
	extern int *nodrank,*nodtyp,*nodout,*nodname,*segtyp,*nk,*ista,*iend;
	extern int **nodseg,**nodnod;
	extern float *q;

	int inod,j,iseg,nod1,nod2,flag;

	for(inod=1; inod<=nnod; inod++){
		nodtyp[inod] = 0;
		nodout[inod] = 0;
	}
	nsegfl = 0;	//added TWS 2010
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		if(q[iseg] >= 0){
			nod1 = ista[iseg];
			nod2 = iend[iseg];
		}
		else{
			nod1 = iend[iseg];
			nod2 = ista[iseg];
		}
		nodtyp[nod1]++;
		nodseg[nodtyp[nod1]][nod1] = iseg;
		nodnod[nodtyp[nod1]][nod1] = nod2;
		nodout[nod1]++;
		nsegfl++;
	}
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		if(q[iseg] >= 0){
			nod1 = ista[iseg];
			nod2 = iend[iseg];
		}
		else{
			nod1 = iend[iseg];
			nod2 = ista[iseg];
		}
		nodtyp[nod2]++;
		nodseg[nodtyp[nod2]][nod2] = iseg;
		nodnod[nodtyp[nod2]][nod2] = nod1;
	}	
//assign low ranks to inflow nodes
	nnodfl = 0;
	for(inod=1; inod<=nnod; inod++){
		nk[inod] = 0;
 		if(nodtyp[inod] == 1 && nodout[inod] == 1){
			nnodfl++;
			nk[inod] = 1;
			nodrank[nnodfl] = inod;
		}
	}
//assign increasing ranks to downstream connected nodes
	flag = 1;
	while(flag == 1){
		flag = 0;
		for(inod=1; inod<=nnod; inod++)	if(nk[inod] == 0 && nodtyp[inod] > 0){
			for(j=nodout[inod]+1; j<=nodtyp[inod]; j++){
				iseg = nodseg[j][inod]; 
				if(inod == iend[iseg] && (nk[ista[iseg]] == 0 || q[iseg] <= 0.)) goto skipnode;
				if(inod == ista[iseg] && (nk[iend[iseg]] == 0 || q[iseg] >= 0.)) goto skipnode;
			}
			nnodfl++;
			nk[inod] = 1;
			nodrank[nnodfl] = inod;
			flag = 1;
			skipnode:;
		}
	}
//check for unprocessed nodes--should be none.
	for(inod=1; inod<=nnod; inod++)	if(nodtyp[inod] != 0 && nk[inod] == 0)
		printf("*** Error: unprocessed node %i in putrank\n",inod);
}
