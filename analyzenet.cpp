/************************************************************************
analyzenet - for NetFlowV1.  
Set up nodtyp, nodseg, nodnod arrays based on flowing segments if allsegs = 0,
all segments if allsegs = 1
These values are subsequently overridden by putrank
TWS 2012.
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void analyzenet()
{
	extern int nseg,nnod,nnodbc,nodsegm;
	extern int *bcnod,*bcnodname,*nodtyp,*segtyp,*nodname,*segname,*ista,*iend;
	extern int **segnodname,**nodnod,**nodseg;
	extern float *diam,*lseg,*q;
	extern float **cnode;

	int k,iseg,inod,inod1,inod2,inodbc;

	for(inod=1; inod<=nnod; inod++) nodtyp[inod] = 0;
	for(iseg=1; iseg<=nseg; iseg++)	if(segtyp[iseg] == 4 || segtyp[iseg] == 5) {
//Search for nodes corresponding to this segment
		inod = segnodname[1][iseg];	//if node names are sequential, then no search is needed
		if (inod <= nnod) {
			if (nodname[inod] == inod) {
				ista[iseg] = inod;
				goto foundit1;
			}
		}
		for(inod=1; inod<=nnod; inod++) if(nodname[inod] == segnodname[1][iseg]){
			ista[iseg] = inod;
			goto foundit1;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[1][iseg]);
	foundit1:;
		inod = segnodname[2][iseg];	//if node names are sequential, then no search is needed
		if (inod <= nnod) {
			if (nodname[inod] == inod) {
				iend[iseg] = inod;
				goto foundit2;
			}
		}
		for(inod=1; inod<=nnod; inod++) if(nodname[inod] == segnodname[2][iseg]){
			iend[iseg] = inod;
			goto foundit2;
		}
		printf("*** Error: No matching node found for nodname %i\n", segnodname[2][iseg]);
	foundit2:;
//Setup nodtyp, nodseg and nodnod
		inod1 = ista[iseg];
		inod2 = iend[iseg];
		nodtyp[inod1]++;
		nodtyp[inod2]++;
		if(nodtyp[inod1] > nodsegm)
			printf("*** Error: Too many segments connected to node %i\n", inod1);
		if(nodtyp[inod2] > nodsegm)
			printf("*** Error: Too many segments connected to node %i\n", inod2);
		nodseg[nodtyp[inod1]][inod1] = iseg;
		nodseg[nodtyp[inod2]][inod2] = iseg;
		nodnod[nodtyp[inod1]][inod1] = inod2;
		nodnod[nodtyp[inod2]][inod2] = inod1;
	}
	for (inodbc=1; inodbc<=nnodbc; inodbc++){
		bcnod[inodbc] = 0;
		for(inod=1; inod<=nnod; inod++) if(nodname[inod] == bcnodname[inodbc]){	//Search for node corresponding to this node name
			bcnod[inodbc] = inod;
			if(nodtyp[inod] != 1)
				printf("*** Error: Boundary node %i is not a 1-segment node\n", inod);
			goto foundit;
		}
		printf("*** Error: No matching node found for nodname %i, %i\n", nodname[inod],bcnodname[inodbc]);
		foundit:;
	}
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){	//calculate vessel length
		lseg[iseg] = 0.;
		for(k=1; k<=3; k++)	lseg[iseg] += SQR(cnode[k][ista[iseg]] - cnode[k][iend[iseg]]);
		lseg[iseg] = sqrt(lseg[iseg]);
	}
}