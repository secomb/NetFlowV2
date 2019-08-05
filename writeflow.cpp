/************************************************************************
writeflow - for NetFlowV1
rewrite network.dat with updated flows
TWS 2012
Pressures in mmHg, flows in um^3/s, viscosity in cP
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void writeflow()
{
	extern int nseg,*segname,*segtyp,*nodname;
	extern int **segnodname;
	extern float *diam,*q,*qinput,*hd,*hdinput;

	int i,iseg,max=200,length,idum;
	float fdum;
	FILE *ifp,*ofp;
	char bb[200];

//network data file
	ifp = fopen("Network.dat", "r");
	ofp = fopen("Current/NetworkNew.dat", "w");
	for(i=1; i<=8; i++){ 
		fgets(bb,max,ifp);
		fprintf(ofp,"%s",bb);
	}
	for(iseg=1; iseg<=nseg; iseg++){
		fscanf(ifp, "%i %i %i %i %f %f %f", 
			&idum,&idum,&idum,&idum,&fdum,&fdum,&fdum);
		fgets(bb,max,ifp);
		fprintf(ofp, "%i %i %i %i %f %f %f", 
			segname[iseg],segtyp[iseg],segnodname[1][iseg],segnodname[2][iseg],diam[iseg],q[iseg],hd[iseg]);
		fprintf(ofp,"%s",bb);
	}
//copy rest of file
	do{ 
		if(fgets(bb,max,ifp) != NULL){
			length = strlen(bb);
			fprintf(ofp,"%s",bb);			
		}
		else length = 1;
	}
	while(length > 1);
	fclose(ifp);
	fclose(ofp);
}