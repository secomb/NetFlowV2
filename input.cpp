/************************************************************************
input - for NetFlowV1
TWS 2012
Pressures in mmHg, flows in um^3/s, viscosity in cP
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void input()
{
	extern int mxx,myy,mzz,nseg,nnod,nnodbc,nodsegm;
	extern int nitmax1,nitmax,varyviscosity,phaseseparation;
	extern int *segname,*segtyp,*bcnodname,*bcnod,*bctyp,*nodname;
	extern int **segnodname;
 
	extern float pi1,facfp,vplas,mcvcorr;
	extern float alx,aly,alz,lb,maxl;
	extern float tol,qtol,hdtol,omega,optw,optlam,constvisc,consthd;
	extern float *diam,*lseg,*q,*qinput,*hd,*hdinput,*bcprfl,*bchd;
	extern float *bifpar,*cpar,*viscpar,*xsl0, *xsl1, *xsl2;
	extern float **cnode;

	float mcv;
	int i,inodbc,iseg,max=200;
	FILE *ifp;
	char bb[200];

	bifpar = vector(1,3);
	cpar = vector(1,4);
	viscpar = vector(1,6);

	ifp = fopen("RheolParams.dat", "r");
	fscanf(ifp,"%f %f %f%*[^\n]", &bifpar[1],&bifpar[2],&bifpar[3]);
	fscanf(ifp,"%f %f %f %f%*[^\n]", &cpar[1],&cpar[2],&cpar[3],&cpar[4]);
	fscanf(ifp,"%f %f %f %f %f %f%*[^\n]", &viscpar[1],&viscpar[2],&viscpar[3],&viscpar[4],&viscpar[5],&viscpar[6]);
	fscanf(ifp,"%i %f %f%*[^\n]", &nitmax,&tol,&omega);
	fscanf(ifp,"%i %f %f%*[^\n]", &nitmax1,&qtol,&hdtol);
	fscanf(ifp,"%f %f%*[^\n]", &optw,&optlam);
	fscanf(ifp,"%f %f %f%*[^\n]", &constvisc,&vplas,&mcv);
	fscanf(ifp,"%f%*[^\n]", &consthd);
	fscanf(ifp,"%i%*[^\n]", &varyviscosity);
	fscanf(ifp,"%i", &phaseseparation);
	fclose(ifp);

//conversion factor for flow from pressure based on Poiseuille flows measured in nl/min
	facfp = pi1*1333./128./0.01*60./1.e6;
	mcvcorr = pow(92./mcv,0.33333);

//network data file
	ifp = fopen("Network.dat", "r");
	fgets(bb,max,ifp);
	printf("%s\n",bb);
	fscanf(ifp, "%f %f %f%*[^\n]", &alx,&aly,&alz);
	fscanf(ifp, "%i %i %i%*[^\n]", &mxx,&myy,&mzz);
	fscanf(ifp, "%f%*[^\n]", &lb);
	fscanf(ifp, "%f%*[^\n]", &maxl);
	fscanf(ifp, "%i%*[^\n]", &nodsegm);
	fscanf(ifp, "%i%*[^\n]", &nseg);
	fgets(bb,max,ifp);
	fgets(bb,max,ifp);
	segname = ivector(1,nseg);
	segtyp = ivector(1,nseg);
	segnodname = imatrix(1,2,1,nseg);
	diam = vector(1,nseg);
	qinput = vector(1,nseg);
	q = vector(1,nseg);
	hd = vector(1,nseg);
	hdinput = vector(1,nseg);
	for(iseg=1; iseg<=nseg; iseg++)	fscanf(ifp, "%i %i %i %i %f %f %f%*[^\n]", 
		&segname[iseg],&segtyp[iseg],&segnodname[1][iseg],&segnodname[2][iseg],&diam[iseg],&qinput[iseg],&hdinput[iseg]);
//number of nodes in vessel network
	fgets(bb,max,ifp);
	fscanf(ifp,"%i%*[^\n]", &nnod);
	fgets(bb,max,ifp);
	fgets(bb,max,ifp);
//coordinates of nodes
	nodname = ivector(1,nnod);
	cnode = matrix(1,3,1,nnod);
	for(i=1; i<=nnod; i++) fscanf(ifp, "%i %f %f %f%*[^\n]", &nodname[i],&cnode[1][i],&cnode[2][i],&cnode[3][i]);
//boundary nodes
	fscanf(ifp,"%i%*[^\n]", &nnodbc);
	fgets(bb,max,ifp);
	fgets(bb,max,ifp);
	bcnodname = ivector(1,nnodbc);
	bcnod = ivector(1,nnodbc);
	bctyp = ivector(1,nnodbc);
	bcprfl = vector(1,nnodbc);
	bchd = vector(1,nnodbc);
	for(inodbc=1; inodbc<=nnodbc; inodbc++)	fscanf(ifp,"%i %i %f %f%*[^\n]", &bcnodname[inodbc],&bctyp[inodbc],&bcprfl[inodbc],&bchd[inodbc]);
	fclose(ifp);

	//Read parameters for drawing segments and nodes
	xsl0 = vector(1, 3);
	xsl1 = vector(1, 3);
	xsl2 = vector(1, 3);
	ifp = fopen("ContourParams.dat", "r");
	fscanf(ifp, "%f %f %f%*[^\n]", &xsl0[1], &xsl0[2], &xsl0[3]);
	fscanf(ifp, "%f %f %f%*[^\n]", &xsl1[1], &xsl1[2], &xsl1[3]);
	fscanf(ifp, "%f %f %f%*[^\n]", &xsl2[1], &xsl2[2], &xsl2[3]);
	fclose(ifp);
}