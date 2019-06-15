/**********************************************************
picturenetwork.cpp - project network on a defined plane
Uses parameters from ContourParams.dat
Labels nodes with nodvar and segments with segvar (must be float)
Colors segments according to segvar
Generates a postscript file
Version for NetFlowV2, TWS June 2017
***********************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void picturenetwork(float *nodvar, float *segvar, const char fname[])
{
	extern int nseg,nnod;
	extern int *segtyp,*ista,*iend;
	extern float *diam, **cnode, *xsl0, *xsl1, *xsl2;

	int i, k, iseg, inod;
	float xmin,xmax,ymin,ymax,xs,ys,picfac,red,green,blue,xz,xzmin,xzmax;
	FILE *ofp;

	xmin = 0.;
	xmax = sqrt(SQR(xsl1[1] - xsl0[1]) + SQR(xsl1[2] - xsl0[2]) + SQR(xsl1[3] - xsl0[3]));
	ymin = 0.;
	ymax = sqrt(SQR(xsl2[1] - xsl0[1]) + SQR(xsl2[2] - xsl0[2]) + SQR(xsl2[3] - xsl0[3]));

	picfac = FMIN(500. / xmax, 700. / ymax);//updated April 2010
	ofp = fopen(fname, "w");
	fprintf(ofp, "%%!PS-Adobe-2.0\n");
	fprintf(ofp, "%%%%Pages: 1\n");
	fprintf(ofp, "%%%%EndComments\n");
	fprintf(ofp, "%%%%Page: 1 1\n");
	fprintf(ofp, "/mx {%g mul 50 add} def\n", picfac);
	fprintf(ofp, "/my {%g mul 50 add} def\n", picfac);
	fprintf(ofp, "/cf {closepath fill} def\n");
	fprintf(ofp, "/cs {closepath stroke} def\n");
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/n {newpath} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/sl {setlinewidth} def\n");
	fprintf(ofp, "/sc {setrgbcolor} def\n");
	fprintf(ofp, "/s {stroke} def\n");
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "4 scalefont\n");
	fprintf(ofp, "setfont\n");
	//plot vessels according to segvar
	xzmin = 1.e6;
	xzmax = -1.e6;
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		xzmin = FMIN(xzmin,segvar[iseg]);
		xzmax = FMAX(xzmax,segvar[iseg]);
	}
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		fprintf(ofp, "%g sl\n", picfac*diam[iseg]);
		//Set up colors using Matlab 'jet' scheme based on z value
				if(xzmin != xzmax) xz = (segvar[iseg] - xzmin)/(xzmax - xzmin);
				else xz = 0.75;
		blue = FMIN(FMAX(1.5-4.*fabs(xz-0.25), 0.), 1.);
				green= FMIN(FMAX(1.5-4.*fabs(xz-0.5), 0.), 1.);
				red  = FMIN(FMAX(1.5-4.*fabs(xz-0.75), 0.), 1.);
				fprintf(ofp,"%f %f %f sc\n",red,green,blue);
		xs = 0.;
		ys = 0.;
		for (i = 1; i <= 3; i++) {
			xs += (cnode[i][ista[iseg]] - xsl0[i])*(xsl1[i] - xsl0[i]) / xmax;
			ys += (cnode[i][ista[iseg]] - xsl0[i])*(xsl2[i] - xsl0[i]) / ymax;
			}
		fprintf(ofp, "%g mx %g my m ", xs, ys);
		xs = 0.;
		ys = 0.;
		for (i = 1; i <= 3; i++) {
			xs += (cnode[i][iend[iseg]] - xsl0[i])*(xsl1[i] - xsl0[i]) / xmax;
			ys += (cnode[i][iend[iseg]] - xsl0[i])*(xsl2[i] - xsl0[i]) / ymax;
		}
		fprintf(ofp, "%g mx %g my l s \n", xs, ys);
	}

//label nodes in black
	fprintf(ofp,"0 0 0 setrgbcolor\n");//black
	for(inod=1; inod<=nnod; inod++){
		xs = 0.;
		ys = 0.;
		for (i = 1; i <= 3; i++) {
			xs += (cnode[i][inod] - xsl0[i])*(xsl1[i] - xsl0[i]) / xmax;
			ys += (cnode[i][inod] - xsl0[i])*(xsl2[i] - xsl0[i]) / ymax;
		}
		//comment out next two lines to remove node numbers
		//fprintf(ofp, "%g mx %g my m ", xs + 0.5 / picfac, ys);
		//fprintf(ofp, "(%g) show\n", nodvar[inod]);
	}
//label segments in blue
	//fprintf(ofp, "0 0 1 setrgbcolor\n");//blue
	fprintf(ofp, "0 0 0 setrgbcolor\n");//black
	for(iseg=1; iseg<=nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
		xs = 0.;
		ys = 0.;
		for (i = 1; i <= 3; i++) {
			xs += ((cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2. - xsl0[i])*(xsl1[i] - xsl0[i]) / xmax;
			ys += ((cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2. - xsl0[i])*(xsl2[i] - xsl0[i]) / ymax;
		}
		//comment out next two lines to remove segment numbers
		fprintf(ofp, "%g mx %g my m ", xs + 0.5*picfac,ys);
		fprintf(ofp, "(%g) show\n",segvar[iseg]);
	}
//create a color bar
	float c;
	float cbbox = 15.; //size of boxes
	float cbx = 540; //origin of color bar
	float cby = 100;//origin of color bar
	fprintf(ofp, "0.5 setlinewidth\n");
	fprintf(ofp, "/Times-Roman findfont\n");
	fprintf(ofp, "8 scalefont\n");
	fprintf(ofp, "setfont\n");
	for(k=0; k<=10; k++){
		xz = k*0.1;
		c = xzmin + (xzmax - xzmin)*xz;
		blue = FMIN(FMAX(1.5-4.*fabs(xz-0.25), 0.), 1.);//Set up colors using Matlab 'jet' scheme
		green= FMIN(FMAX(1.5-4.*fabs(xz-0.5), 0.), 1.);
		red  = FMIN(FMAX(1.5-4.*fabs(xz-0.75), 0.), 1.);
		fprintf(ofp, "%f %f %f setrgbcolor\n",red,green,blue);
		fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
			cbx,cby+k*cbbox,cbx+cbbox,cby+k*cbbox,cbx+cbbox,cby+(k+1)*cbbox,cbx,cby+(k+1)*cbbox);
		if(k>0) fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%g) show\n",cbx+cbbox*1.1,cby+cbbox*(k-0.1),c);
	}
	fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
		cbx,cby,cbx+cbbox,cby,cbx+cbbox,cby+cbbox*11,cbx,cby+cbbox*11);
	fprintf(ofp, "showpage\n");
	fclose(ofp);
}