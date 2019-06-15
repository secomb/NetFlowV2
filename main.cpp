/************************************************************************
NetFlowV2
Program to compute flow in microvascular networks, including in-vivo
viscosity law and hematocrit partition in bifurcations.
See Pries and Secomb 2008, Handbook of Physiology, for formulas.
Reads Network.dat file, computes flows, and writes NetworkNew.dat.
Note that only segments with segtyp = 4 or 5 in network.dat are included.
This allows segments to be excluded from the network easily.
TWS, September 2012
Updated June 2019 to include generation of histograms,
and for consistency with FlowEstimateV1
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#if defined(__linux__)
	// Requires c++17 support, should be included in all current linux releases
	#include <experimental/filesystem> 
	namespace fs = std::experimental::filesystem::v1;
#elif defined(__APPLE__)
	// Requires removal of the -lstdc++fs flag from makefile
	#include <filesystem>
	namespace fs = std::filesystem;
#elif defined(_WIN32)    //Windows version
	#include <Windows.h>
#endif

void input();
void analyzenet();
void setuparrays1(int nseg, int nnod);
void flow();
void cmgui(float *segvar);
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void writeflow();
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void analyzeresults();

int max=100,nseg,nnod,nnodbc,niter,nnodfl,nsegfl,mxx,myy,mzz,nodsegm;
int nitmax1,nitmax,varyviscosity,phaseseparation;
int *segtyp,*segname,*bcnodname,*bcnod,*bctyp,*nodname,*nodtyp;
int *nodrank,*nodout,*nk,*ista,*iend;
int **nodnod,**nodseg,**segnodname;

float pi1 = atan(1.0)*4.0;
float facfp,vplas,mcvcorr,alx,aly,alz,lb,maxl;
float constvisc,consthd,tol,qtol,hdtol,omega,optw,optlam;
float *diam,*q,*qinput,*hd,*hdinput,*bcprfl,*bchd,*lseg,*qold,*hdold,*cond,*segpress;
float *bifpar,*cpar,*viscpar,*segvar,*nodvar, *xsl0, *xsl1, *xsl2;
float *qq, *tau, *length_weight, *histogramdisplay, *histogramweight;
float **cnode;
double *nodpress;

int main(int argc, char *argv[])
{
	int iseg,inod;

	//Create a Current subdirectory if needed. Copy data files to it.
#if defined(__unix__)
	if (!fs::exists("Current")) fs::create_directory("Current");
	fs::copy_file("ContourParams.dat", fs::path("Current/ContourParams.dat"), fs::copy_options::overwrite_existing);
	fs::copy_file("Network.dat", fs::path("Current/Network.dat"), fs::copy_options::overwrite_existing);
	fs::copy_file("RheolParams.dat", fs::path("Current/RheolParams.dat"), fs::copy_options::overwrite_existing);
#elif defined(_WIN32)
	BOOL NoOverwrite = FALSE;
	DWORD ftyp = GetFileAttributesA("Current\\");
	if (ftyp != FILE_ATTRIBUTE_DIRECTORY) system("mkdir Current");
	CopyFile("ContourParams.dat", "Current\\ContourParams.dat", NoOverwrite);
	CopyFile("Network.dat", "Current\\Network.dat", NoOverwrite);
	CopyFile("RheolParams.dat", "Current\\RheolParams.dat", NoOverwrite);
#endif

	input();	//read data files

	setuparrays1(nseg,nnod);

	analyzenet();

	flow();

	writeflow();	//creat new network.dat file called 'Current/NetworkNew.dat'

	analyzeresults();	//statistics and histograms of resulting flows

	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = segpress[iseg];

	cmgui(segvar);		//3D image of network

	picturenetwork(nodvar, segvar, "Current/SegPress.ps");	//2D projection of network

	return 0;
}