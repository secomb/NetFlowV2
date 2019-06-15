/************************************************************************
cmgui - for Greens.  TWS December 2011
Based on code provided by Gib Bogle
Produces greens.exelem and greens.exnode
To visualize, use CMGUI (e.g. cmgui_install.msi) at:
http://sourceforge.net/projects/cmiss/files/cmgui/cmgui-wx-2.8.0/
Start CMGUI, File > Open > Com File
Navigate to the directory contaning the data files and select greens.com.txt (as listed below).
At the bottom of the window that pops up, click the All button.  Change the view with the Layout,
Up and Front selections on the left hand side of the Graphics window.  Rotate by holding the left
button down and moving the mouse, pan with the middle button and zoom with the right button.
Use Graphics > Scene editor to change the display.  Note also the Spectrum editor.
++++++++++++++++++greens.com.txt+++++++++++++++++++++++++++++++++
# Create a material in addition to the default.
gfx cre mat gold ambient 1 0.7 0 diffuse 1 0.7 0 specular 0.5 0.5 0.5 shininess 0.8

gfx create spectrum jet
gfx modify spectrum jet clear overwrite_colour
gfx modify spectrum jet linear range 0 1 red   colour_range 0 1 ambient diffuse component 1
gfx modify spectrum jet linear range 0 1 green colour_range 0 1 ambient diffuse component 2
gfx modify spectrum jet linear range 0 1 blue  colour_range 0 1 ambient diffuse component 3

# Read in the reticular mesh (group vessels) and hide the axes.
gfx read nodes greens001.exnode
gfx read elements greens001.exelem

# The radius of the vessel is stored in component 1 of field
# 'vessel_radius', defined over the elements in the vessels group.

# Destroy the default lines.
gfx modify g_element vessels lines delete

gfx destroy node all
gfx modify g_element vessels general clear;
gfx modify g_element vessels cylinders coordinate coordinates tessellation default local circle_discretization 12 radius_scalar vessel_radius scale_factor 1 native_discretization NONE data node_colour spectrum jet
gfx modify g_element vessels node_points coordinate coordinates local glyph sphere general size "0*0*0" centre 0,0,0 font default orientation vessel_radius scale_factors "2*2*2" data node_colour spectrum jet

# Open the graphics window and turn on perspective (if desired).
gfx cre win 1
gfx mod win 1 view perspective
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <cstdio>
#include <math.h>
#include <string.h>
#include "nrutil.h"

void WriteExnodeHeader(FILE *exnode) // Write initial section of .exnode file
{
	//    fprintf(exnode, "Region: /vessels\n");
	fprintf(exnode, "Group name: vessels\n");
	fprintf(exnode, " #Fields=3\n");
	fprintf(exnode, " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n");
	fprintf(exnode, "  x.  Value index=1, #Derivatives=0, #Versions=1\n");
	fprintf(exnode, "  y.  Value index=2, #Derivatives=0, #Versions=1\n");
	fprintf(exnode, "  z.  Value index=3, #Derivatives=0, #Versions=1\n");
	fprintf(exnode, " 2) vessel_radius, coordinate, rectangular cartesian, #Components=1\n");
	fprintf(exnode, "  1.  Value index=4, #Derivatives=0, #Versions=1\n");
	fprintf(exnode, " 3) node_colour, coordinate, rectangular cartesian, #Components=3\n");
	fprintf(exnode, "  1.  Value index=5, #Derivatives=0, #Versions=1\n");
	fprintf(exnode, "  2.  Value index=6, #Derivatives=0, #Versions=1\n");
	fprintf(exnode, "  3.  Value index=7, #Derivatives=0, #Versions=1\n");
}

void WriteExelemHeader(FILE *exelem)  // Write initial section of .exelem file
{
	//   fprintf(exelem, "Region: /vessels\n");
	fprintf(exelem, "Group name: vessels\n");
	fprintf(exelem, " Shape.  Dimension=1\n");
	fprintf(exelem, " #Scale factor sets= 1\n");
	fprintf(exelem, "  l.Lagrange, #Scale factors= 2\n");
	fprintf(exelem, " #Nodes= 2\n #Fields=3\n");
	fprintf(exelem, " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n");
	fprintf(exelem, "   x.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
	fprintf(exelem, "      1.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
	fprintf(exelem, "      2.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
	fprintf(exelem, "   y.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
	fprintf(exelem, "      1.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
	fprintf(exelem, "      2.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
	fprintf(exelem, "   z.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
	fprintf(exelem, "      1.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
	fprintf(exelem, "      2.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
	fprintf(exelem, " 2) vessel_radius, coordinate, rectangular cartesian, #Components=1\n");
	fprintf(exelem, "   1.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
	fprintf(exelem, "      1.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
	fprintf(exelem, "      2.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
	fprintf(exelem, " 3) node_colour, coordinate, rectangular cartesian, #Components=3\n");
	fprintf(exelem, "   1.  l.Lagrange, no modify, standard node based.\n");
	fprintf(exelem, "     #Nodes= 2\n");
	fprintf(exelem, "      1.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n");
	fprintf(exelem, "       Scale factor indices:   1\n");
	fprintf(exelem, "      2.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n");
	fprintf(exelem, "       Scale factor indices:   2\n");
	fprintf(exelem, "   2.  l.Lagrange, no modify, standard node based.\n");
	fprintf(exelem, "     #Nodes= 2\n");
	fprintf(exelem, "      1.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n");
	fprintf(exelem, "       Scale factor indices:   1\n");
	fprintf(exelem, "      2.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n");
	fprintf(exelem, "       Scale factor indices:   2\n");
	fprintf(exelem, "   3.  l.Lagrange, no modify, standard node based.\n");
	fprintf(exelem, "     #Nodes= 2\n");
	fprintf(exelem, "      1.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n");
	fprintf(exelem, "       Scale factor indices:   1\n");
	fprintf(exelem, "      2.  #Values=1\n");
	fprintf(exelem, "       Value indices:     1\n");
	fprintf(exelem, "       Scale factor indices:   2\n");
}

void cmgui(float *segvar)
{
	extern int nseg;
	extern int *segtyp, *ista, *iend, *segname;
	extern float *diam, **cnode;

	int iseg, is, in;
	float red, green, blue, xz, xzmin, xzmax;
	char fname[80];
	FILE *exelem, *exnode;

	strcpy(fname, "Current/network.exelem");
	exelem = fopen(fname, "w");
	strcpy(fname, "Current/network.exnode");
	exnode = fopen(fname, "w");
	WriteExelemHeader(exelem);
	WriteExnodeHeader(exnode);

	xzmin = 1.e6;
	xzmax = -1.e6;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		xzmin = FMIN(xzmin, segvar[iseg]);
		xzmax = FMAX(xzmax, segvar[iseg]);
	}

	is = 0;
	in = 0;
	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		if (xzmin != xzmax) xz = (segvar[iseg] - xzmin) / (xzmax - xzmin);
		else xz = 0.75;
		blue = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.25), 0.), 1.);	//Set up colors using Matlab 'jet' scheme
		green = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.5), 0.), 1.);
		red = FMIN(FMAX(1.5 - 4.*fabs(xz - 0.75), 0.), 1.);
		is++;
		in++;
		if (fabs(segvar[iseg]) < .00000001) { red = 1; blue = 1; green = 1; }	//white
		//if(iseg == 1836 || iseg == 3431) { red = 0; blue = 0; green = 0; }	//black
		//write to elements file
		fprintf(exelem, "Element: %d 0 0\n", is);
		fprintf(exelem, "  Nodes: %d %d\n", in, in + 1);
		fprintf(exelem, "  Scale factors: 1 1\n");
		//write to nodes file
		fprintf(exnode, "Node: %d\n", in);
		fprintf(exnode, "%6.1f %6.1f %6.1f\n", cnode[1][ista[iseg]], cnode[2][ista[iseg]], cnode[3][ista[iseg]]);
		fprintf(exnode, "%6.2f\n", diam[iseg] / 2.);
		fprintf(exnode, "%6.2f %6.2f %6.2f\n", red, green, blue);
		in++;
		fprintf(exnode, "Node: %d\n", in);
		fprintf(exnode, "%6.1f %6.1f %6.1f\n", cnode[1][iend[iseg]], cnode[2][iend[iseg]], cnode[3][iend[iseg]]);
		fprintf(exnode, "%6.2f\n", diam[iseg] / 2.);
		fprintf(exnode, "%6.2f %6.2f %6.2f\n", red, green, blue);
	}
	fclose(exelem);
	fclose(exnode);
}