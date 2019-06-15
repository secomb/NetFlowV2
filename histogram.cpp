/*****************************************************
Evaluate histograms of solute levels.  TWS December 07.
Version 2.0, May 1, 2010.
Version 3.0, May 17,2011.
Revised Bohan Li, 2018.
Added weighting of data, June 2019.
******************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include "nrutil.h"
#include <stdio.h>
#include <string.h>

void histogram(float *var, float *weight, int n, const char filename[])
{
	int i, j, nctop, ncbot, numbins = 101, imaxstat;
	float step = 1e6, dev, xmax, ymax, scalefacx, scalefacy, xshift;
	float *stepmax, *stat, *cumu, *mstat;
	float mean = 0, min = 1e15, max = 0, maxstat = 0, weightsum = 0;
	char histogramimagename[100];
	FILE *ofp;

	strcpy(histogramimagename, filename);
	strcat(histogramimagename, ".ps");

	for (i = 1; i <= n; i++) {
		mean += var[i] * weight[i];
		weightsum += weight[i];
		if (min > var[i]) min = var[i];
		if (max < var[i]) max = var[i];
	}
	mean = mean / weightsum;

	for (i = 0; i <= 36; i++) {
		numbins = (int(floor(max / step)) - int(floor(min / step)) + 1);
		if (numbins >= 20) goto done;	//found number of bins
		if (i % 3 != 1) step = step * 0.5;
		else step = step * 0.4;
	}
	printf("Error: bin size not found for %s\n", filename);
	return;
done:;
	nctop = floor(max / step) + 1.;
	ncbot = floor(min / step);
	if (ncbot <= 2 && ncbot > 0) ncbot = 0;
	if (nctop >= -2 && nctop < 0) nctop = 0;
	numbins = nctop - ncbot;
	dev = 0.;
	stepmax = vector(1, numbins);
	stat = vector(1, numbins);
	cumu = vector(1, numbins);
	mstat = vector(1, numbins);
	for (i = 1; i <= numbins; i++) {
		stepmax[i] = step * (i + ncbot);
		mstat[i] = 0;
	}
	for (i = 1; i <= n; i++) {
		dev = dev + SQR(mean - var[i]) * weight[i];
		for (j = 1; j <= numbins; j++)	if (var[i] < stepmax[j]) {
			mstat[j] += weight[i];
			goto binned;
		}
	binned:;
	}

	dev = sqrt(dev / weightsum);
	for (i = 1; i <= numbins; i++) {
		stat[i] = mstat[i] * 100. / weightsum;
		maxstat = FMAX(maxstat, stat[i]);
	}
	imaxstat = (maxstat + 5.) / 5.;
	cumu[1] = stat[1];
	for (i = 2; i <= numbins; i++) cumu[i] = cumu[i - 1] + stat[i];
	ofp = fopen(filename, "w");
	fprintf(ofp, "Histogram data\n");
	fprintf(ofp, "value  %% cumul. %%\n");
	fprintf(ofp, "%g %7.2f %7.2f\n", step*ncbot, 0., 0.);
	for (i = 1; i <= numbins; i++) fprintf(ofp, "%g %7.2f %7.2f\n", stepmax[i], stat[i], cumu[i]);
	fprintf(ofp, "Mean = %f deviation = %f min = %f max  = %f\n", mean, dev, min, max);
	fclose(ofp);
	
	xmax = max;
	scalefacx = 500. / (numbins*step);
	ymax = imaxstat * 5.;
	scalefacy = 500. / ymax;
	xshift = (ncbot*step)*scalefacx;
	ofp = fopen(histogramimagename, "w");
	fprintf(ofp, "%%!PS\n");
	fprintf(ofp, "/mx {%g mul %f add} def\n", scalefacx, 50 - xshift);
	fprintf(ofp, "/my {%g mul %f add} def\n", scalefacy, 50.);
	fprintf(ofp, "/m {moveto} def\n");
	fprintf(ofp, "/l {lineto} def\n");
	fprintf(ofp, "/cf {1 0 360 arc closepath fill} def\n");
	fprintf(ofp, "newpath\n");
	fprintf(ofp, "%g mx %g my m\n", ncbot*step, ymax);
	fprintf(ofp, "%g mx %g my l\n", ncbot*step, 0.);
	fprintf(ofp, "%g mx %g my l\n", nctop*step, 0.);
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "/Times-Roman findfont 10 scalefont setfont\n");
	//fprintf(ofp, "%d %d (%s) show\n", 600, 250, title);
	for (i = 0; i <= imaxstat; i++) fprintf(ofp, "30 %f my m (%4.1f) show\n", i * 5., i * 5.);	//y-axis labels
	for (i = 0; i <= numbins; i++)	if (i % 2 == 0) fprintf(ofp, "%f mx -9 add 40 m (%3.1f) show\n", step*(ncbot + i), step*(ncbot + i));	//x-axis labels (used -9 add to center labels across bars)
	fprintf(ofp, "stroke\n");	//needed to avoid artifact in plot
	fprintf(ofp, "0 1 0 setrgbcolor\n");
	fprintf(ofp, "newpath\n");
	for (i = 1; i <= numbins; i++) {
		fprintf(ofp, "%g mx %g my m\n", (i - 1 + ncbot)*step, 0.);
		fprintf(ofp, "%g mx %g my l\n", (i - 1 + ncbot)*step, stat[i]);
		fprintf(ofp, "%g mx %g my l\n", (i + ncbot)*step, stat[i]);
		fprintf(ofp, "%g mx %g my l\n", (i + ncbot)*step, 0.);
	}
	fprintf(ofp, "stroke\n");
	fprintf(ofp, "showpage\n");
	fclose(ofp);

	free_vector(stepmax, 1, numbins);
	free_vector(stat, 1, numbins);
	free_vector(cumu, 1, numbins);
	free_vector(mstat, 1, numbins);
}