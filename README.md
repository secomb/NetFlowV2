# NetFlowV2
Simulation of microvascular network hemodynamics

This program is designed to simulate the dynamics of blood flow in networks of microvessels. Such simulations are a necessary component in theoretical models for several aspects of the microcirculation, including transport of oxygen and other solutes [1], structural adaptation of microvessels [2–4], and regulation of blood flow [5].

The simulation is based on a series of studies, in which we developed empirical equations to describe the dynamics of blood flow in microvascular networks. The key components are:
- equations to describe the apparent viscosity of blood in microvessels in vivo as a function of vessel diameter and discharge hematocrit [6,7];
- equations to describe the partition of hematocrit in diverging bifurcations, as a function of the flow rates in the branches, the diameters of the branches, and the hematocrit arriving at the parent vessel [7,8].
  
The parameters describing these relationships are given in RheolParams.dat. This gives parameters needed for estimating the apparent viscosity in each segment. Set varyviscosity = 1 to get diameter-dependent viscosity. Set phaseseparation = 1 to compute phase separation in diverging bifurcations. This can lead to nonconvergence of the method in some cases. Also, the number of iterations for linear and nonlinear loops are specified. These may need to be increased for very large networks.

The equations used in the simulations are the versions given in our article in the Handbook of Physiology [9]. These equations allow the prediction of flow rates, pressures and hematocrits throughout the network, provided that the following information is available:
- the network topology;
- the lengths and diameters of each segment;
- the boundary conditions at every boundary segment (i.e. segment flowing in or out of the network), in the form of specified pressure or flow, including at least one pressure boundary condition;
- specified hematocrits at each inflowing boundary segment.

In the implementation provided here, it is assumed that the spatial coordinates of each node are known. The length of each segment is computed as a straight-line distance. Curved segments can be approximated by combining several straight segments.

The above information is conveniently represented in a file typically named network.dat (see example). On output, the file is renamed NetworkNew.dat, and contains the computed flow and hematocrit values for each segment. This file consists of four parts:
- header and global parameters. Except for nseg (number of segments) these parameters are ignored by NetFlow; they are used for simulation of solute trabnsport (see GreensV4).
- list of segments with the following details for each segment:
	- name - integer, not necessarily sequential
	- type - integer, if type is not 4 or 5, the segment is ignored
	- from - integer, 'from' node
	- to - integer, 'to' node
	- diam. - real, diameter in microns
	- flow - flow in nl/min.  On input, the value is ignored. On output, the computed value is given
	- hem. - hematocrit from 0 to 1.  On input, the value is ignored. On output, the computed value is given
- list of nodes with the following details for each node:
	- name - integer, not necessarily sequential
	- x,y,z - Cartesian coordinates of node in micron
- list of boundary nodes with the following details for each node:
	- name - integer, as in node list
	- bctyp - integer, 0 for defined pressure, 2 for defined inflow rate
	- press/flow - boundary pressure in mmHg or flow in nl/min
	- HD - inflow discharge hematocrit. Ignored for outflow nodes
	- PO2, solute 3, etc. - solute data, ignored by this program

The program generates several postscript files:
- histograms of node pressure (mmHg), segment flow velocities (mm/s), segment flow rates (nl/min) and shear stress (dyn/cm2)
- network geometry visualized as projected on a plane. The position and orientation of the plane are specified by ContourParams.dat, giving the coordinates of three corners. Each vessel segment is color-coded according to midpoint pressure.

The program generates greens.exelem and greens.exenode, allowing 3D visualization of the network using cmgui. When cmgui is started, use File - Open - com file, select greens.com.txt and hit “All” for the visualization. For details of cmgui, see: http://sourceforge.net/projects/cmiss/files/cmgui/cmgui-wx-2.8.0/

This code is free to use at your own risk. Feedback and/or acknowledgement are appreciated. Email secomb@u.arizona.edu.

Note: This code makes use of nrutil.h and nrutil.cpp as placed in the public domain by Numerical Recipes at http://apps.nrbook.com/c/index.html.

1. Secomb TW, Hsu R, Park EY, Dewhirst MW (2004) Green's function methods for analysis of oxygen delivery to tissue by microvascular networks. Ann Biomed Eng 32: 1519-1529.

2. Pries AR, Secomb TW, Gaehtgens P (1998) Structural adaptation and stability of microvascular networks: theory and simulations. Am J Physiol 275: H349-H360.

3. Pries AR, Reglin B, Secomb TW (2001) Structural adaptation of microvascular networks: functional roles of adaptive responses. Am J Physiol Heart Circ Physiol 281: H1015-H1025.

4. Pries AR, Reglin B, Secomb TW (2005) Remodeling of blood vessels: responses of diameter and wall thickness to hemodynamic and metabolic stimuli. Hypertension 46: 725-731.

5. Roy TK, Pries AR, Secomb TW (2012) Theoretical comparison of wall-derived and erythrocyte-derived mechanisms for metabolic flow regulation in heterogeneous microvascular networks. Am J Physiol Heart Circ Physiol 302: H1945-H1952.

6. Pries AR, Secomb TW, Gessner T, Sperandio MB, Gross JF, Gaehtgens P (1994) Resistance to blood flow in microvessels in vivo. Circ Res 75: 904-915.

7. Pries AR, Secomb TW (2005) Microvascular blood viscosity in vivo and the endothelial surface layer. Am J Physiol Heart Circ Physiol 289: H2657-H2664.

8. Pries AR, Ley K, Claassen M, Gaehtgens P (1989) Red cell distribution at microvascular bifurcations. Microvasc Res 38: 81-101.

9. Pries AR, Secomb TW (2008) Blood Flow in Microvascular Networks. In: Tuma RF, Duran WN, Ley K, editors. Handbook of Physiology: Microcirculation. Second Edition. San Diego: Academic Press. pp. 3-36.