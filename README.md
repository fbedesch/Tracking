# Tracking
Track and vertex handling code

These set of classes is used for a fast simulation of tracking resolution inside a solenoidal magnetic field. A vertex fitting class is also included.

- Geometry generation
A tracking geometry is described as a list of cylindrical layers or disks at constant z. Inert material layers are also included. The class SolGeomIDEA.cc is an example that generates the IDEA detector for FCC/CEPC. The geometry can be stored and read out in a .txt files. Examplesincluded here are:
GeoCLD.txt, GeoIDEA_BASE.txt, GeoIDEA_GT.txt

-  Track covariance matrix generation
SolTrack is a class that given a track point of origin and its correponding momentum vector calculates the helix parameters and their covariance matrix starting from the previously defined geometry class.

- Track covariance grid
The previous step of covariance calculation can be slow, so the class SolGridCov generates and reads a root file with a (pt, theta) grid of track covariances (theta is the polar angle). It can then return any covariance by doing a bilinear interpolation among the grid nodes. Examples of these root files are:
CovIDEA-BASE.root, CovCLD.root

- Smearing track parameters
ObsTrk is a class to generate resolution smeared track parameters according to the covariance matrix obtained from SolGridCov. It inherits from the class TrkUtil that contains some useful track manipulation methods.
	A test program is provided: TestObsTrk.c

- Vertex fitting 
A class VertexFit is provided to find a common vertex and its error matrix given a list of track parameters and their covariance matrices. A test program is provided: TestVertex.c

- Track utilities
A class TrkUtil is provided to help with common track related operations. It now includes the generation of the number of ionization clusters generated in a gas volume defined by two cylinders and two planes. A program to test the track length: TestTrLen.c and one to test the cluster generation: TestClCount.c are provided.  
