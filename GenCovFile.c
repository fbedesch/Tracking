#include <TMath.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>
#include <iostream>
#include "SolGeom.h"
#include "SolGridCov.h"

//
void GenCovFile(char * GeoFile, TString CovFile)
{
	//
	//	Init geometry
	//
	SolGeom *G = new SolGeom(GeoFile);
	//
	// Write grid
	SolGridCov* GC = new SolGridCov();
	GC->Write(CovFile, G);
}