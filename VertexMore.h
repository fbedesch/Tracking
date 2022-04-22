//
#ifndef G__VERTEXMORE_H
#define G__VERTEXMORE_H
//
#include <TMath.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include "TrkUtil.h"
#include "VertexFit.h"
#include <vector>
#include <iostream>
//
// Class to include found vertices in vertex fitting

class VertexMore: public TrkUtil
{
	//
	// Vertex fitting with vertices
	// Author: F. Bedeschi, INFN-Pisa, Italy
	// April 7, 2022
	//
private:
	//
	// Inputs
	VertexFit* fV;						// Input vertex
	Bool_t fUnits;						// Default is m, can switch to mm
	Double_t fa;						// conversion constant C = a/(2Pt)
	Int_t fNtr;							// NUmber of tracks in vertex
	TVectorD fPar;						// Vertex track parameters
	TMatrixDSym fCov;					// Vertex track covariance
	std::vector<Double_t> fQ;			// Track Charges
	std::vector<TVector3*> fpi;			// Track Momenta
	std::vector<TMatrixDSym*> fCpi;		// Track Momentum errors
	TVector3 fP;						// Total momentum
	Double_t fQtot;						// Vertex charge
	TMatrixDSym fCp;					// Total momentum errors
	//	
	void CalcParCov();					// Calculate parameters and covariance
	TMatrixD dPdAlf(Int_t i);			// Derivatives of momentum wrt parameters
public:
	//
	// Constructors
	VertexMore(VertexFit *V);			// Initialize with found vertex (use meter as units)
	VertexMore(VertexFit* V, Bool_t Opt);		// Initialize with unit option TRUE = use mm, FALSE = use meters
	// Destructor
	~VertexMore();
	//
	TVectorD GetVpar();				// Get vertex track parameters
	TMatrixDSym GetVcov();				// Get vertex track covariance
	Double_t GetCharge(Int_t i) { return fQ[i]; };
	TVector3 GetMomentum(Int_t i) { return *fpi[i]; };		// Momentum of track i at vertex
	TMatrixDSym GetMomentumC(Int_t i) { return *fCpi[i]; };		// Momentum errors of track i at vertex
	TVector3 GetTotalP() { return fP; };				// Total vertex momentum
	Double_t GetTotalQ() { return fQtot; };				// Total vertex charge
	TMatrixDSym GetTotalPcov() { return fCp; };			// Total vertex momentum errors
};

#endif
