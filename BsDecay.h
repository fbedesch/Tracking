//
#ifndef G__BSDECAY_H
#define G__BSDECAY_H
//
#include <vector>
#include <iostream>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TRandom.h>

//
// Generate and store Bs-->Ds- pi+ (Ds- --> phi pi-) (phi --> K+ K-)

class BsDecay
 {
	//
	// 
	// Author: F. Bedeschi, INFN-Pisa, Italy
	// May 16, 2021
	//
private:
	//	
	// Variables
	Int_t fNtrk;			// particles from Bs chain decay
	std::vector<TVector3*> fPlist;	// List of track momenta in Track collections
	std::vector<Int_t>fNlist;	// List of particle types
	std::vector<Double_t>fMlist;	// List of particle masses
	std::vector<Double_t>fLlist;	// List of particle lifepath (meters)
	Int_t fBsInd;			// Index of Bs
	Int_t fPiBsInd;			// Index of pion from Bs
	Int_t fDsInd;			// Index of Ds
	Int_t fPhiInd;			// Index of Phi
	Int_t fPiDsInd;			// Index of Pi from Ds
	Int_t fKpPhiInd;		// Index of K+ from Phi
	Int_t fKmPhiInd;		// Index of K- from phi
public:
	//
	// Constructors
	BsDecay(TVector3 pBs);
	// Destructor
	~BsDecay();
	// 
	// Particle locators
	//
	// Methods
	// Plot
	void PlotTrack(TVectorD pr, Double_t phiMx, Int_t Npx, Double_t *x, Double_t *y);
	// Decay model
	TLorentzVector Boosted(TLorentzVector Booster, TLorentzVector InMoment);
	void TwoBodyPS(TLorentzVector pin, Double_t m1, Double_t m2, TLorentzVector**& pout);
	// Locators:
	Int_t GetBs(){return fBsInd;};		// Index of Bs
	Int_t GetPiBs(){return fPiBsInd;};		// Index of pion from Bs
	Int_t GetDs(){return fDsInd;};		// Index of Ds
	Int_t GetPhi(){return fPhiInd;};		// Index of Phi
	Int_t GetPiDs(){return fPiDsInd;};		// Index of Pi from Ds
	Int_t GetKpPhi(){return fKpPhiInd;};	// Index of K+ from Phi
	Int_t GetKmPhi(){return fKmPhiInd;};	// Index of K- from phi
	//
	TVector3 GetMomentum(Int_t i){return *fPlist[i];};
	Int_t GetPtype(Int_t i){return fNlist[i];};
	Double_t GetPmass(Int_t i){return fMlist[i];};
	Double_t GetPpath(Int_t i){return fLlist[i];};
};

#endif
