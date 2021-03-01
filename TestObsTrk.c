#include <TMath.h>
#include <TH1.h>
#include <TCanvas.h> 
#include <iostream>
#include "SolGridCov.h"
#include "ObsTrk.h"
//
void TestObsTrk(Double_t ph0 = 0.0, Double_t Q = 1.0, Int_t ind = 1)
{
	//
	// Test track conversion
	//
	// Tracks from origin
	//
	Double_t txv[3] = { -0.20, 0.03, 0.10 };
	const Int_t Npt = 5;
	Double_t tpt[Npt] = { 0.5, 1.0, 5.0, 10., 100. };
	const Int_t Nct = 3;
	Double_t tct[Nct] = { 0.0, 1.0, 2.0 };
	TVector3 X(txv);
	TVector3 P;
	Double_t B = 2.0;
	SolGridCov *GC = new SolGridCov();
	GC->Read("CovIDEA-BASE.root");
	const Int_t Npar = 5;
	//
	for (Int_t ip = 0; ip < Npt; ip++)
	{
		for (Int_t ic = 0; ic < Nct; ic++)
		{
			P(0) = tpt[ip] * TMath::Cos(ph0);
			P(1) = tpt[ip] * TMath::Sin(ph0);
			P(2) = tct[ic] * tpt[ip];
			ObsTrk *Trk = new ObsTrk(X, P, Q, B, GC);
			TVectorD gPar = Trk->GetGenPar();
			for (Int_t i = 0; i < Npar; i++)cout << "  gPar(" << i << ") = " << gPar(i);
			cout << endl;
			TVector3 gX = Trk->GetGenX();
			TVector3 gP = Trk->GetGenP();
			Double_t gQ = Trk->GetGenQ();
			TVector3 gpX = Trk->ParToX(gPar);
			TVector3 gpP = Trk->ParToP(gPar, B);
			Double_t gpQ = Trk->ParToQ(gPar);
			cout << "In  trk:";
			for (Int_t i = 0; i < 3; i++)cout << " x(" << i << ") = " << gX(i);
			for (Int_t i = 0; i < 3; i++)cout << " p(" << i << ") = " << gP(i);
			cout << " Q = "<<gQ<<endl;
			cout << "Out trk:";
			for (Int_t i = 0; i < 3; i++)cout << " x(" << i << ") = " << gpX(i);
			for (Int_t i = 0; i < 3; i++)cout << " p(" << i << ") = " << gpP(i);
			cout << " Q = " << gpQ << endl;
		}
	}
	//
	// Test Randomization
	//
	// Set track parameters
	P(0) = tpt[ind] * TMath::Cos(ph0);
	P(1) = tpt[ind] * TMath::Sin(ph0);
	P(2) = tct[1] * tpt[ind];
	// 
	// Init histograms
	TH1D *h_D   = new TH1D("h_D", "Impact parameter resolution", 100, -5., 5.);
	TH1D *h_ph0 = new TH1D("h_ph0", "Phi0 resolution"          , 100, -5., 5.);
	TH1D *h_C   = new TH1D("h_C", "Curvature resolution"       , 100, -5., 5.);
	TH1D *h_z0  = new TH1D("h_z0", "z0 resolution"             , 100, -5., 5.);
	TH1D *h_ct  = new TH1D("h_ct", "Cot(theta) resolution"     , 100, -5., 5.);
	// Fill histograms
	Int_t Nev = 10000;
	TMatrixDSym Cri(5);
	TMatrixDSym Cv(5); 
	for (Int_t n = 0; n < Nev; n++)
	{
		ObsTrk *Trk = new ObsTrk(X, P, Q, B, GC);
		TVectorD gPar = Trk->GetGenPar();
		TVectorD oPar = Trk->GetObsPar();
		Cv = Trk->GetCov();
		h_D  ->Fill((gPar(0) - oPar(0)) / TMath::Sqrt(Cv(0, 0)));
		h_ph0->Fill((gPar(1) - oPar(1)) / TMath::Sqrt(Cv(1, 1)));
		h_C  ->Fill((gPar(2) - oPar(2)) / TMath::Sqrt(Cv(2, 2)));
		h_z0 ->Fill((gPar(3) - oPar(3)) / TMath::Sqrt(Cv(3, 3)));
		h_ct ->Fill((gPar(4) - oPar(4)) / TMath::Sqrt(Cv(4, 4)));
		//
		// Check correlations
		TVectorD pDiff = oPar - gPar;
		Cri.Rank1Update(pDiff);
	}
	TMatrixDSym Dinv(5); Dinv.Zero();
	for (Int_t i = 0; i < 5; i++)Dinv(i, i) = 1.0 / TMath::Sqrt(Cv(i, i));
	Cri = (1.0/(Double_t)Nev)*Cri;
	cout << "Covariance matrix:" << endl;
	Cv.Similarity(Dinv);
	Cv.Print();
	cout << "Correlation found:" << endl;
	Cri.Similarity(Dinv);
	Cri.Print();
	// Plot histograms
	TCanvas *cc = new TCanvas("cc", "Normalized resolutions", 50, 50, 800, 500);
	cc->Divide(2, 3);
	cc->cd(1);
	h_D->Fit("gaus");
	h_D->Draw();
	cc->cd(2);
	h_ph0->Fit("gaus");
	h_ph0->Draw();
	cc->cd(3);
	h_C->Fit("gaus");
	h_C->Draw();
	cc->cd(4);
	h_z0->Fit("gaus");
	h_z0->Draw();
	cc->cd(5);
	h_ct->Fit("gaus");
	h_ct->Draw();
}
