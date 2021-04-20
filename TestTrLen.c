#include <iostream>
#include <TMath.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TF1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLine.h>
#include <TRandom.h>
#include <TCanvas.h>
#include "TrkUtil.h"

Double_t fRz(Double_t* x, Double_t* Par)
{
	// Track pararameters
	Double_t D = Par[0];		// Transverse impact parameter
	Double_t phi0 = Par[1];		// Transverse direction at minimum approach
	Double_t C = Par[2];		// Half curvature
	Double_t z0 = Par[3];		// Z at minimum approach
	Double_t ct = Par[4];		// cot(theta)
	//
	Double_t z = x[0];
	//
	Double_t R = 0;
	if ((z - z0) / ct > 0.0)
	{
		Double_t ph = 2 * C * (z - z0) / ct;
		//
		Double_t x0 = -D * TMath::Sin(phi0);
		Double_t y0 = D * TMath::Cos(phi0);
		Double_t xtr = x0 + (TMath::Sin(ph + phi0) - TMath::Sin(phi0)) / (2 * C);
		Double_t ytr = y0 - (TMath::Cos(ph + phi0) - TMath::Cos(phi0)) / (2 * C);
		//
		R = TMath::Sqrt(xtr * xtr + ytr * ytr);
	}
	//
	return R;
}
void TestTrLen(Double_t ang = 90., Double_t z0 = 0.0)
{
	//
	// Constants
	Double_t angRad = ang*TMath::Pi() / 180.;	// To radians
	Double_t Bz = 2.0;				// Field
	TVector3 x(0.0, 0.0, z0);		// Track origin
	//
	// Setup chamber
	//
	TrkUtil *TU = new TrkUtil(Bz);
	//cout << "TU done" << endl;
	Double_t Rmin = 0.35;
	Double_t Rmax = 2.0;
	Double_t Zmin = -2.0;
	Double_t Zmax = 2.0;
	TU->SetDchBoundaries(Rmin, Rmax, Zmin, Zmax);
	//
	// Functions
	//
	const Int_t Ntr = 9;
	Double_t ptr[Ntr] = { 0.2,0.3,0.4,0.5,0.6,0.8,1.0,10.,50.};
	Double_t tln[Ntr];
	TF1* fTrkRz[Ntr];
	for (Int_t i = 0; i < Ntr; i++)
	{
		char fName[10]; 
		Int_t status = sprintf(fName, "fRz%d",i);
		fTrkRz  [i] = new TF1(fName, fRz  , 2*Zmin, 2*Zmax, 5);	// R vs z
		Double_t px = ptr[i] * TMath::Sin(angRad);
		Double_t pz = ptr[i] * TMath::Cos(angRad);
		TVector3 p(px, 0.0, pz);
		Double_t Q = 1.0;
		TVectorD Par = TrkUtil::XPtoPar(x, p, Q, Bz);
		Double_t* pr = new Double_t[5];;
		pr = Par.GetMatrixArray();
		fTrkRz[i]->SetParameters(pr);
		fTrkRz[i]->SetMinimum(0.0);
		fTrkRz[i]->SetMaximum(Rmax+0.1);
		fTrkRz[i]->SetLineColor(i+1);
		Double_t len = TU->TrkLen(Par);
		std::cout << "Track momentum = " << ptr[i] << " --> Length = " << len << std::endl;
		tln[i] = len;
	}
	//
	// Plots
	//
	TCanvas *cnv = new TCanvas("cnv", "Track length display", 50, 50, 800, 500);
	cnv->Divide(1, 1);
	//
	// Draw trajectories
	fTrkRz[0]->Draw();
	for (Int_t i = 1; i < Ntr; i++) fTrkRz[i]->Draw("SAME");
	TGraph* gTrkl = new TGraph(Ntr, ptr, tln);
	//
	// Draw R-z box
	TLine* lRin  = new TLine(Zmin, Rmin, Zmax, Rmin);
	TLine *lRout = new TLine(Zmin, Rmax, Zmax, Rmax);
	TLine* lZin  = new TLine(Zmin, Rmin, Zmin, Rmax);
	TLine* lZout = new TLine(Zmax, Rmin, Zmax, Rmax);
	lRin ->Draw("SAME");
	lRout->Draw("SAME");
	lZin ->Draw("SAME");
	lZout->Draw("SAME");
	//
	TCanvas* cnv1 = new TCanvas("cnv1", "Track length plot", 100, 100, 800, 500);
	cnv1->Divide(1, 1);
	cnv1->cd(1);
	gPad->SetLogx();
	gTrkl->SetTitle("Track length vs. p");
	gTrkl->GetXaxis()->SetTitle("Momentum (GeV)");
	gTrkl->GetYaxis()->SetTitle("Track length (meter)");
	gTrkl->Draw("APL");
}