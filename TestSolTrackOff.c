#include <TMath.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h> 
#include <iostream>
#include "SolGeom.h"
#include "TrkUtil.h"
#include "SolTrack.h"
#include "ObsTrk.h"
//
void TestSolTrackOff(Double_t theta = 90)
{
	//
	//
	// Load geometry
	//
	SolGeom* G = new SolGeom("GeoIDEA_BASE.txt");
	std::cout << "Min radius = " << G->GetRmin()<<", ";
	std::cout << "Negative Z min = " << G->GetZminNeg()<<", ";
	std::cout << "Positive Z min = " << G->GetZminPos() << std::endl;
	//
	// Setup arrays
	const Int_t Npt = 17;
	Double_t dlen[Npt] = { 0.0, 1.0e-2,2.0e-2,3.0e-2,5.0e-2,8.e-2,10.e-2,20.e-2,
		40.e-2,50.e-2,60.e-2,70.e-2,100.e-2,130.e-2,150.e-2,200.e-2,250.e-2};
	Double_t Nhits[Npt];
	Double_t NmHits[Npt];
	Double_t Rh[Npt];
	Double_t Zh[Npt];
	Double_t sD0[Npt];
	Double_t sPt[Npt];
	Double_t th = TMath::Pi() * theta / 180.;
	Double_t pTot = 100.;
	TVector3 p(0.0, pTot * TMath::Sin(th), pTot * TMath::Cos(th));
	//
	Int_t Nct = 0;
	for (Int_t ip = 0; ip < Npt; ip++)
	{
		Rh[ip] = dlen[ip] * TMath::Sin(th);
		Zh[ip] = dlen[ip] * TMath::Cos(th);
		TVector3 x(0.0,Rh[ip], Zh[ip]);
		SolTrack* trk = new SolTrack(x, p, G);
		Nhits[ip] = (Double_t)trk->nHit();
		NmHits[ip] = (Double_t)trk->nmHit();
		if (NmHits[ip] > 6) {
			Nct++;
			Bool_t Res = kTRUE;
			Bool_t MS = kTRUE;
			trk->CovCalc(Res, MS);
			sD0[ip] = trk->s_D() * 1.e6;
			sPt[ip] = trk->s_C() / TMath::Abs(trk->C());
		}
	}
	TCanvas *cHit = new TCanvas("cHit", "Hit count", 50, 50, 800, 500);
	cHit->Divide(2, 2);
	cHit->cd(1);
	TGraph* gHit = new TGraph(Nct, Rh, Nhits);	// All hits
	gHit->SetTitle("Total nr. hits");
	gHit->Draw("APL");
	cHit->cd(2);				
	TGraph* gHitM = new TGraph(Nct, Rh, NmHits);	// Measure hits
	gHitM->SetTitle("Total nr. measure hits");
	gHitM->Draw("APL");
	cHit->cd(3);
	TGraph* gHitSd = new TGraph(Nct, Rh, sD0);	// Measure hits
	gHitSd->SetTitle("D error");
	gHitSd->Draw("APL");
	cHit->cd(4);
	TGraph* gHitSpt = new TGraph(Nct, Rh, sPt);	// Measure hits
	gHitSpt->SetTitle("Pt relative error");
	gHitSpt->Draw("APL");
}
