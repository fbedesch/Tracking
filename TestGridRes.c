//
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>
#include <THStack.h>
#include <TFile.h>
#include <TSystem.h>
#include <iostream>
#include "SolGeom.h"
#include "SolTrack.h"
#include "SolGridCov.h"
//
void TestGridRes(Double_t Perp, Double_t Ang)
{
	//
	// Track to draw
	Double_t x[3] = { 0.0, 0.0, 0.0 };		// Track starting point
	Double_t ppt = 10.;				// Track pt
	Double_t th = Ang * TMath::Pi() / 180.;
	Double_t ppz = ppt / TMath::Tan(th);		// Track pz
	Double_t p[3] = { ppt, 0.0, ppz };		// Track momentum
	//
	//**************************
	//	Initialize geometry    *
	//**************************
	//
	//
	SolGeom* Gid = new SolGeom("GeoIDEA_BASE.txt");
	SolGridCov *Gidea = new SolGridCov();	// Initialize IDEA geometry
	Gidea->Read("CovIDEA_BASE.root");	// Geometry IDEA
	//
	SolGridCov *Gcld = new SolGridCov();	// Initialize CLD  geometry	
	Gcld->Read("CovCLD.root");		// Geometry CLD
	//
	//******************************************************
	// Compare track parameter resolutions vs pt           *
	//******************************************************
	//
	TCanvas *resol = new TCanvas("resol", "Comparison of resolutions vs pt ", 100, 100, 500, 500);
	resol->Divide(2, 2);
	// Define graphs for IDEA
	TGraph *grpt_id;				// pt resolution graphs
	TGraph *grd0_id;				// D resolution graphs
	TGraph *grz0_id;				// z0 resolution graphs
	TGraph *grth_id;				// theta resolution
	// Define graphs for CLD
	TGraph *grpt_cl;				// pt resolution graphs
	TGraph *grd0_cl;				// D resolution graphs
	TGraph *grz0_cl;				// z0 resolution graphs
	TGraph *grth_cl;				// theta resolution
	// Setup graph arrays
	Int_t Npt = 200;				// Nr. of points per graph
	Double_t * pt = new Double_t[Npt];
	Double_t * pp = new Double_t[Npt];
	Double_t *spt_id = new Double_t[Npt];
	Double_t *sd0_id = new Double_t[Npt];
	Double_t *sz0_id = new Double_t[Npt];
	Double_t *sth_id = new Double_t[Npt];
	Double_t *spt_cl = new Double_t[Npt];
	Double_t *sd0_cl = new Double_t[Npt];
	Double_t *sz0_cl = new Double_t[Npt];
	Double_t *sth_cl = new Double_t[Npt];
	// Fill graph arrays
	Double_t ptmin = 2.0;
	Double_t ptmax = 100;
	Double_t pts = (ptmax - ptmin) / (Double_t)(Npt-1);
	for (Int_t k = 0; k < Npt; k++)	// Loop on pt
	{
		Double_t x[3]; Double_t p[3];
		x[0] = 0; x[1] = 0; x[2] = 0;			// Set origin
		pt[k] = ptmin+k*pts;				// Set transverse momentum
		p[0] = pt[k]; p[1] = 0;	
		pp[k] = pt[k]/TMath::Sin(th);			// Set momentum
		p[2] = TMath::Sqrt(pp[k] * pp[k] - p[0] * p[0]);
		// Fill IDEA arrays
		Double_t theta = TMath::ATan2(pt[k],p[2])*180./TMath::Pi();
		TMatrixDSym cov_id = Gidea->GetCov(pt[k],theta);
		Double_t sD_id  = TMath::Sqrt(cov_id(0,0));
		Double_t sPh_id = TMath::Sqrt(cov_id(1,1));
		Double_t sC_id  = TMath::Sqrt(cov_id(2,2));
		Double_t sPt_id = 2 * sC_id*pt[k] / (0.2998*Gid->B());
		Double_t sZ_id  = TMath::Sqrt(cov_id(3,3));
		Double_t sCt_id = TMath::Sqrt(cov_id(4,4));
		spt_id[k] = sPt_id;				// Dpt/pt
		sd0_id[k] = sD_id*1e6;				// D  res. - change to microns
		sz0_id[k] = sZ_id*1e6;				// z0 res. - change to microns
		sth_id[k] = sCt_id / (1 + pow(p[2]/pt[k], 2));	// theta resolution
		//
		// Fill CLD arrays
		TMatrixDSym cov_cl = Gcld->GetCov(pt[k],theta);
		Double_t sD_cl  = TMath::Sqrt(cov_cl(0,0));
		Double_t sPh_cl = TMath::Sqrt(cov_cl(1,1));
		Double_t sC_cl  = TMath::Sqrt(cov_cl(2,2));
		Double_t sPt_cl = 2 * sC_cl*pt[k] / (0.2998*Gid->B());
		Double_t sZ_cl  = TMath::Sqrt(cov_cl(3,3));
		Double_t sCt_cl = TMath::Sqrt(cov_cl(4,4));
		spt_cl[k] = sPt_cl;				// Dpt/pt
		sd0_cl[k] = sD_cl*1e6;				// D  res. - change to microns
		sz0_cl[k] = sZ_cl*1e6;				// z0 res. - change to microns
		sth_cl[k] = sCt_cl / (1 + pow(p[2]/pt[k], 2));	// theta resolution
	}
	//
	// Compare pt resolution
	resol->cd(1);
	grpt_cl = new TGraph(Npt, pt, spt_cl);			// pt resolution
	grpt_cl->SetLineColor(kRed);
	grpt_cl->SetMarkerColor(kRed);
	grpt_cl->SetTitle("#sigma_{pt}/pt");
	grpt_cl->SetMinimum(0.0);
	grpt_cl->GetXaxis()->SetTitle("pt (GeV)");
	grpt_cl->Draw("APL");
	grpt_id = new TGraph(Npt, pt, spt_id);			// pt resolution
	grpt_id->SetLineColor(kBlue);
	grpt_id->SetMarkerColor(kBlue);
	grpt_id->SetTitle("#sigma_{pt}/pt");
	grpt_id->SetMinimum(0.0);
	grpt_id->GetXaxis()->SetTitle("pt (GeV)");
	grpt_id->Draw("SAME");
	Int_t iang = TMath::Nint(Ang);
	TLegend *lgpt = new TLegend(0.2, 0.9, 0.6, 0.70);
	TString LgTitle; 
	LgTitle.Form("Track angle %d deg.",iang);
	lgpt->SetHeader(LgTitle);
	lgpt->AddEntry(grpt_id, "IDEA", "L");
	lgpt->AddEntry(grpt_cl, "CLD", "L");
	lgpt->Draw();
	// Compare d0 resolution
	resol->cd(2);
	grd0_id = new TGraph(Npt, pp, sd0_id);			// D resolution
	grd0_id->SetLineColor(kBlue);
	grd0_id->SetMarkerColor(kBlue);
	grd0_id->SetTitle("D_{0} (#mum)");
	grd0_id->SetMinimum(0.0);
	grd0_id->GetXaxis()->SetTitle("p (GeV)");
	grd0_id->Draw("APL");
	grd0_cl = new TGraph(Npt, pp, sd0_cl);			// D resolution
	grd0_cl->SetLineColor(kRed);
	grd0_cl->SetMarkerColor(kRed);
	grd0_cl->SetTitle("D_{0} (#mum)");
	grd0_cl->SetMinimum(0.0);
	grd0_cl->GetXaxis()->SetTitle("p (GeV)");
	grd0_cl->Draw("SAME"); 
	TLegend *lgd0 = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgd0->SetHeader(LgTitle);
	lgd0->AddEntry(grpt_id, "IDEA", "L");
	lgd0->AddEntry(grpt_cl, "CLD", "L");
	lgd0->Draw();
	// Compare z0 resolution
	resol->cd(3);
	grz0_id = new TGraph(Npt, pp, sz0_id);			// z0 resolution
	grz0_id->SetLineColor(kBlue);
	grz0_id->SetMarkerColor(kBlue);
	grz0_id->SetTitle("Z_{0} (#mum)");
	grz0_id->GetXaxis()->SetTitle("p (GeV)");
	grz0_id->Draw("APL");
	grz0_cl = new TGraph(Npt, pp, sz0_cl);			// z0 resolution
	grz0_cl->SetLineColor(kRed);
	grz0_cl->SetMarkerColor(kRed);
	grz0_cl->SetTitle("Z_{0} (#mum)");
	grz0_cl->GetXaxis()->SetTitle("p (GeV)");
	grz0_cl->Draw("SAME");			// Compare theta resolution
	TLegend *lgz0 = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgz0->SetHeader(LgTitle);
	lgz0->AddEntry(grpt_id, "IDEA", "L");
	lgz0->AddEntry(grpt_cl, "CLD", "L");
	lgz0->Draw();
	resol->cd(4);
	grth_id = new TGraph(Npt, pp, sth_id);			// theta resolution
	grth_id->SetLineColor(kBlue);
	grth_id->SetMarkerColor(kBlue);
	grth_id->SetTitle("#theta (rad)");
	grth_id->SetMinimum(0.0);
	grth_id->GetXaxis()->SetTitle("p (GeV)");
	grth_id->Draw("APL");
	grth_cl = new TGraph(Npt, pp, sth_cl);			// theta resolution
	grth_cl->SetLineColor(kRed);
	grth_cl->SetMarkerColor(kRed);
	grth_cl->SetTitle("#theta (rad)");
	grth_cl->SetMinimum(0.0);
	grth_cl->GetXaxis()->SetTitle("p (GeV)");
	grth_cl->Draw("SAME");
	TLegend *lgth = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgth->SetHeader(LgTitle);
	lgth->AddEntry(grpt_id, "IDEA", "L");
	lgth->AddEntry(grpt_cl, "CLD", "L");
	lgth->Draw();
	//
}

