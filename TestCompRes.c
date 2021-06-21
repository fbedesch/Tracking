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
void TestCompRes(Double_t Ang)
{
	//
	// Track to draw
	Double_t x[3] = { 0.0, 0.0, 0.0 };		// Track starting point
	Double_t ppt = 10.;						// Track pt
	Double_t th = Ang * TMath::Pi() / 180.;
	Double_t ppz = ppt / TMath::Tan(th);	// Track pz
	Double_t p[3] = { ppt, 0.0, ppz };		// Track momentum
	//
	//**************************
	//	Initialize geometry    *
	//**************************
	//
	//
	SolGeom *Gidea;			// Initialize IDEA geometry
	Gidea = new SolGeom("GeoIDEA_BASE.txt");	// Geometry IDEA
	Gidea->Draw();			// Draw IDEA geometry
	TCanvas *cc_id = Gidea->cnv();					// Get canvas with geo display
	cc_id->cd(1);
	SolTrack *trk_id = new SolTrack(x, p, Gidea);	// Initialize track
	TGraph *gr_id = trk_id->TrkPlot();			// graph intersection with layers
	gr_id->Draw("PLSAME");						// plot track
	TCanvas *cnv_id = new TCanvas("cnv_id", "IDEA Geometry");
	cnv_id->Divide(1, 1); cnv_id->cd(1); cc_id->DrawClonePad();
	//
	SolGeom *Gcld;			// Initialize CLD  geometry	
	Gcld  = new SolGeom("GeoCLD.txt");			// Geometry CLD
	Gcld->Draw();			// Draw CLD  geometry
	TCanvas *cc_cl = Gcld->cnv();
	cc_cl->cd(1);
	SolTrack *trk_cl = new SolTrack(x, p, Gcld);	// Initialize track
	TGraph *gr_cl = trk_cl->TrkPlot();			// graph intersection with layers
	gr_cl->Draw("PLSAME");						// plot track
	TCanvas *cnv_cl = new TCanvas("cnv_cl", "CLD Geometry");
	cnv_cl->Divide(1, 1); cnv_cl->cd(1); cc_cl->DrawClonePad();
	cc_cl->Close(); gSystem->ProcessEvents(); delete cc_cl; cc_cl = 0;
	//
	//
	// Geometry initialized and drawn
	//********************************************************
	//	
	//******************************************************
	// Compare track parameter resolutions vs pt and theta *
	//******************************************************
	//
	TCanvas *resol = new TCanvas("resol", "Comparison of resolutions", 100, 100, 500, 500);
	resol->Divide(2, 2);
	// Define graphs for IDEA
	TGraph *grpt_id;				// pt resolution graphs
	TGraph *grptms_id;				// pt resolution graphs MS only
	TGraph *grd0_id;				// D resolution graphs
	TGraph *grz0_id;				// z0 resolution graphs
	TGraph *grth_id;				// theta resolution
	// Define graphs for CLD
	TGraph *grpt_cl;				// pt resolution graphs
	TGraph *grptms_cl;				// pt resolution graphs MS only
	TGraph *grd0_cl;				// D resolution graphs
	TGraph *grz0_cl;				// z0 resolution graphs
	TGraph *grth_cl;				// theta resolution
	// Setup graph arrays
	Int_t Npt = 200;			// Nr. of points per graph
	Double_t * pt = new Double_t[Npt];
	Double_t * pp = new Double_t[Npt];
	Double_t *spt_id = new Double_t[Npt];
	Double_t *sptms_id = new Double_t[Npt];
	Double_t *sd0_id = new Double_t[Npt];
	Double_t *sz0_id = new Double_t[Npt];
	Double_t *sth_id = new Double_t[Npt];
	Double_t *spt_cl = new Double_t[Npt];
	Double_t *sptms_cl = new Double_t[Npt];
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
		pt[k] = ptmin+k*pts;					// Set transverse momentum
		p[0] = pt[k]; p[1] = 0;	
		pp[k] = pt[k]/TMath::Sin(th);			// Set momentum
		p[2] = TMath::Sqrt(pp[k] * pp[k] - p[0] * p[0]);
		// Fill IDEA arrays
		SolTrack *tr_id = new SolTrack(x, p, Gidea);	// Initialize track
		Bool_t Res = kTRUE;		// Enable detector resolution effects
		Bool_t MS  = kTRUE;		// Enable multiple scattering
		tr_id->CovCalc(Res,MS);					// Calculate covariance
		spt_id[k] = tr_id->s_pt();							// Dpt/pt
		sd0_id[k] = tr_id->s_D()*1e6;							// D  res. - change to microns
		sz0_id[k] = tr_id->s_z0()*1e6;						// z0 res. - change to microns
		sth_id[k] = tr_id->s_ct() / (1 + pow(tr_id->ct(), 2));	// theta resolution
		//
		Res = kFALSE;
		MS = kTRUE;
		tr_id->CovCalc(Res,MS);					// Calculate covariance with only MS
		sptms_id[k] = tr_id->s_pt();							// Dpt/pt
		// Fill CLD arrays
		SolTrack *tr_cl = new SolTrack(x, p, Gcld);	// Initialize track
		Res = kTRUE;		// Enable detector resolution effects
		MS  = kTRUE;		// Enable multiple scattering
		tr_cl->CovCalc(Res,MS);					// Calculate covariance
		spt_cl[k] = tr_cl->s_pt();							// Dpt/pt
		sd0_cl[k] = tr_cl->s_D()*1e6;							// D  res. - change to microns
		sz0_cl[k] = tr_cl->s_z0()*1e6;						// z0 res. - change to microns
		sth_cl[k] = tr_cl->s_ct() / (1 + pow(tr_cl->ct(), 2));	// theta resolution
		Res = kFALSE;
		tr_cl->CovCalc(Res,MS);					// Calculate covariance with only MS
		sptms_cl[k] = tr_cl->s_pt();							// Dpt/pt
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
	grptms_cl = new TGraph(Npt, pt, sptms_cl);			// pt resolution MS only
	grptms_cl->SetLineColor(kRed);
	grptms_cl->SetMarkerColor(kRed);
	grptms_cl->SetLineStyle(7);
	grptms_cl->SetTitle("#sigma_{pt}/pt");
	grptms_cl->SetMinimum(0.0);
	grptms_cl->GetXaxis()->SetTitle("pt (GeV)");
	grptms_cl->Draw("SAME");
	grpt_id = new TGraph(Npt, pt, spt_id);			// pt resolution
	grpt_id->SetLineColor(kBlue);
	grpt_id->SetMarkerColor(kBlue);
	grpt_id->SetTitle("#sigma_{pt}/pt");
	grpt_id->SetMinimum(0.0);
	grpt_id->GetXaxis()->SetTitle("pt (GeV)");
	grpt_id->Draw("SAME");
	grptms_id = new TGraph(Npt, pt, sptms_id);			// pt resolution MS only
	grptms_id->SetLineColor(kBlue);
	grptms_id->SetMarkerColor(kBlue);
	grptms_id->SetLineStyle(7);
	grptms_id->SetTitle("#sigma_{pt}/pt");
	grptms_id->SetMinimum(0.0);
	grptms_id->GetXaxis()->SetTitle("pt (GeV)");
	grptms_id->Draw("SAME");
	Int_t iang = TMath::Nint(Ang);
	TLegend *lgpt = new TLegend(0.2, 0.9, 0.6, 0.70);
	TString LgTitle; 
	LgTitle.Form("Track angle %d deg.",iang);
	lgpt->SetHeader(LgTitle);
	lgpt->AddEntry(grpt_id, "IDEA", "L");
	lgpt->AddEntry(grptms_id, "IDEA MS only", "L");
	lgpt->AddEntry(grpt_cl, "CLD", "L");
	lgpt->AddEntry(grptms_cl, "CLD MS only", "L");
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

