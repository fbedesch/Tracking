#include <TMath.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TRandom.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include "SolGridCov.h"
#include "ObsTrk.h"
#include "VertexFit.h"
//
//
Double_t fchi2(Double_t *x, Double_t *p)
{
	Double_t Ndof = p[0];
	Double_t Norm = p[1];
	Double_t z = x[0];
	//
	Double_t ex = Ndof / 2.0 - 1.0;
	Double_t chi2 = pow(z, ex)*TMath::Exp(-z / 2.0);
	Double_t bin = 5.*Ndof / 100.;
	chi2 = Norm * chi2 * bin / (pow(2, Ndof / 2.0)*TMath::Gamma(Ndof / 2.0));
	return chi2;
}
//
// Testing program
//
void TestVertex(Int_t Nvtx = 100, Int_t Ntr = 2)
{
	//
	// Init geometry
	Double_t Bfield = 2.0;	// Magnetic field
	SolGridCov *GC = new SolGridCov();
	GC->Read("CovIDEA-BASE.root");			// Read in covariance array
	//
	// Ranges
	//
	Double_t ThDegMin = 40.0;
	Double_t ThDegMax = 140.0;
	Double_t Lmin = 0.020;
	Double_t Lmax = 0.020;
	Double_t dTheta = 0.20;
	Double_t dPhi = 0.20;
	Double_t Pmin = 0.5;
	Double_t Pmax = 1.;
	//
	// Histograms
	Int_t Nbin = 100;
	//hTry = new TH1F("hTry", "Number of iterations", Nbin, 0., 100.);
	TH1D *hXpull = new TH1D("hXpull", "Pull X vertex component", Nbin, -10., 10.);
	TH1D *hYpull = new TH1D("hYpull", "Pull Y vertex component", Nbin, -10., 10.);
	TH1D *hZpull = new TH1D("hZpull", "Pull Z vertex component", Nbin, -10., 10.);
	Double_t Ndof = 2.0 * Ntr - 3.0;
	TH1D *hChi2 = new TH1D("hChi2", "Vertex #chi^{2}", 100, 0., 5 * Ndof);
	TF1 *fch = new TF1("fch", fchi2, 0., 5 * Ndof, 2);
	fch->SetParameter(0, Ndof);
	fch->SetParameter(1, Nvtx);
	//
	// Loop on # vertices
	//
	Int_t NbadX = 0.0;
	Int_t NbadY = 0.0;
	Int_t NbadZ = 0.0;
	//
	for (Int_t n = 0; n < Nvtx; n++)
	{
		//
		// Extract a direction
		Double_t ThMin = ThDegMin * TMath::Pi() / 180.;
		Double_t ThMax = ThDegMax * TMath::Pi() / 180.;
		Double_t rnTh = gRandom->Rndm();
		Double_t Th = ThMin + rnTh * (ThMax - ThMin);
		Double_t rnPh = gRandom->Rndm();
		Double_t Ph = TMath::TwoPi()*rnPh;
		//
		// Extract a flight distance (m)
		Double_t rnL = gRandom->Rndm();
		Double_t Lvtx = Lmin + rnL * (Lmax - Lmin);
		//
		TVector3 x;
		Double_t zz0 = 0.0;
		x(0) = Lvtx * TMath::Sin(Th)*TMath::Cos(Ph);
		x(1) = Lvtx * TMath::Sin(Th)*TMath::Sin(Ph);
		x(2) = Lvtx * TMath::Cos(Th) + zz0;
		//
		//cout << "Event # " << n << ", Rvertex = " << sqrt(x(0) * x(0) + x(1) * x(1)) << endl;
		//cout << "True vertex: x = " << x(0) << ", y = " << x(1) << ", z = " << x(2) << endl;
		//
		// Loop on tracks
		ObsTrk **tracks = new ObsTrk*[Ntr];	// Smear tracks according to covariance matrix
		Double_t *ppa = new Double_t[Ntr];
		Double_t *tha = new Double_t[Ntr];
		Double_t *pha = new Double_t[Ntr];
		for (Int_t i = 0; i < Ntr; i++)
		{
			Double_t rnP = gRandom->Rndm();
			Double_t Ptot = Pmin + rnP * (Pmax - Pmin);
			ppa[i] = Ptot;
			Double_t ThP = gRandom->Gaus(Th, dTheta);
			tha[i] = ThP;
			Double_t PhP = gRandom->Gaus(Ph, dPhi);
			if (PhP > TMath::TwoPi())PhP -= TMath::TwoPi();
			if (PhP < 0)PhP += TMath::TwoPi();
			pha[i] = PhP;
			//
			TVector3 P;
			P(0) = Ptot * TMath::Sin(ThP)*TMath::Cos(PhP);
			P(1) = Ptot * TMath::Sin(ThP)*TMath::Sin(PhP);
			P(2) = Ptot * TMath::Cos(ThP);
			Double_t Q = 1.0;
			if (gRandom->Rndm() > 0.5)Q = -1.0;
			if(GC->IsAccepted(P))tracks[i] = new ObsTrk(x, P, Q, Bfield, GC);
			else i--;
			//TVectorD par = tracks[i]->GetObsPar();
			//cout << "i = " << i << ", loading par = " << endl; par.Print();
		}
		//Double_t xa[3];
		//x.GetXYZ(xa);
		//cout << "xvtx = " << endl; xvtx.Print();
		//cout << "Event # " << n << endl;
		//
		// Invocation of vertexfit
		//********************************************
		// Option 1
		/*
		VertexFit* Vtx = new VertexFit(Ntr, tracks); 
		TVectorD xvtx = Vtx->GetVtx();
		TMatrixDSym covX = Vtx->GetVtxCov();
		Double_t Chi2 = Vtx->GetVtxChi2();
		*/
		// Option 2
		TVectorD** pr = new TVectorD * [Ntr];
		TMatrixDSym** cv = new TMatrixDSym * [Ntr];
		for (Int_t i = 0; i < Ntr; i++)
		{
			pr[i] = new TVectorD(tracks[i]->GetObsPar());
			cv[i] = new TMatrixDSym(tracks[i]->GetCov());
		}
		VertexFit* Vtx = new VertexFit(Ntr, pr, cv);
		TVectorD xvtx = Vtx->GetVtx();
		TMatrixDSym covX = Vtx->GetVtxCov();
		Double_t Chi2 = Vtx->GetVtxChi2();
		//********************************************
		//
		//cout << "Fit vertex: x = " << xvtx(0) << ", y = " << xvtx(1) << ", z = " << xvtx(2) << endl;
		//cout << "Cov(x):" << endl; covX.Print();
		//
		// X pulls
		//
		Double_t PullX = (xvtx(0) - x(0)) / TMath::Sqrt(covX(0, 0));
		hXpull->Fill(PullX);
		Int_t Nover  = hXpull->GetBinContent(Nbin + 1);
		if (Nover > NbadX)
		{
			cout << "Event: " << n << ", X pull overflow. PullX = " << PullX
				<< ", Xv = "<<xvtx(0)<<", cov(0,0) = " << covX(0, 0) << endl;
			NbadX = Nover;
		}
		if (PullX > 10. || PullX < -10.)cout << "Event " << n << ", PullX = " << PullX << endl;
		//
		// Y pulls
		//
		Double_t PullY = (xvtx(1) - x(1)) / TMath::Sqrt(covX(1, 1));
		hYpull->Fill(PullY);
		Nover = hYpull->GetBinContent(Nbin + 1);
		if (Nover > NbadY)
		{
			cout << "Event: " << n << ", Y pull overflow. PullY = " << PullY
				<< ", Yv = " << xvtx(1) << ", cov(1,1) = " << covX(1, 1) << endl;
			NbadY = Nover;
		}
		if (PullY > 10 || PullY < -10.)cout << "Event " << n << ", PullY = " << PullY << endl;
		//
		// Z pulls
		//
		Double_t PullZ = (xvtx(2) - x(2)) / TMath::Sqrt(covX(2, 2));
		hZpull->Fill(PullZ);
		Nover = hZpull->GetBinContent(Nbin + 1);
		if (Nover > NbadZ)
		{
			cout << "Event: " << n << ", Z pull overflow. PullZ = " << PullZ 
				<< ", Zv = " << xvtx(2) <<", cov(2,2) = "<< covX(2,2)<< endl;

			NbadZ = Nover;
		}
		if (PullZ > 10. || PullZ < -10.)cout << "Event " << n << ", PullZ = " << PullZ << endl;
		hChi2->Fill(Chi2);
		//
		// Debug
		//
		if (TMath::Abs(PullX) > 10. || TMath::Abs(PullY) > 10. || TMath::Abs(PullZ) > 10.)
		{
			cout << "Event number: " << n << endl;
			cout << "True vertex: x = " << x(0) << ", y = " << x(1) << ", z = " << x(2) << endl;
			//TVectorD x00(3);
			//x00 = Vertex0(Ntr, tracks);
			//cout << "Fit vertex0: x = " << x00(0) << ", y = " << x00(1) << ", z = " << x00(2) << endl;
			cout << "Fit vertex: x = " << xvtx(0) << ", y = " << xvtx(1) << ", z = " << xvtx(2)
				<< ", Chi square = " << Chi2 << endl;
			cout << "Pulls: x = " << PullX << ", y = " << PullY << ", z = " << PullZ << endl;
			cout << "Cov matrix: " << endl; covX.Print();
			for (Int_t i = 0; i < Ntr; i++)
			{
				cout << "track # " << i;
				cout << ", Ptot = " << ppa[i] << ", Theta = " << tha[i] << ", Phi = " << pha[i] << endl;
			}
		}
		//
		// Cleanup
		//
		for (Int_t i = 0; i < Ntr; i++) delete tracks[i];
		for (Int_t i = 0; i < Ntr; i++) delete pr[i];
		for (Int_t i = 0; i < Ntr; i++) delete cv[i];
		delete[] tracks;
		delete[] pr;
		delete[] cv;
		delete[] ppa;
		delete[] tha;
		delete[] pha;
	}
	//
	// Plots
	TCanvas *cnv = new TCanvas("cnv", "Normalized residuals", 20, 20, 500, 500);
	cnv->Divide(2, 2);
	cnv->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hXpull->Fit("gaus");
	hXpull->Draw();
	cnv->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hYpull->Fit("gaus");
	hYpull->Draw();
	cnv->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hZpull->Fit("gaus");
	hZpull->Draw();
	cnv->cd(4);
	hChi2->Fit("fch");
	hChi2->Draw();
	//
	
}

