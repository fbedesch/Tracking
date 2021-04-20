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
void TestVertexEC(Int_t Nvtx = 100, Int_t Ntr = 2)
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
	Double_t Lmin = 0.002;
	Double_t Lmax = 0.004;
	Double_t dTheta = 0.20;
	Double_t dPhi = 0.20;
	Double_t Pmin = 0.5;
	Double_t Pmax = 1.;
	//
	// Constraint
	//
	TMatrixDSym V(3); V.Zero();
	V(0,0) = 1.e-5*1.e-5;	// 10 microns
	V(1,1) = 1.e-6*1.e-6;	// 1 microns
	V(2,2) = 5.e-4*5.e-4;	// 500 microns
	//
	// Histograms
	Int_t Nbin = 100;
	//hTry = new TH1F("hTry", "Number of iterations", Nbin, 0., 100.);
	TH1D *hXpull = new TH1D("hXpull", "Pull X vertex component", Nbin, -10., 10.);
	TH1D *hYpull = new TH1D("hYpull", "Pull Y vertex component", Nbin, -10., 10.);
	TH1D *hZpull = new TH1D("hZpull", "Pull Z vertex component", Nbin, -10., 10.);
	//Double_t Ndof = 2.0 * Ntr - 3.0;
	Double_t Ndof = 2.0 * Ntr;
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
		Double_t zz0 = 0.01;
		TVectorD xc(3);
		xc(0) = Lvtx * TMath::Sin(Th)*TMath::Cos(Ph);
		x(0) = gRandom->Gaus(xc(0), TMath::Sqrt(V(0,0)));
		xc(1) = Lvtx * TMath::Sin(Th)*TMath::Sin(Ph);
		x(1) = gRandom->Gaus(xc(1), TMath::Sqrt(V(1,1)));
		xc(2) = Lvtx * TMath::Cos(Th) + zz0;
		x(2) = gRandom->Gaus(xc(2), TMath::Sqrt(V(2,2)));
		//
		if(n%5000 == 0)cout << "Event # " << n << ", Rvertex = " << sqrt(x(0) * x(0) + x(1) * x(1)) << endl;
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
		//********************************************
		// Invocation of vertexfit    ****************
		//********************************************
		// Option 1
		//std::cout << "TestVertex: before constructor" << std::endl;
		
		VertexFit* Vtx = new VertexFit(Ntr, tracks);
		//std::cout << "TestVertex: after constructor" << std::endl;
		Vtx->AddVtxConstraint(xc, V);
		TVectorD xvtx = Vtx->GetVtx();
		//std::cout << "TestVertex: after GetVtx()" << std::endl;
		TMatrixDSym covX = Vtx->GetVtxCov();
		Double_t Chi2 = Vtx->GetVtxChi2();
		//
		// Option 2
		/*
		TVectorD** pr = new TVectorD * [Ntr];
		TMatrixDSym** cv = new TMatrixDSym * [Ntr];
		for (Int_t i = 0; i < Ntr; i++)
		{
			pr[i] = new TVectorD(tracks[i]->GetObsPar());
			cv[i] = new TMatrixDSym(tracks[i]->GetCov());
		}
		//std::cout << "TestVertex: before call to VertexFit" << std::endl;
		VertexFit* Vtx = new VertexFit(Ntr, pr, cv);
		//std::cout << "TestVertex: after call to VertexFit" << std::endl;
		TVectorD xvtx = Vtx->GetVtx();
		//std::cout << "TestVertex: after call to GetVtx" << std::endl;
		TMatrixDSym covX = Vtx->GetVtxCov();
		//std::cout << "TestVertex: after call to GetVtxCov" << std::endl;
		Double_t Chi2 = Vtx->GetVtxChi2();
		//std::cout << "TestVertex: after call to GetVtxChi2" << std::endl;
		*/
		//
		// Option 3
		/*
		VertexFit* Vtx = new VertexFit();
		for (Int_t i = 0; i < Ntr; i++)
		{
			TVectorD *par = new TVectorD(tracks[i]->GetObsPar());
			TMatrixDSym *Cov = new TMatrixDSym(tracks[i]->GetCov());
			Vtx->AddTrk(par, Cov);
		}
		TVectorD xvtx = Vtx->GetVtx();
		TMatrixDSym covX = Vtx->GetVtxCov();
		Double_t Chi2 = Vtx->GetVtxChi2();
		*/
		//
		// Mixed option 2 and 3
		// Add two tracks at the end
		/*
		TVectorD** pr = new TVectorD * [Ntr-2];
		TMatrixDSym** cv = new TMatrixDSym * [Ntr-2];
		for (Int_t i = 0; i < Ntr-2; i++)
		{
			pr[i] = new TVectorD(tracks[i]->GetObsPar());
			cv[i] = new TMatrixDSym(tracks[i]->GetCov());
		}
		VertexFit* Vtx = new VertexFit(Ntr - 2, pr, cv);
		//std::cout << "TestVertex: after new Vtx" << std::endl;
		TVectorD xvtx = Vtx->GetVtx();		// Ntr-2 fit
		//std::cout << "TestVertex: first fit completed" << std::endl;
		for (Int_t i = Ntr-2; i < Ntr; i++)
		{
			TVectorD* par = new TVectorD(tracks[i]->GetObsPar());
			TMatrixDSym* Cov = new TMatrixDSym(tracks[i]->GetCov());
			Vtx->AddTrk(par, Cov);
		}
		xvtx = Vtx->GetVtx();					// Updated fit
		//std::cout << "TestVertex: updated fit completed" << std::endl;
		TMatrixDSym covX = Vtx->GetVtxCov();
		Double_t Chi2 = Vtx->GetVtxChi2();
		*/
		//
		// Option 4 
		// Track removal
		// Remove two random tracks
		/*
		TVectorD** pr = new TVectorD *[Ntr];
		TMatrixDSym** cv = new TMatrixDSym *[Ntr];
		for (Int_t i = 0; i < Ntr; i++)
		{
			pr[i] = new TVectorD(tracks[i]->GetObsPar());
			cv[i] = new TMatrixDSym(tracks[i]->GetCov());
		}
		VertexFit* Vtx = new VertexFit(Ntr, pr, cv);
		//std::cout << "TestVertex: after new Vtx" << std::endl;
		TVectorD xvtx = Vtx->GetVtx();		// Ntr-2 fit
		//std::cout << "TestVertex: first fit completed. Ntr = "<<Ntr << std::endl;
		Double_t ran1 = gRandom->Rndm();
		Int_t it1 = TMath::Nint(ran1*(Ntr - 1));
		Int_t it2 = -1;
		while (it2 < 0 || it2 == it1){
			Double_t ran2 = gRandom->Rndm();
			it2 = TMath::Nint(ran2*(Ntr - 2));
		}
		//std::cout << "Removing tracks: " << it1 << " and " << it2 << std::endl;
		Vtx->RemoveTrk(it1);
		Vtx->RemoveTrk(it2);
		//std::cout << "Tracks removed" << std::endl;
		//
		xvtx = Vtx->GetVtx();					// Updated fit
		//std::cout << "TestVertex: updated fit completed" << std::endl;
		TMatrixDSym covX = Vtx->GetVtxCov();
		Double_t Chi2 = Vtx->GetVtxChi2();
		*/
		delete Vtx;	// Cleanup
		//std::cout << "TestVertex: after delete Vtx" << std::endl;
		//
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
		//std::cout << "TestVertex: before final cleanup" << std::endl;
		for (Int_t i = 0; i < Ntr; i++) delete tracks[i];
		/*
		for (Int_t i = 0; i < Ntr; i++) {
			pr[i]->Clear(); delete pr[i];
		}
		for (Int_t i = 0; i < Ntr; i++) {
			cv[i]->Clear(); delete cv[i];
		}
		*/
		delete[] tracks;
		//delete[] pr; 
		//delete[] cv; 
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

