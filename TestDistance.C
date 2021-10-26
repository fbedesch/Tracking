#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TVectorD.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TGraph.h>
#include <iostream>
#include "SolGeom.h"
#include "TrkUtil.h"
#include "SolTrack.h"
#include "SolGridCov.h"
#include "ObsTrk.h"
#include "VertexFit.h"

Double_t mpip = 0.13957061;
Double_t mpi0 = 0.1349770;
Double_t mrhop = 0.7754;
Double_t mtau = 1.77682;
Double_t mk0s = 0.497611;


const Int_t Npx = 100;

void PlotTrack(TVectorD pr, Double_t phiMx, Double_t *x, Double_t *y)
{
	//

	for (Int_t i = 0; i < Npx; i++)
	{
		Double_t ph = i * phiMx / (Double_t) Npx;
		x[i] = -pr(0) * sin(pr(1)) + (sin(pr(1) + ph) - sin(pr(1))) / (2 * pr(2));
		y[i] =  pr(0) * cos(pr(1)) - (cos(pr(1) + ph) - cos(pr(1))) / (2 * pr(2));
	}
	//
}

TLorentzVector Boosted(TLorentzVector Booster, TLorentzVector InMoment)
{
	//
	// Get boost features from 4-vector
	//
	Double_t gamma = Booster.Gamma();
	Double_t beta = Booster.Beta();
	TVector3 nBoost = (1.0 / Booster.Rho()) * Booster.Vect();	// Boost versor 
	Double_t pl = nBoost.Dot(InMoment.Vect());				// Pre-boost longitudinal component
	TVector3 pt = InMoment.Vect() - pl * nBoost;			// Transverse component
	Double_t pll = gamma * (pl + beta * InMoment.E());		// Boosted logitudinal component
	TVector3 pout = pt + pll * nBoost;						// Boosted momentum
	Double_t Eout = gamma * (InMoment.E() + beta * pl);		// Boosted energy
	TLorentzVector Pboosted(pout, Eout);					// Boosted LorentzVector
	//
	return Pboosted;
}
void TwoBodyPS(TLorentzVector pin, Double_t m1, Double_t m2, TLorentzVector**& pout)
{
	Double_t Min2 = pin.M2();			// Incoming particle mass squared
	Double_t Min = TMath::Sqrt(Min2);	// Incoming particle mass
	if (Min < m1 + m2) {
		std::cout << "TwoBodyPS: incoming mass small than sum of outgoing. Abort." << std::endl;
		exit(EXIT_FAILURE);
	}
	//
	// Decay particle momenta in CM frame
	Double_t E1 = (Min2 + m1 * m1 - m2 * m2) / (2 * Min);	// Particle 1 energy
	Double_t E2 = (Min2 + m2 * m2 - m1 * m1) / (2 * Min);	// Particle 2 energy
	Double_t p = TMath::Sqrt((Min2 - (m1 - m2) * (m1 - m2)) * (Min2 - (m1 + m2) * (m1 + m2))) / (2 * Min); // momentum of outgoing particles
	// Assume flat in phi and cos(theta)
	Double_t CosTh = 2 * gRandom->Rndm() - 1.0;
	Double_t SinTh = TMath::Sqrt(1.0 - CosTh * CosTh);
	Double_t Phi = TMath::TwoPi() * gRandom->Rndm();
	TLorentzVector p1cm(p * SinTh * TMath::Cos(Phi), p * SinTh * TMath::Sin(Phi), p * CosTh, E1);	// CM momentum of part. 1
	TLorentzVector p2cm(-p * SinTh * TMath::Cos(Phi), -p * SinTh * TMath::Sin(Phi), -p * CosTh, E2);	// CM momentum of part. 1
	//
	// Boost to lab frame
	//
	pout[0] = new TLorentzVector(Boosted(pin, p1cm));
	pout[1] = new TLorentzVector(Boosted(pin, p2cm));
	//
}

void TestDistance(Double_t dist = 20.)
{
	//
	// Init geometry
	SolGeom* G = new SolGeom("GeoIDEA_BASE.txt");	// Read in geometry
	Double_t Bfield = G->B();				// Magnetic field
	SolGridCov* GC = new SolGridCov();
	GC->Read("CovIDEA-BASE.root");			// Read in covariance array
	//
	// Decaying particle momentum
	//
	TVector3 p0(1., 5., 0.5);
	Double_t Ein = sqrt(p0.Mag2() + mk0s * mk0s);
	TLorentzVector pin(p0, Ein);
	// Decay position
	TVector3 x = (dist / p0.Mag()) * p0;
	// Inner box boundaries
	Double_t Rin = G->GetRmin();
	Double_t ZinPos = G->GetZminPos();
	Double_t ZinNeg = G->GetZminNeg();
	Bool_t inside = TrkUtil::IsInside(x, Rin, ZinNeg, ZinPos); // Check if in inner box
	//
	// Event loop
	//
	Bool_t ask = kTRUE;
	while (ask)
	{
		TLorentzVector** pout = new TLorentzVector * [2];
		TwoBodyPS(pin, mpip, mpip, pout);
		cout << "Input vector:"; pin.Print();
		cout << "pout[0]:"; pout[0]->Print();
		cout << "pout[1]:"; pout[1]->Print();
		cout << "M = " << pin.M() << ", m1 = " << pout[0]->M() << ", m2 = " << pout[1]->M() << endl;
		Int_t Nacc = 0;	// # of tracks accepted	
		TVectorD** pr = new TVectorD * [2];
		TMatrixDSym** cv = new TMatrixDSym * [2];
		Double_t Q = 1.0;
		if (gRandom->Rndm() > 0.5)Q = -1.0;
		for (Int_t i = 0; i < 2; i++)
		{
			Q *= -1.;
			TVector3 po = pout[i]->Vect();
			TVectorD prGen = TrkUtil::XPtoPar(x, po, Q, Bfield);	// Generated parameters
			TMatrixDSym Cov(5);
			Bool_t accpt = kFALSE;
			if (inside)
			{
				if (GC->IsAccepted(po))			// Check track simplified acceptance
				{
					Nacc++; accpt = kTRUE;
					// Observed track parameters
					Double_t pt = po.Pt();
					Double_t angd = po.Theta() * 180. / TMath::Pi();
					Cov = GC->GetCov(pt, angd);				// Track covariance
				}
			}
			else
			{
				SolTrack* trk = new SolTrack(x, po, G);
				if (trk->nmHit() >= GC->GetMinHits())		// Check precise track acceptance
				{
					cout << "SolTrack (" << i << "). #hits = " << trk->nmHit() << endl;
					//cout << "Measured hits = " << trk->nmHit() << endl;
					Nacc++; accpt = kTRUE;
					// Track covariance
					Bool_t Res = kTRUE; Bool_t MS = kTRUE;
					trk->CovCalc(Res, MS);					// Calculate covariance matrix
					//cout << "1 ==" << endl;
					Cov = trk->Cov();						// Track covariance
					//cout << "2 ==" << endl;
				}
				//cout << "3 ==" << endl;
				delete trk;
			}
			//cout << "3.5 ==" << endl;
			if (accpt)
			{
				TVectorD prObs = TrkUtil::CovSmear(prGen, Cov);		// Observed parameters
				//cout << "4 ==" << endl;
				// Fill vertexing arrays
				pr[i] = new TVectorD(prObs);						// Track list for vertexing
				cv[i] = new TMatrixDSym(Cov);						// Covariance list for vertexing
				if (TMath::IsNaN(prObs(0)))
				{
					cout << "prObs is NaN. i= " << i << endl;
					cout << "po: "; po.Print();
					cout << "prGen: "; prGen.Print();
					cout << "Cov:"; Cov.Print();
				}
			}
			// clean decay particles
			pout[i]->Clear(); delete pout[i];
		}
		delete[] pout;
		//cout << "Nacc = " << Nacc << endl;
		//
		if (Nacc == 2)
		{
			//
			//********************************************
			// Invocation of vertexfit    ****************
			//********************************************
			VertexFit* Vtx = new VertexFit(Nacc, pr, cv);
			//cout << "1 ==" << endl;
			TVectorD xvtx = Vtx->GetVtx();
			//cout << "2 ==" << endl;
			TMatrixDSym covX = Vtx->GetVtxCov();
			//cout << "3 ==" << endl;
			Double_t Chi2 = Vtx->GetVtxChi2();
			//cout << "4 ==" << endl;
			// clean vertex
			delete Vtx;
			//
			cout << "x gen = (" << x(0) << ", " << x(1) << ", " << x(2) << endl;
			cout << "x found = (" << xvtx(0) << ", " << xvtx(1) << ", " << xvtx(2) << endl;
			if (TMath::IsNaN(xvtx(0)))
			{
				cout << "pr[0]:"; pr[0]->Print(); cv[0]->Print();
				cout << "pr[1]:"; pr[1]->Print(); cv[1]->Print();
			}
		}
		//
		// Graphs
		//
		TCanvas* CC = new TCanvas("CC", "K0s decay", 10, 10, 500, 500);
		Double_t xmin = 1000;
		Double_t xmax = -1000.;
		Double_t ymin = 1000.;
		Double_t ymax = -1000.;
		Double_t x1[Npx];
		Double_t y1[Npx];
		Double_t x2[Npx];
		Double_t y2[Npx];
		for (Int_t i = 0; i < 2; i++)
		{
			TVectorD par = *pr[i];
			Double_t Rm = 2.0;
			Double_t phiMx = 2. * asin(par(2) * sqrt((Rm * Rm - par(0) * par(0)) / (1. + 2 * par(0) * par(2))));
			Double_t xx[Npx];
			Double_t yy[Npx];
			PlotTrack(par, phiMx, xx, yy);
			if (i == 0)for (Int_t j = 0; j < Npx; j++)x1[j] = xx[j];
			if (i == 0)for (Int_t j = 0; j < Npx; j++)y1[j] = yy[j];
			if (i == 1)for (Int_t j = 0; j < Npx; j++)x2[j] = xx[j];
			if (i == 1)for (Int_t j = 0; j < Npx; j++)y2[j] = yy[j];
			Double_t x0n = -par(0) * sin(par(1));
			if (x0n < xmin) xmin = x0n;
			Double_t y0n = par(0) * cos(par(1));
			if (y0n < ymin) ymin = y0n;
			Double_t x0x = x0n + (sin(phiMx + par(1)) - sin(par(1))) / (2 * par(2));
			if (x0x > xmax)xmax = x0x;
			Double_t y0x = y0n - (cos(phiMx + par(1)) - cos(par(1))) / (2 * par(2));
			if (y0x > ymax)ymax = y0x;
		}
		TGraph* g[2];
		g[0] = new TGraph(Npx, x1, y1);
		g[1] = new TGraph(Npx, x2, y2);
		g[0]->SetLineColor(kBlue);
		g[1]->SetLineColor(kRed);
		Double_t margin = 0.1;
		for (Int_t i = 0; i < 2; i++) {
			g[i]->GetXaxis()->SetLimits(xmin - margin, xmax + margin);
			g[i]->GetYaxis()->SetLimits(ymin - margin, ymax + margin);
		}
		//
		cout << "New event?y/n" << endl;
		CC->Draw();
		g[0]->Draw("AC");
		g[1]->Draw("SAME");
		CC->Update();
		char resp[2];
		cin >> resp;
		char no[2] = "n";
		if (strcmp(resp, no) == 0)ask = kFALSE;
		//
		//
		// Cleanup
		//
		//clean track parameters and covariance
		for (Int_t i = 0; i < Nacc; i++) {
			pr[i]->Clear(); delete pr[i];
		}
		for (Int_t i = 0; i < Nacc; i++) {
			cv[i]->Clear(); delete cv[i];
		}
		delete[] pr;
		delete[] cv;
		CC->Clear();
	}
}

