//
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
#include "SolGeom.h"
#include "TrkUtil.h"
#include "SolGridCov.h"
#include "ObsTrk.h"
#include "VertexFit.h"
#include "VertexMore.h"
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
void TestVtxMom(Int_t Nvtx = 100, Int_t Ntr = 2)
{
	//
	// Init geometry
	SolGeom *G = new SolGeom("GeoIDEA_BASE.txt");
	Double_t Bz = G->B();
	SolGridCov *GC = new SolGridCov();
	GC->Read("CovIDEA_BASE.root");			// Read in covariance array
	//
	// Ranges
	//
	Double_t ThDegMin = 40.0;
	Double_t ThDegMax = 140.0;
	Double_t Lmin = 0.002;
	Double_t Lmax = 0.010;
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
	// new parameter pulls
	TH1D *hDpull = new TH1D("hDpull", "Pull new D", Nbin, -10., 10.);
	TH1D *hP0pull = new TH1D("hP0pull", "Pull new #varphi_0", Nbin, -10., 10.);
	TH1D *hCpull = new TH1D("hCpull", "Pull new C", Nbin, -10., 10.);
	TH1D *hZ0pull = new TH1D("hZ0pull", "Pull new #Z_0", Nbin, -10., 10.);
	TH1D *hCtpull = new TH1D("hCtpull", "Pull new #lambda", Nbin, -10., 10.);
	//
	// old parameter pulls
	TH1D *hDpull0 = new TH1D("hDpull0", "Pull old D", Nbin, -10., 10.);
	TH1D *hP0pull0 = new TH1D("hP0pull0", "Pull old #varphi_0", Nbin, -10., 10.);
	TH1D *hCpull0 = new TH1D("hCpull0", "Pull old C", Nbin, -10., 10.);
	TH1D *hZ0pull0 = new TH1D("hZ0pull0", "Pull old #Z_0", Nbin, -10., 10.);
	TH1D *hCtpull0 = new TH1D("hCtpull0", "Pull old #lambda", Nbin, -10., 10.);
	//
	// Track momenta pulls
	TH1D* hpXpull = new TH1D("hpXpull", "Pull track momentum x component", Nbin, -10., 10.);
	TH1D* hpYpull = new TH1D("hpYpull", "Pull track momentum y component", Nbin, -10., 10.);
	TH1D* hpZpull = new TH1D("hpZpull", "Pull track momentum z component", Nbin, -10., 10.);
	//
	// Total momentum pulls
	TH1D* hptotXpull = new TH1D("hptotXpull", "Pull total momentum x component", Nbin, -10., 10.);
	TH1D* hptotYpull = new TH1D("hptotYpull", "Pull total momentum y component", Nbin, -10., 10.);
	TH1D* hptotZpull = new TH1D("hptotZpull", "Pull total momentum z component", Nbin, -10., 10.);
	//
	// vertex parameter pulls
	TH1D *hDpullV = new TH1D("hDpullV", "Pull vtx D", Nbin, -10., 10.);
	TH1D *hP0pullV = new TH1D("hP0pullV", "Pull vtx #varphi_0", Nbin, -10., 10.);
	TH1D *hCpullV = new TH1D("hCpullV", "Pull vtx C", Nbin, -10., 10.);
	TH1D *hZ0pullV = new TH1D("hZ0pullV", "Pull vtx #Z_0", Nbin, -10., 10.);
	TH1D *hCtpullV = new TH1D("hCtpullV", "Pull vtx #lambda", Nbin, -10., 10.);
	//
	// Loop on # vertices
	//
	Int_t NbadX = 0;
	Int_t NbadY = 0;
	Int_t NbadZ = 0;
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
		if(n%500 == 0)cout << "Event # " << n << ", Rvertex = " << sqrt(x(0) * x(0) + x(1) * x(1)) << endl;
		//cout << "True vertex: x = " << x(0) << ", y = " << x(1) << ", z = " << x(2) << endl;
		//
		// Loop on tracks
		ObsTrk **tracks = new ObsTrk*[Ntr];	// Smear tracks according to covariance matrix
		Double_t *ppa = new Double_t[Ntr];
		Double_t *tha = new Double_t[Ntr];
		Double_t *pha = new Double_t[Ntr];
		TVector3* pg = new TVector3[Ntr];
		Int_t i = 0;
		TVector3 p_in;
		Double_t Qgen = 0.0;
		Double_t MinR = 10000.0;
		while (i < Ntr)
		{
			Double_t rnP = gRandom->Rndm();
			Double_t Ptot = Pmin + rnP * (Pmax - Pmin);
			ppa[i] = Ptot;
			Double_t ThP = gRandom->Gaus(Th, dTheta);
			if(ThP > TMath::Pi())ThP = TMath::TwoPi()-ThP;
			if(ThP < 0.0)ThP = -ThP;		
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
			pg[i] = P;
			p_in += P;
			Double_t Q = 1.0;
			if (gRandom->Rndm() > 0.5)Q = -1.0;
			if(GC->IsAccepted(P)){
				Qgen += Q;
				tracks[i] = new ObsTrk(x, P, Q, GC, G);	
				TVector3 Xfrst = tracks[i]->GetFirstHit();
				if(Xfrst.Pt() < MinR) MinR = Xfrst.Pt();
				//std::cout << "i=" << i << ", Gen p= (" << pg[i].X() << ", " << pg[i].Y() << ", " << pg[i].Z() << "), Q= "<<Q << std::endl;
				//std::cout<<"Track "<<i<<" is accepted"<<std::endl;
				i++;
			}
			else{
				std::cout<<"Event "<<n<<", Track "<<i<<" is NOT accepted"<<std::endl;
			}
		}
		//
		//********************************************
		// Invocation of vertexfit    ****************
		//********************************************
		//
		//std::cout << "TestVertex: before constructor" << std::endl;
		
		VertexFit* Vtx = new VertexFit(Ntr, tracks);
		//std::cout << "TestVertex: after constructor" << std::endl;
		Vtx->SetStartR(MinR);
		Double_t Chi2 = Vtx->GetVtxChi2();
		if(Chi2/(2.0*Ntr-3.0) > 25.)Vtx->SetStartR(MinR);
		TVectorD xvtx = Vtx->GetVtx();
		//std::cout << "TestVertex: after GetVtx()" << std::endl;
		TMatrixDSym covX = Vtx->GetVtxCov();
		Chi2 = Vtx->GetVtxChi2();
		//
		// VertexMore stuff
		VertexMore* VM = new VertexMore(Vtx);
		TVector3 pTotG(0.,0.,0.);
		TVector3 pTotR(0.,0.,0.);
		for (Int_t i = 0; i < Ntr; i++) {
			TVector3 pr = VM->GetMomentum(i);
			Double_t Qr = VM->GetCharge(i);
			//std::cout << "i=" << i << ", Rec p= (" << pr(0) << ", " << pr(1) << ", " << pr(2) << "), Q= " << Qr << std::endl;
			//
			// Track momentum pulls
			TVector3 dp = pr - pg[i];
			//cout << "Momentum difference:"; dp.Print();
			TMatrixDSym pCov = VM->GetMomentumC(i);
			Double_t pXpull = dp(0) / TMath::Sqrt(pCov(0, 0));
			hpXpull->Fill(pXpull);
			Double_t pYpull = dp(1) / TMath::Sqrt(pCov(1, 1));
			hpYpull->Fill(pYpull);
			Double_t pZpull = dp(2) / TMath::Sqrt(pCov(2, 2));
			hpZpull->Fill(pZpull);
			//
			pTotG += pg[i];
			pTotR += pr;
			//
		}
		pTotR = VM->GetTotalP();
		TVector3 dptot = pTotR - pTotG;
		TMatrixDSym Ctot = VM->GetTotalPcov();
		Double_t ptotXpull = dptot(0) / TMath::Sqrt(Ctot(0, 0));
		hptotXpull->Fill(ptotXpull);
		Double_t ptotYpull = dptot(1) / TMath::Sqrt(Ctot(1, 1));
		hptotYpull->Fill(ptotYpull);
		Double_t ptotZpull = dptot(2) / TMath::Sqrt(Ctot(2, 2));
		hptotZpull->Fill(ptotZpull);
		//
		// Vertex parameter pulls
		TVectorD Vpar = VM->GetVpar();
		TMatrixDSym VparC = VM->GetVcov();
		TVectorD Gpar = TrkUtil::XPtoPar(x,p_in, Qgen, Bz);
		hDpullV->Fill((Vpar(0)-Gpar(0))/TMath::Sqrt(VparC(0,0)));
		hP0pullV->Fill((Vpar(1)-Gpar(1))/TMath::Sqrt(VparC(1,1)));
		hCpullV->Fill((Vpar(2)-Gpar(2))/TMath::Sqrt(VparC(2,2)));
		hZ0pullV->Fill((Vpar(3)-Gpar(3))/TMath::Sqrt(VparC(3,3)));
		hCtpullV->Fill((Vpar(4)-Gpar(4))/TMath::Sqrt(VparC(4,4)));
		//
		TMatrixDSym Mbig = VM->GetBigCov();
		//std::cout << "Big covariance:"; Mbig.Print();
		Double_t Mass = 5.28;
		Double_t *masses = new Double_t[Ntr];
		for(Int_t im=0; im<Ntr; im++)masses[i] = 0.135;
		VM->MassConstraint(Mass, masses);
		delete VM;
		//
		delete[] pg;
		//
		//
		// Parameter pulls	
		for(Int_t i=0; i<Ntr; i++){
			TVectorD truPar = tracks[i]->GetGenPar();
		// Old
			TVectorD obsPar = tracks[i]->GetObsPar();
			TMatrixDSym obsCov = tracks[i]->GetCov();
			hDpull0 ->Fill((obsPar(0)-truPar(0))/TMath::Sqrt(obsCov(0,0)));
			hP0pull0->Fill((obsPar(1)-truPar(1))/TMath::Sqrt(obsCov(1,1)));
			hCpull0 ->Fill((obsPar(2)-truPar(2))/TMath::Sqrt(obsCov(2,2)));
			hZ0pull0->Fill((obsPar(3)-truPar(3))/TMath::Sqrt(obsCov(3,3)));
			hCtpull0->Fill((obsPar(4)-truPar(4))/TMath::Sqrt(obsCov(4,4)));
		// New	
			TVectorD newPar = Vtx->GetNewPar(i);
			TMatrixDSym newCov = Vtx->GetNewCov(i);
			hDpull ->Fill((newPar(0)-truPar(0))/TMath::Sqrt(newCov(0,0)));
			hP0pull->Fill((newPar(1)-truPar(1))/TMath::Sqrt(newCov(1,1)));
			hCpull ->Fill((newPar(2)-truPar(2))/TMath::Sqrt(newCov(2,2)));
			hZ0pull->Fill((newPar(3)-truPar(3))/TMath::Sqrt(newCov(3,3)));
			hCtpull->Fill((newPar(4)-truPar(4))/TMath::Sqrt(newCov(4,4)));
		}
		delete Vtx;
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
		delete[] tracks;
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
	TCanvas *cnv1 = new TCanvas("cnv1", "New parameter pulls", 100, 100, 500, 500);
	cnv1->Divide(3, 2);
	cnv1->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hDpull->Fit("gaus");
	hDpull->Draw();
	cnv1->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hP0pull->Fit("gaus");
	hP0pull->Draw();
	cnv1->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hCpull->Fit("gaus");
	hCpull->Draw();
	cnv1->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hZ0pull->Fit("gaus");
	hZ0pull->Draw();
	cnv1->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hCtpull->Fit("gaus");
	hCtpull->Draw();
	//
	TCanvas *cnv2 = new TCanvas("cnv2", "Old parameter pulls", 200, 100, 500, 500);
	cnv2->Divide(3, 2);
	cnv2->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hDpull0->Fit("gaus");
	hDpull0->Draw();
	cnv2->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hP0pull0->Fit("gaus");
	hP0pull0->Draw();
	cnv2->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hCpull0->Fit("gaus");
	hCpull0->Draw();
	cnv2->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hZ0pull0->Fit("gaus");
	hZ0pull0->Draw();
	cnv2->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hCtpull0->Fit("gaus");
	hCtpull0->Draw();
	//
	//
	TCanvas* cnv3 = new TCanvas("cnv3", "Track momenta pulls", 250, 150,900, 600);
	cnv3->Divide(3, 2);
	cnv3->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hpXpull->Fit("gaus");
	hpXpull->Draw();
	cnv3->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hpYpull->Fit("gaus");
	hpYpull->Draw();
	cnv3->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hpZpull->Fit("gaus");
	hpZpull->Draw();
	//
	cnv3->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hptotXpull->Fit("gaus");
	hptotXpull->Draw();
	cnv3->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hptotYpull->Fit("gaus");
	hptotYpull->Draw();
	cnv3->cd(6); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hptotZpull->Fit("gaus");
	hptotZpull->Draw();
	//	
	//
	TCanvas *cnv4 = new TCanvas("cnv2", "Vertex parameter pulls", 300, 300, 500, 500);
	cnv4->Divide(3, 2);
	cnv4->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hDpullV->Fit("gaus");
	hDpullV->Draw();
	cnv4->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hP0pullV->Fit("gaus");
	hP0pullV->Draw();
	cnv4->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hCpullV->Fit("gaus");
	hCpullV->Draw();
	cnv4->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hZ0pullV->Fit("gaus");
	hZ0pullV->Draw();
	cnv4->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(1111);
	hCtpullV->Fit("gaus");
	hCtpullV->Draw();
}

