#include "SolGeom.h"
#include "TrkUtil.h"
#include "SolGridCov.h"
#include "ObsTrk.h"
#include "VertexFit.h"
#include "VertexMore.h"
#include "BsDecay.h"
#include <TRandom.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>

void TestBsDecay(Double_t p = 10., Int_t Nev = 2)
{
	//
	// Init geometry
	SolGeom *G = new SolGeom("GeoIDEA_BASE.txt");
	Double_t Bz = G->B();
	SolGridCov *GC = new SolGridCov();
	GC->Read("CovIDEA_BASE.root");			// Read in covariance array
	//
	// Define histograms
	// Ds vertex
	TH1D* h_DsXv = new TH1D("h_DsXv","Ds X-vertex pull",100,-10.,10.);
	TH1D* h_DsYv = new TH1D("h_DsYv","Ds Y-vertex pull",100,-10.,10.);
	TH1D* h_DsZv = new TH1D("h_DsZv","Ds Z-vertex pull",100,-10.,10.);
	// Ds vertex momentum
	TH1D* h_DsPx = new TH1D("h_DsPx","Ds X-momentum pull",100,-10.,10.);
	TH1D* h_DsPy = new TH1D("h_DsPy","Ds Y-momentum pull",100,-10.,10.);
	TH1D* h_DsPz = new TH1D("h_DsPz","Ds Z-momentum pull",100,-10.,10.);
	// Bs vertex
	TH1D* h_BsXv = new TH1D("h_BsXv","Bs X-vertex pull",100,-10.,10.);
	TH1D* h_BsYv = new TH1D("h_BsYv","Bs Y-vertex pull",100,-10.,10.);
	TH1D* h_BsZv = new TH1D("h_BsZv","Bs Z-vertex pull",100,-10.,10.);
	// Bs vertex momentum
	TH1D* h_BsPx = new TH1D("h_BsPx","Bs X-momentum pull",100,-10.,10.);
	TH1D* h_BsPy = new TH1D("h_BsPy","Bs Y-momentum pull",100,-10.,10.);
	TH1D* h_BsPz = new TH1D("h_BsPz","Bs Z-momentum pull",100,-10.,10.);
	
	//
	// Event loop
	Double_t rx, ry, rz;
	TVector3 pBs;
	for(Int_t n=0; n<Nev; n++){
		Double_t Theta = 0.;
		while(Theta <30. || Theta > 60.){
			gRandom->Sphere(rx, ry, rz, p);
			pBs.SetXYZ(rx, ry, rz);
			Theta = pBs.Theta()*180./TMath::Pi();
		}
		BsDecay Bs(pBs);
		// 
		// Check printout
		if(n<3){
			//
			// Bs decay
			Int_t ind = Bs.GetBs();
			std::cout<<"Bs: ("<<Bs.GetPtype(ind)<<"), mass = "<<
			Bs.GetPmass(ind)<<", Path = "<<Bs.GetPpath(ind)<<std::endl;
			ind = Bs.GetDs();
			std::cout<<"Ds: ("<<Bs.GetPtype(ind)<<"), mass = "<<
			Bs.GetPmass(ind)<<", Path = "<<Bs.GetPpath(ind)<<std::endl;
			ind = Bs.GetPiBs();
			std::cout<<"PiBs: ("<<Bs.GetPtype(ind)<<"), mass = "<<
			Bs.GetPmass(ind)<<", Path = "<<Bs.GetPpath(ind)<<std::endl;
			//
			std::cout<<"Bs momentum";
			ind = Bs.GetBs();
			TVector3 pp = Bs.GetMomentum(ind);
			Double_t mass = Bs.GetPmass(ind);
			Double_t EE = TMath::Sqrt(pp.Mag2()+mass*mass);
			TLorentzVector ptot(pp,EE);
			pp.Print();
			std::cout<<"Ds momentum";
			ind = Bs.GetDs();
			pp = Bs.GetMomentum(ind);
			mass = Bs.GetPmass(ind);
			EE = TMath::Sqrt(pp.Mag2()+mass*mass);
			TLorentzVector p1(pp,EE);
			pp.Print();
			std::cout<<"PiBs momentum";
			ind = Bs.GetPiBs();
			pp = Bs.GetMomentum(ind);
			mass = Bs.GetPmass(ind);
			EE = TMath::Sqrt(pp.Mag2()+mass*mass);
			TLorentzVector p2(pp,EE);
			pp.Print();
			// invariant mass test
			TLorentzVector BsLV = p1+p2;
			std::cout<<"Ds-Pi+ invariant mass = "<<BsLV.M()<<std::endl;
			//
			// Ds decay
			ind = Bs.GetDs();
			std::cout<<"Ds: ("<<Bs.GetPtype(ind)<<"), mass = "<<
			Bs.GetPmass(ind)<<", Path = "<<Bs.GetPpath(ind)<<std::endl;
			ind = Bs.GetPhi();
			std::cout<<"Phi: ("<<Bs.GetPtype(ind)<<"), mass = "<<
			Bs.GetPmass(ind)<<", Path = "<<Bs.GetPpath(ind)<<std::endl;
			ind = Bs.GetPiDs();
			std::cout<<"PiDs: ("<<Bs.GetPtype(ind)<<"), mass = "<<
			Bs.GetPmass(ind)<<", Path = "<<Bs.GetPpath(ind)<<std::endl;
			//
			std::cout<<"Ds momentum";
			ind = Bs.GetDs();
			pp = Bs.GetMomentum(ind);
			pp.Print();
			std::cout<<"Phi momentum";
			ind = Bs.GetPhi();
			pp = Bs.GetMomentum(ind);
			pp.Print();
			std::cout<<"PiDs momentum";
			ind = Bs.GetPiDs();
			pp = Bs.GetMomentum(ind);
			pp.Print();
			//
			// Phi decay
			ind = Bs.GetPhi();
			std::cout<<"Phi: ("<<Bs.GetPtype(ind)<<"), mass = "<<
			Bs.GetPmass(ind)<<", Path = "<<Bs.GetPpath(ind)<<std::endl;
			ind = Bs.GetKpPhi();
			std::cout<<"K+: ("<<Bs.GetPtype(ind)<<"), mass = "<<
			Bs.GetPmass(ind)<<", Path = "<<Bs.GetPpath(ind)<<std::endl;
			ind = Bs.GetKmPhi();
			std::cout<<"K-: ("<<Bs.GetPtype(ind)<<"), mass = "<<
			Bs.GetPmass(ind)<<", Path = "<<Bs.GetPpath(ind)<<std::endl;
			//
			std::cout<<"Phi momentum";
			ind = Bs.GetPhi();
			pp = Bs.GetMomentum(ind);
			pp.Print();
			std::cout<<"K+ momentum";
			ind = Bs.GetKpPhi();
			pp = Bs.GetMomentum(ind);
			pp.Print();
			std::cout<<"K- momentum";
			ind = Bs.GetKmPhi();
			pp = Bs.GetMomentum(ind);
			pp.Print();
		}	// end check printout
		//	
		// Store generated vertices
		Int_t iBs = Bs.GetBs();
		TVector3 xBs = (Bs.GetPpath(iBs)/pBs.Mag())*pBs;	// Bs vertex
		Int_t iDs = Bs.GetDs();
		TVector3 pDs = Bs.GetMomentum(iDs);
		TVector3 xDs = (Bs.GetPpath(iDs)/pDs.Mag())*pDs + xBs;	// Ds vertex
		//
		// Tracks for Ds vertex fit
		const Int_t NtrDs = 3;
		Int_t iTrDs[3];
		iTrDs[0] = Bs.GetPiDs();	// pi- from Ds decay
		iTrDs[1] = Bs.GetKpPhi();	// K+ from Phi
		iTrDs[2] = Bs.GetKmPhi();	// K- from Phi
		//
		ObsTrk **tracks = new ObsTrk*[NtrDs];
		Int_t NtPass = 0;
		for(Int_t i=0; i<NtrDs; i++){
			TVector3 pTr = Bs.GetMomentum(iTrDs[i]);
			if(GC->IsAccepted(xDs, pTr, G)){
				Double_t Q = TMath::Sign(1.0,(Double_t)Bs.GetPtype(iTrDs[i]));
				tracks[NtPass] = new ObsTrk(xDs, pTr, Q, GC, G);
				NtPass++;
			}else std::cout<<"Ds track out of acceptance"<<std::endl;
		} // End Ds track loop
		// Fit Ds
		if(NtPass == NtrDs){
			VertexFit * vDs = new VertexFit(NtrDs, tracks);
			//Double_t rDs = xDs.Pt();
			//vDs->SetStartR(rDs);
			TVectorD xvDs = vDs->GetVtx();
			TMatrixDSym covXvDs = vDs->GetVtxCov();
			Double_t Chi2Ds = vDs->GetVtxChi2();
			//
			// Tracks for Bs vertex fit
			const Int_t NtrBs = 2;
			TVectorD **trPar = new TVectorD*[NtrBs];
			TMatrixDSym **trCov = new TMatrixDSym*[NtrBs];
			// Set Ds vertex track
			VertexMore* VM = new VertexMore(vDs);
			trPar[0] = new TVectorD(VM->GetVpar());
			trCov[0] = new TMatrixDSym(VM->GetVcov());
			Int_t NtPassB = 0;
			Int_t iPiBs = Bs.GetPiBs();
			TVector3 piTr = Bs.GetMomentum(iPiBs);
			if(GC->IsAccepted(xBs, piTr, G)){
				Double_t Q = TMath::Sign(1.0,(Double_t)Bs.GetPtype(iPiBs));
				ObsTrk *TrPi = new ObsTrk(xBs, piTr, Q, GC, G);
				trPar[1] = new TVectorD(TrPi->GetObsPar());
				trCov[1] = new TMatrixDSym(TrPi->GetCov());
				NtPassB++;
				}else std::cout<<"Bs track out of acceptance"<<std::endl;
			// Fit Bs
			if(NtPassB == 1){
				//
				// Fill Ds histograms
				// Vertex
				h_DsXv->Fill((xvDs(0)-xDs(0))/TMath::Sqrt(covXvDs(0,0)));
				h_DsYv->Fill((xvDs(1)-xDs(1))/TMath::Sqrt(covXvDs(1,1)));
				h_DsZv->Fill((xvDs(2)-xDs(2))/TMath::Sqrt(covXvDs(2,2)));
				// Momentum
				TVector3 pDsRec = VM->GetTotalP();
				TMatrixDSym covDsP = VM->GetTotalPcov();
				h_DsPx->Fill((pDsRec(0)-pDs(0))/TMath::Sqrt(covDsP(0,0)));
				h_DsPy->Fill((pDsRec(1)-pDs(1))/TMath::Sqrt(covDsP(1,1)));
				h_DsPz->Fill((pDsRec(2)-pDs(2))/TMath::Sqrt(covDsP(2,2)));
				//
				// Fit Bs vertex
				VertexFit* vBs = new VertexFit(NtrBs, trPar, trCov);
				//Double_t rBs = xBs.Pt();
				//vBs->SetStartR(rBs);
				TVectorD xvBs = vBs->GetVtx();
				TMatrixDSym covXvBs = vBs->GetVtxCov();
				Double_t Chi2Bs = vBs->GetVtxChi2();
				VertexMore* VMs = new VertexMore(vBs);
				//
				// Fill Bs histograms
				// Vertex
				h_BsXv->Fill((xvBs(0)-xBs(0))/TMath::Sqrt(covXvBs(0,0)));
				h_BsYv->Fill((xvBs(1)-xBs(1))/TMath::Sqrt(covXvBs(1,1)));
				h_BsZv->Fill((xvBs(2)-xBs(2))/TMath::Sqrt(covXvBs(2,2)));
				// Momentum
				TVector3 pBsRec = VMs->GetTotalP();
				TMatrixDSym covBsP = VMs->GetTotalPcov();
				h_BsPx->Fill((pBsRec(0)-pBs(0))/TMath::Sqrt(covBsP(0,0)));
				h_BsPy->Fill((pBsRec(1)-pBs(1))/TMath::Sqrt(covBsP(1,1)));
				h_BsPz->Fill((pBsRec(2)-pBs(2))/TMath::Sqrt(covBsP(2,2)));
				//
				// Clean
				delete vBs;
				delete VMs;
			}
			//
			// Clean
			for (Int_t i = 0; i < NtPassB+1; i++){
				trPar[i]->Clear();
				trCov[i]->Clear();
			}
			delete[] trPar;
			delete[] trCov;
			delete vDs;
			delete VM;
		}
		//
		// Clean
		// std::cout<<"Event # "<<n<<", NtPass = "<<NtPass<<std::endl;
		for (Int_t i = 0; i < NtPass; i++) delete tracks[i];
		delete[] tracks;
		
	}	// End event loop
	//
	// Display results
	//
	// Ds plots
	TCanvas * cnv = new TCanvas("cnv","Ds vertex pulls",10,10, 900,600);
	cnv->Divide(3,2);
	// X
	cnv->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_DsXv->Fit("gaus"); h_DsXv->Draw();
	// Y
	cnv->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_DsYv->Fit("gaus"); h_DsYv->Draw();
	// Z
	cnv->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_DsZv->Fit("gaus"); h_DsZv->Draw();
	// Px
	cnv->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_DsPx->Fit("gaus"); h_DsPx->Draw();
	// Py
	cnv->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_DsPy->Fit("gaus"); h_DsPy->Draw();
	// Pz
	cnv->cd(6); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_DsPz->Fit("gaus"); h_DsPz->Draw();
	//
	// Bs plots
	TCanvas * cnv1 = new TCanvas("cnv1","Ds vertex pulls",100,100, 900,600);
	cnv1->Divide(3,2);
	// X
	cnv1->cd(1); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_BsXv->Fit("gaus"); h_BsXv->Draw();
	// Y
	cnv1->cd(2); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_BsYv->Fit("gaus"); h_BsYv->Draw();
	// Z
	cnv1->cd(3); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_BsZv->Fit("gaus"); h_BsZv->Draw();
	// Px
	cnv1->cd(4); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_BsPx->Fit("gaus"); h_BsPx->Draw();
	// Py
	cnv1->cd(5); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_BsPy->Fit("gaus"); h_BsPy->Draw();
	// Pz
	cnv1->cd(6); gPad->SetLogy(1);
	gStyle->SetOptStat(111111); gStyle->SetOptFit(1111);
	h_BsPz->Fit("gaus"); h_BsPz->Draw();
}
