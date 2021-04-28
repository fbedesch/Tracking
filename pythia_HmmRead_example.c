
/// Note that before executing this script, for some Pythia8 builds:
///
///  - the env variable PYTHIA8 must point to the pythia8100 (or newer) directory
///  - the env variable PYTHIA8DATA must be defined and it must point to $PYTHIA8/xmldoc
///
/// Description:
/// Read in HZ-> mu mu events and calculate recoil mass
///
/// \author: Franco Bedeschi

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include "SolGeom.h"
#include "SolGridCov.h"
#include "ObsTrk.h"

void HmmRead(Int_t nev  = 100, TString fname = "Hmm.root", Int_t ndeb = 1)
{
	//
	// Load libraries
//	gSystem->Load("SolGeomIDEA");
//	gSystem->Load("SolTrack");
//	gSystem->Load("AcceptanceClx");
//	gSystem->Load("SolGridCov");
//	gSystem->Load("TrkUtil");
//	gSystem->Load("ObsTrk");

//
//	Relevant particle codes
//
	Int_t Muon = 13;		// Muon PDG code
	Int_t H0   = 25;		// Z0   PDG code
// Initialize B field
	Double_t Bfield = 2.0;			// Set B field in Tesla
// Initialize tracking resolution
	SolGridCov *GCid = new SolGridCov();		// IDEA geometry
	GCid->Read("CovIDEA-BASE.root");			// Read in covariance array
	SolGridCov *GCcl = new SolGridCov();		// CLD geometry
	GCcl->Read("CovCLD.root");					// Read in covariance array

//
// Histograms
   TH1D* etaMu  = new TH1D("etaMu",  "Pseudorapidity",                             120,  -3.,   3.);
   TH1D* ptMu   = new TH1D("ptMu",   "pt",                                         100,   0., 100.);
   TH1D* gMinv  = new TH1D("Minv",   "Generated di-muon invariant mass",           100,  115., 135.);
   TH1D* oMinv_id  = new TH1D("Minv_id", "IDEA: Observed  di-muon invariant mass", 200, 115., 135.);
   TH1D* oMinv_cl  = new TH1D("Minv_cl", "CLD:  Observed  di-muon invariant mass", 200, 115., 135.);

//
// Setup File and Tree for input
   TFile *fin = new TFile(fname,"READ");
   TTree *T = (TTree*)fin->Get("tID");
// Array of particles
   TClonesArray* particles = new TClonesArray("TParticle", 1000);
// Link to TTree
   T->SetBranchAddress("Particles", &particles);

//
   Int_t MaxEv = T->GetEntries();			// Number of event in file
   Int_t Nevents = TMath::Min(nev, MaxEv);	// reduce event if you want
// Event loop
   for (Int_t iev = 0; iev < Nevents; iev++)	// Main event loop
   {
	  Int_t nRead = T->GetEntry(iev);		// Get event
	  if (nRead <= 0) continue;
      Int_t np = particles->GetEntriesFast();	// Get # particles
// Particle loop
	  TLorentzVector p1(0., 0., 0., 0.);		// Initialize 4-momenta
	  TLorentzVector p2(0., 0., 0., 0.);
	  Double_t Q1 = 0; Double_t M1 = 0;
	  Double_t Q2 = 0; Double_t M2 = 0;
      Int_t Nmu = 0;
      for (Int_t ip = 0; ip < np; ip++)			// Main loop on event particles
	  {
         TParticle* part = (TParticle*) particles->At(ip);	// Get Prticle
         Int_t ist = part->GetStatusCode();					// Get PDG code
         // Positive codes are final particles.
         if (ist <= 0) continue;							// Skip if not final state
         Int_t pdg = part->GetPdgCode();
		 Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/3.0;
		 Float_t mass   = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
         if (charge == 0.) continue;						// Skip if neutral
		 // Select muons from Higgs
		 if (TMath::Abs(pdg) == Muon)		// Final state muon
		 {
			 Int_t Mother = part->GetFirstMother();
			 while (Mother > 0)				// Climb decay tree looking for Higgs
			 {
				 TParticle* mpart = (TParticle*)particles->At(Mother);
				 Int_t mpdg = mpart->GetPdgCode();
				 if (mpdg == H0)				// Found Higgs in decay tree
				 {
					 // Fill control histograms
					 Float_t eta = part->Eta();
					 Float_t pt  = part->Pt();
					 etaMu->Fill(eta);
					 ptMu ->Fill(pt);
					 // fill muon vectors
					 Nmu++;
					 if (Nmu == 1)				// Store first muon info
					 {
						 Q1 = charge; M1 = mass;
						 p1.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
					 }
					 if (Nmu == 2)				// Store second muon info
					 {
						 Q2 = charge; M2 = mass;
						 p2.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->Energy());
					 }
					 break;
				 }
				 Mother = mpart->GetFirstMother();
			 }
		 }
      }			// End main loop on particles in the event
	  //
	  if (Nmu == 2)	// Found two muons from H0. Now fill histograms.
	  {
		  TLorentzVector ptot = p1 + p2;
		  Double_t Mi = ptot.M();						// Generated mu mu invariant mass
		  gMinv->Fill(Mi);								// Fill generated invariant mass histogram
		  //
		  // Now apply track resolution effects
		  //
		  TVector3 tP1 = p1.Vect();		// Get generated momentum first  muon
		  TVector3 tP2 = p2.Vect();		// Get generated momentum second muon
		  // Acceptance above 10 degrees and pt>1 GeV
		  Double_t AngMin = 10.0;
		  Double_t th1 = TMath::ACos(TMath::Abs(tP1.CosTheta()))*180. / TMath::Pi();
		  Double_t th2 = TMath::ACos(TMath::Abs(tP2.CosTheta()))*180. / TMath::Pi();
		  //
		  if (th1 < AngMin || th2 < AngMin)continue;		// Angular acceptance cuts
		  if (tP1.Pt() < 1. || tP2.Pt() < 1.)continue;		// Pt cuts

		  //
		  // Fill IDEA Histograms
		  TVector3 tX(0.0, 0.0, 0.0);						// Set origin to (0,0,0)
		  ObsTrk *Tr1id = new ObsTrk(tX, tP1, Q1, Bfield, GCid);	// Apply track resolution
		  ObsTrk *Tr2id = new ObsTrk(tX, tP2, Q2, Bfield, GCid);
		  TVector3 obsP1id = Tr1id->GetObsP();					// Get observed momenta
		  TVector3 obsP2id = Tr2id->GetObsP();
		  Double_t E1id = TMath::Sqrt(M1*M1 + obsP1id.Mag2());	// Observed energies
		  Double_t E2id = TMath::Sqrt(M2*M2 + obsP2id.Mag2());
		  TLorentzVector oP1id(obsP1id, E1id);						// Fill observed Lorentz vectors
		  TLorentzVector oP2id(obsP2id, E2id);
		  TLorentzVector oPtotid = oP1id + oP2id;					// Total momentum 4-vector
		  Float_t MiObsid = oPtotid.M();							// Invariant mass
		  oMinv_id->Fill(MiObsid);								// Fill invariant mass histogram

		  //
		  // Fill CLD Histograms
		  ObsTrk *Tr1cl = new ObsTrk(tX, tP1, Q1, Bfield, GCcl);	// Apply track resolution
		  ObsTrk *Tr2cl = new ObsTrk(tX, tP2, Q2, Bfield, GCcl);
		  TVector3 obsP1cl = Tr1cl->GetObsP();					// Get observed momenta
		  TVector3 obsP2cl = Tr2cl->GetObsP();
		  Double_t E1cl = TMath::Sqrt(M1*M1 + obsP1cl.Mag2());	// Observed energies
		  Double_t E2cl = TMath::Sqrt(M2*M2 + obsP2cl.Mag2());
		  TLorentzVector oP1cl(obsP1cl, E1cl);						// Fill observed Lorentz vectors
		  TLorentzVector oP2cl(obsP2cl, E2cl);
		  TLorentzVector oPtotcl = oP1cl + oP2cl;					// Total momentum 4-vector
		  Float_t MiObscl = oPtotcl.M();							// Invariant mass
		  oMinv_cl->Fill(MiObscl);									// Fill invariant mass histogram
	  }
   }
   //
   fin->Close();
   //
   // Plot histograms
   TCanvas* c1 = new TCanvas("c1","Pythia8 HZ (H{#rightarrow}#mu #mu",10,10,800,800);
   c1->Divide(2, 2);
   c1->cd(1);
   //etaMu->Scale(5./Float_t(nev));
   etaMu->Draw();
   etaMu->SetXTitle("#eta");
   etaMu->SetYTitle("dN/d#eta");
   //
   c1->cd(2);
   //gPad->SetLogy();
   //ptMu->Scale(5./Float_t(nev));
   ptMu->Draw();
   ptMu->SetXTitle("p_{t} [GeV/c]");
   ptMu->SetYTitle("dN/dp_{t}^{2} [GeV/c]^{-2}");
   //
   c1->cd(3);
   // Invariant mass plot
   gMinv->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV)");
   gMinv->Draw();
   c1->cd(4);
   // Invariant mass plot
   oMinv_id->SetXTitle("#mu^{+}#mu^{-} invariant mass (GeV)");
   oMinv_id->SetLineColor(kBlue);
   oMinv_id->Draw();
   oMinv_cl->SetLineColor(kRed);
   oMinv_cl->Draw("SAME");
   //
   TLegend *lg = new TLegend(0.1, 0.9, 0.3, 0.70);
   TString LgTitle = "Detectors:";
   lg->SetHeader(LgTitle);
   lg->AddEntry(oMinv_id, "IDEA", "L");
   lg->AddEntry(oMinv_cl, "CLD", "L");
   lg->Draw();
 }
