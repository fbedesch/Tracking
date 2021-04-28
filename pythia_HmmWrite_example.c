/// \file
/// \ingroup tutorial_pythia
/// pythia8 basic example
///
/// to run, do:
///
/// ~~~{.cpp}
///  root > .x pythia8.C
/// ~~~
///
/// Note that before executing this script, for some Pythia8 builds:
///
///  - the env variable PYTHIA8 must point to the pythia8100 (or newer) directory
///  - the env variable PYTHIA8DATA must be defined and it must point to $PYTHIA8/xmldoc
///
/// \macro_code
///
/// \author Andreas Morsch, modifications to write root file be F. Bedeschi
#ifdef __CLING__
R__ADD_INCLUDE_PATH($PYTHIA8 / include)
R__ADD_LIBRARY_PATH($PYTHIA8/lib)
R__ADD_LIBRARY_PATH($ROOTSYS / lib)
R__LOAD_LIBRARY(libEG)
R__LOAD_LIBRARY(libEGPythia8)
#endif

#include "TSystem.h"
#include "TROOT.h"
#include "TH1.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

void HmmWrite(Int_t nev = 100, TString fname = "PyOut.root", Int_t ndeb = 1)
{
	//
	// Load libraries
	gSystem->Load("libEG");
	gSystem->Load("libEGPythia8"); 
	//
	// Setup File and Tree for output
	TFile *fout = new TFile(fname, "RECREATE");
	TTree *T = new TTree("tID", "Pythia particles");
	// Array of particles
	TClonesArray* particles = new TClonesArray("TParticle", 1000);
	// Link to TTree
	T->Branch("Particles", "TClonesArray", &particles, 64000, 0);

	//
	// Pythia8 initializations
	TPythia8* pythia8 = new TPythia8();		// Create pythia8 object
	pythia8->ReadConfigFile("config_ee_zh_Hmumu.cmd");	// Configure Pythia generation
	pythia8->Initialize(11 /* e- */, -11 /* e+ */, 240 /* GeV */);	// Initialization

	//
	// Histograms
	TH1F* etaMu = new TH1F("etaMu", "Pseudorapidity", 120, -3., 3.);
	TH1F* ptMu = new TH1F("ptMu", "pt", 100, 0., 100.);

	// Event loop
	for (Int_t iev = 0; iev < nev; iev++) {
		pythia8->GenerateEvent();			// Event generation
		if (iev < ndeb) pythia8->EventListing();	// Event printout
		pythia8->ImportParticles(particles, "All");	// Import pythia output into a TParticle collection
		//std::cout << "Particles imported" << std::endl;
		Int_t np = particles->GetEntriesFast();		// # of TParticle
		//std::cout << "# particles "<<np << std::endl;
		// Particle loop
		Int_t Nmu = 0;
		Int_t NmuZ = 0;
		Int_t Hmth = 0;
		for (Int_t ip = 0; ip < np; ip++) {
			TParticle* part = (TParticle*)particles->At(ip);	// Get TParticle
			//std::cout << "Got TParticle. ip = "<<ip << std::endl;
			Int_t ist = part->GetStatusCode();					// Pythia status code
			//std::cout << "Status code = "<<ist << std::endl;
			// Positive codes are final state particles.
			if (ist <= 0) continue;							// skip if not final state
			Int_t pdg = part->GetPdgCode();					// PDG particle code
			Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
			if (charge == 0.) continue;						// Skip if neutral
			Float_t eta = part->Eta();						// Get eta and pt
			Float_t pt = part->Pt();
			Int_t iMuon = 13;				// Muon PDG code
			Int_t iHiggs = 25;				// Higgs PDG code
			if (TMath::Abs(pdg) == iMuon)	// Final state muon
			{
				Nmu++;						// Count muons
				//
				// Look for muons from Higgs decays
				Int_t Mother = part->GetFirstMother();	// Muon mother
				while (Mother > 0)
				{
					TParticle* mpart = (TParticle*)particles->At(Mother);
					Int_t mpdg = mpart->GetPdgCode();
					if (mpdg == iHiggs)		// Mother is Higgs
					{
						Hmth++;				// Count muons from Higgs
						etaMu->Fill(eta);	// Fill histograms
						ptMu->Fill(pt);
		//				std::cout << "eta = " << eta << std::endl;
						break;
					}
					Mother = mpart->GetFirstMother();
				}
			}

		} // end particle loop
		
		//
		// fill tree and write it
		T->Fill();
		fout->Write();

	} // end event loop

	// 
	// Wrap up
	pythia8->PrintStatistics();	// Summary printout
	fout->Close();	// close output file
	std::cout << "Output closed" << std::endl;
	delete fout;
	std::cout << "file handle deleted" << std::endl;
	//
	// Plot histograms
	TCanvas* c1 = new TCanvas("c1", "Pythia8 Z->mu mu", 800, 800);
	std::cout << "Canvas created" << std::endl;
	c1->Divide(1, 2);
	c1->cd(1);
	std::cout << "After cd(1)" << std::endl;
	//etaMu->SetXTitle("#eta");
	//etaMu->SetYTitle("dN/d#eta");
	etaMu->Draw();
	std::cout << "etaMu drawn" << std::endl;
	//
	c1->cd(2);
	std::cout << "After cd(2)" << std::endl;
	//ptMu->SetXTitle("p_{t} [GeV/c]");
	//ptMu->SetYTitle("dN/dp_{t}^{2} [GeV/c]^{-2}");
	ptMu->Draw();
	std::cout << "ptMu drawn" << std::endl;
}
