
#include "BsDecay.h"


//
// constructor
BsDecay::BsDecay(TVector3 pBs)
{
	Double_t cSpeed = TMath::C();	// m/sec
	fNtrk = 0;
	// Add Bs
	fBsInd = fNtrk;
	fPlist.push_back(new TVector3(pBs));
	Int_t BsID = 531;
	fNlist.push_back(BsID);
	Double_t BsMass   = TDatabasePDG::Instance()->GetParticle(BsID)->Mass();
	fMlist.push_back(BsMass);
	//Double_t BsLife = TDatabasePDG::Instance()->GetParticle(BsID)->Lifetime();  // seconds
	//std::cout<<"BsLife = "<<BsLife<<std::endl;
	Double_t BsLife = 1.509e-12;
	Double_t BsEn = TMath::Sqrt(BsMass*BsMass+pBs.Mag2());
	Double_t BsGamma = BsEn/BsMass;
	Double_t BsPath = gRandom->Exp(cSpeed*BsGamma*BsLife);
	//std::cout<<"Bs gammma = "<<BsGamma<<", BsPath = "<<BsPath<<std::endl;
	fLlist.push_back(BsPath);
	fNtrk++;
	//std::cout<<"Done Bs"<<std::endl;
	// Decay Bs --> Ds- pi+
	TLorentzVector** pout = new TLorentzVector * [2];
	TLorentzVector pin(pBs,BsEn);
	Int_t DsID = -431;	// Ds-
	Double_t DsMass = TDatabasePDG::Instance()->GetParticle(DsID)->Mass();
	Int_t PiPlusID	= 211;	// Pi+
	Double_t PiMass = TDatabasePDG::Instance()->GetParticle(PiPlusID)->Mass();
	TwoBodyPS(pin, PiMass, DsMass, pout);
	// Add Pi+ from Bs
	fPiBsInd = fNtrk;
	fPlist.push_back(new TVector3(pout[0]->Vect()));
	fNlist.push_back(PiPlusID);
	fMlist.push_back(PiMass);
	fLlist.push_back(0.);
	fNtrk++;
	//std::cout<<"Done PiBs"<<std::endl;
	// Add Ds-
	fDsInd = fNtrk;
	fPlist.push_back(new TVector3(pout[1]->Vect()));
	fNlist.push_back(DsID);
	fMlist.push_back(DsMass);
	//Double_t DsLife = TDatabasePDG::Instance()->GetParticle(DsID)->Lifetime();  // seconds
	//std::cout<<"DsLife = "<<BsLife<<std::endl;
	Double_t DsLife = 0.5e-12;
	Double_t DsEn = pout[1]->E();
	Double_t DsGamma = DsEn/DsMass;
	Double_t DsPath = gRandom->Exp(cSpeed*DsGamma*DsLife);
	//std::cout<<"Ds gammma = "<<DsGamma<<", DsPath = "<<DsPath<<std::endl;
	fLlist.push_back(DsPath);
	fNtrk++;
	//std::cout<<"Done Ds"<<std::endl;
	// Decay Ds- --> phi pi-
	pin = *pout[1];
	Int_t PhiID = 333;
	Double_t PhiMass = TDatabasePDG::Instance()->GetParticle(PhiID)->Mass();
	TwoBodyPS(pin, PhiMass, PiMass, pout);
	// Add phi
	fPhiInd = fNtrk;
	fPlist.push_back(new TVector3(pout[0]->Vect()));
	fNlist.push_back(PhiID);
	fMlist.push_back(PhiMass);
	fLlist.push_back(0.);
	fNtrk++;
	//std::cout<<"Done Phi"<<std::endl;
	// Add pi- from Ds
	fPiDsInd = fNtrk;
	fPlist.push_back(new TVector3(pout[1]->Vect()));
	fNlist.push_back(-PiPlusID);
	fMlist.push_back(PiMass);
	fLlist.push_back(0.);
	fNtrk++;
	//std::cout<<"Done PiDs"<<std::endl;
	// decay Phi --> K+ K-
	pin = *pout[0];
	Int_t KplusID = 321;
	Double_t KMass = TDatabasePDG::Instance()->GetParticle(KplusID)->Mass();
	TwoBodyPS(pin, KMass, KMass, pout);
	// Add K+ from Phi
	fKpPhiInd = fNtrk;
	fPlist.push_back(new TVector3(pout[0]->Vect()));
	fNlist.push_back(KplusID);
	fMlist.push_back(KMass);
	fLlist.push_back(0.);
	//std::cout<<"Done KpPhi"<<std::endl;
	fNtrk++;
	// Add K- from Phi
	fKmPhiInd = fNtrk;
	fPlist.push_back(new TVector3(pout[1]->Vect()));
	fNlist.push_back(-KplusID);
	fMlist.push_back(KMass);
	fLlist.push_back(0.);
	fNtrk++;
	//std::cout<<"Done KmPhi"<<std::endl;
//
	delete [] pout;
//
	
}

BsDecay::~BsDecay()
{
}

void BsDecay::PlotTrack(TVectorD pr, Double_t phiMx, Int_t Npx, Double_t *x, Double_t *y)
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

TLorentzVector BsDecay::Boosted(TLorentzVector Booster, TLorentzVector InMoment)
{
	//
	// Get boost features from 4-vector
	//
	Double_t gamma = Booster.Gamma();
	Double_t beta = Booster.Beta();
	TVector3 nBoost = (1.0 / Booster.Rho()) * Booster.Vect();	// Boost versor 
	Double_t pl = nBoost.Dot(InMoment.Vect());			// Pre-boost longitudinal component
	TVector3 pt = InMoment.Vect() - pl * nBoost;			// Transverse component
	Double_t pll = gamma * (pl + beta * InMoment.E());		// Boosted logitudinal component
	TVector3 pout = pt + pll * nBoost;				// Boosted momentum
	Double_t Eout = gamma * (InMoment.E() + beta * pl);		// Boosted energy
	TLorentzVector Pboosted(pout, Eout);				// Boosted LorentzVector
	//
	return Pboosted;
}
void BsDecay::TwoBodyPS(TLorentzVector pin, Double_t m1, Double_t m2, TLorentzVector**& pout)
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
	// momentum of outgoing particles
	Double_t p = TMath::Sqrt((Min2 - (m1 - m2) * (m1 - m2)) * (Min2 - (m1 + m2) * (m1 + m2))) / (2 * Min);
	// Assume flat in phi and cos(theta)
	Double_t CosTh = 2 * gRandom->Rndm() - 1.0;
	Double_t SinTh = TMath::Sqrt(1.0 - CosTh * CosTh);
	Double_t Phi = TMath::TwoPi() * gRandom->Rndm();
	// CM momentum of part. 1
	TLorentzVector p1cm(p * SinTh * TMath::Cos(Phi), p * SinTh * TMath::Sin(Phi), p * CosTh, E1);	
	// CM momentum of part. 2
	TLorentzVector p2cm(-p * SinTh * TMath::Cos(Phi), -p * SinTh * TMath::Sin(Phi), -p * CosTh, E2);
	//
	// Boost to lab frame
	//
	pout[0] = new TLorentzVector(Boosted(pin, p1cm));
	pout[1] = new TLorentzVector(Boosted(pin, p2cm));
	//
}

