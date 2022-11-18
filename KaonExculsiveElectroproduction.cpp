#include"KaonExculsiveElectroproduction.h"

#include<iostream>
#include<string.h>
#include<cmath>
#include"TMath.h"

using namespace std;


KaonExculsiveElectroproduction::KaonExculsiveElectroproduction(){
	cout<<"****EicC Meson Form Factor Project"<<endl;
	cout<<"****Coding issues, contact rwang@impcas.ac.cn"<<endl;
	cout<<endl<<endl;
	cout<<"    Simulation starting..."<<endl<<endl;
	cout<<"    Exclusive process: e p --> e Lambda K^+"<<endl;
	cout<<"    Both decays of "<<endl;
	cout<<"        \"Lambda --> p pi^-\" "<<endl;
	cout<<"    and \"Lambda --> n pi^0, pi^0 --> 2 gammas\" are implemented."<<endl;
	cout<<endl;

	beam_cross_angle = 0.05; //// 50 mrad = 0.05rad
	///// the kinematical ranges for MC sampling
	xBmin = 0.0001;
	xBmax = 0.95;
	Q2min = 1;
	Q2max = 50;
	t0 = -0.001;
	t1 = -100;
	ymin = 0.01;
	ymax = 0.99;
	Tmin = 0.01;
	Tmax = 50;
	//// nucleon mass, electron mass and others...
	mN = 0.938272;
	mpi = 0.13957039;
	me = 0.000511;
	mkaon = 0.493677;
	mLambda = 1.115683;
	mp = 0.938272088;
	mn = 0.9395654205;
	pi0mass = 0.1349768;


	cutoff_K = 0.8; //GeV

	PI = TMath::Pi();  //// 3.141592653;
	alpha = 1.0/137.0;

	strFileName = new char[100];
	strcpy(strFileName, "DEMP_pion.root");

	//// EicC optimal collision energy
	eBeamE = 3.5;
	pBeamE = 20;
	eBeam = new TLorentzVector(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam = new TLorentzVector(p_px, 0, p_pz, pBeamE);
	BoostToEIC = new TVector3(0,0,0);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
	eBeam->Boost(*BoostToEIC);
	
	//// initialize the final-state particles
	elec_out = new TLorentzVector(0, 0, 0, 0);
	Lambda_out = new TLorentzVector(0, 0, 0, 0);
	Kp_out = new TLorentzVector(0, 0, 0, 0);
	LambdaDecayProt_out = new TLorentzVector(0, 0, 0, 0);
	LambdaDecayPim_out = new TLorentzVector(0, 0, 0, 0);
	LambdaDecayNeut_out = new TLorentzVector(0, 0, 0, 0);
	LambdaDecayPi0_out = new TLorentzVector(0, 0, 0, 0);
	LambdaDecayPi0Gamma1_out = new TLorentzVector(0, 0, 0, 0);
	LambdaDecayPi0Gamma2_out = new TLorentzVector(0, 0, 0, 0);
	//// initialize the Lambda decay vertex
	LambdaDecayVx = 0;
	LambdaDecayVy = 0;
	LambdaDecayVz = 0;

	/// random seed
	random = new TRandom3(0);
	/// kinematic calculator
	kine = new KineCal();
	/// phase space generator
	eventGenerator = new TGenPhaseSpace();
}
KaonExculsiveElectroproduction::~KaonExculsiveElectroproduction(){
	tree->Write();
	fout->Write();
	fout->Close();
	cout<<"    Data file saved and closed~"<<endl<<endl;
	//delete random;
	//delete tree;
	//delete fout;
	//delete kine;
	//delete eBeam;
	//delete pBeam;
	//delete elec_out;
	//delete Lambda_out;
}

int KaonExculsiveElectroproduction::Generate(int N = 20000){

	MakeROOTFile(strFileName);

	cout<<"    To generate "<<N<<" events..."<<endl;

	eBeam->Boost(-(*BoostToEIC));   /// boost to the proton-rest frame

	for(int i=0; i<N; ){
		xB = random->Uniform(xBmin, xBmax);
		Q2 = random->Uniform(Q2min, Q2max);
		if(Q2>(xB*ymax*(s-mN*mN)))continue;
		W2 = (1.0/xB-1)*Q2 + mN*mN;
		//if(W2<2.6)continue;
		if(W2<4.0)continue;
		y = Q2 / xB / (s-mN*mN);
		epsilon = kine->calEpsilon(y, Q2, s);
		// get the physical range of t and compare to [-Tmax, -Tmin]
		t0 = kine->calTMin(W2, -Q2, mN*mN, mkaon*mkaon, mN*mN);
		if(!(t0==t0))continue;
		if(t0<(-Tmax))continue;
		if(t0>(-Tmin))t0 = -Tmin;
		t1 = kine->calTMax(W2, -Q2, mN*mN, mkaon*mkaon, mN*mN);
		if(!(t1==t1))continue;
		if(t1>(-Tmin))continue;
		if(t1<(-Tmax))t1 = -Tmax;
		//cout<<"t0 "<<t0<<endl<<t1<<endl;  //test code.
		t = random->Uniform(t1, t0);		

		// calculate the kinematics of the final-state electron in the proton-rest frame
		double nv = Q2/2.0/mN/xB;
		double eOutE = eBeamENRest - nv;
		//cout<<eBeamENRest<<"\t"<<nv<<endl;   // test code.
		double costheta_e = 1 - Q2/2.0/eBeamENRest/eOutE;
		double sintheta_e = sqrt(1 - costheta_e*costheta_e);
		double ephi = random->Uniform(0, 2*PI);
		double emom = sqrt(eOutE*eOutE - me*me);
		elec_out->SetXYZT(emom*sintheta_e*cos(ephi), emom*sintheta_e*sin(ephi), emom*costheta_e, eOutE);
		// calculate the kinematics of the final-state Lambda and kaon^+ in the proton-rest frame
		double ELambda = (mN*mN+mLambda*mLambda-t)/2.0/mN;
		double EKp = nv + mN - ELambda;
		double PLambda = sqrt(ELambda*ELambda - mLambda*mLambda);
		double PKp = sqrt(EKp*EKp - mkaon*mkaon);
		TLorentzVector virtualPhoton = (*eBeam) - (*elec_out);
		TVector3 elec_out_v3 = elec_out->Vect(); 
		TVector3 virtualPhoton_v3 = virtualPhoton.Vect();
		TVector3 normal_v3 = virtualPhoton_v3.Cross(elec_out_v3);
		TVector3 leptonplane_transverse_v3 = normal_v3.Cross(virtualPhoton_v3);
		//// check the kinematics
		if(virtualPhoton_v3.Mag()>=(PLambda+PKp))continue;
		if(virtualPhoton_v3.Mag()<=fabs(PLambda-PKp))continue;
		double costheta_Kp = (virtualPhoton_v3.Mag2()+PKp*PKp-PLambda*PLambda) /2.0 /virtualPhoton_v3.Mag() /PKp;
		double costheta_Lambda = (virtualPhoton_v3.Mag2()+PLambda*PLambda-PKp*PKp) /2.0 /virtualPhoton_v3.Mag() /PLambda;
		double sintheta_Kp = sqrt(1 - costheta_Kp*costheta_Kp);
		double sintheta_Lambda = sqrt(1 - costheta_Lambda*costheta_Lambda);
		double Phi = random->Uniform(0, 2*PI); //the angle between the hadronic plane and the leptonic plane
		TVector3 Kp_out_v3 = (PKp*costheta_Kp)*(virtualPhoton_v3.Unit());
		Kp_out_v3 = Kp_out_v3 + (PKp*sintheta_Kp*cos(Phi))*(leptonplane_transverse_v3.Unit());
		Kp_out_v3 = Kp_out_v3 + (PKp*sintheta_Kp*sin(Phi))*(normal_v3.Unit());
		TVector3 Lambda_out_v3 = (PLambda*costheta_Lambda)*(virtualPhoton_v3.Unit());
		Lambda_out_v3 = Lambda_out_v3 + (PLambda*sintheta_Lambda*cos(Phi+PI))*(leptonplane_transverse_v3.Unit());
		Lambda_out_v3 = Lambda_out_v3 + (PLambda*sintheta_Lambda*sin(Phi+PI))*(normal_v3.Unit());
		Kp_out->SetXYZT(Kp_out_v3.Px(), Kp_out_v3.Py(), Kp_out_v3.Pz(), EKp);
		Lambda_out->SetXYZT(Lambda_out_v3.Px(), Lambda_out_v3.Py(), Lambda_out_v3.Pz(), ELambda);
		//cout<<costheta_Kp<<"  "<<costheta_Lambda<<"  "<<Kp_out_v3.Pz()<<"  "<<Lambda_out_v3.Pz()<<endl;  //test code

		//// Boost from proton rest frame to collider frame!
		elec_out->Boost(*BoostToEIC);
		Kp_out->Boost(*BoostToEIC);
		Lambda_out->Boost(*BoostToEIC);

		//// Now performing the Lambda decay
		//// Sampling the decay vertex of the Lambda
		double ctau = 0.07896;  //// c*rest_frame_life, get from PDG
 		double length = Lambda_out->Beta() * Lambda_out->Gamma() * ctau; // mean fly path, in unit of meter
		double lengthMC = random->Exp(length);  //Monte-Carlo sampling
		//// in the unit of meter
 		LambdaDecayVx = lengthMC *  sin(Lambda_out->Theta() )   * cos(Lambda_out->Phi() );
 		LambdaDecayVy = lengthMC *  sin(Lambda_out->Theta() )   * sin(Lambda_out->Phi() );
 		LambdaDecayVz = lengthMC *  cos(Lambda_out->Theta() ) ;
		//// Lambda decay to proton and pi^-
		double masses[2] = {mp, mpi};
		eventGenerator->SetDecay(*Lambda_out, 2, masses);
		eventGenerator->Generate();
		*LambdaDecayProt_out = (*(eventGenerator->GetDecay(0)));
		*LambdaDecayPim_out = (*(eventGenerator->GetDecay(1)));
		//// Lambda decay to neutron and pi0, pi0-->2gammas
		masses[0]=mn;  masses[1]=pi0mass;
		eventGenerator->SetDecay(*Lambda_out, 2, masses);
		eventGenerator->Generate();
		*LambdaDecayNeut_out = (*(eventGenerator->GetDecay(0)));
		*LambdaDecayPi0_out = (*(eventGenerator->GetDecay(1)));
		masses[0]=0;  masses[1]=0;
		eventGenerator->SetDecay(*LambdaDecayPi0_out, 2, masses);
		eventGenerator->Generate();
		*LambdaDecayPi0Gamma1_out = (*(eventGenerator->GetDecay(0)));
		*LambdaDecayPi0Gamma2_out = (*(eventGenerator->GetDecay(1)));


		//// calculate the differential cross sections
		d4sigma = d4sigma_dQ2dxBdtdPhi(Q2, xB, t, 0);
		d3sigma = d3sigma_dQ2dxBdt(Q2, xB, t);

		tree->Fill();
		i++;
	}

	eBeam->Boost(*BoostToEIC);   ///the elec. beam boost back to the collider frame!!!

	cout<<"    Event generation done! "<<endl;
	return N;
}

//// Create a ROOT file and a TTree.
void KaonExculsiveElectroproduction::MakeROOTFile(char *filename){
	//// create the output file and the output TTree
	cout<<"    Creating the output file: "<<filename<<endl;
	fout = new TFile(filename,"recreate");
	tree = new TTree("tree","Exclusive kaon electroproduction");
	tree->Branch("xB", &xB, "xB/D");
	tree->Branch("Q2", &Q2, "Q2/D");
	tree->Branch("t", &t, "t/D");
	tree->Branch("y", &y, "y/D");
	tree->Branch("W2", &W2, "W2/D");
	tree->Branch("epsilon", &epsilon, "epsilon/D");
	tree->Branch("s", &s, "s/D");
	tree->Branch("d4sigma", &d4sigma, "d4sigma/D");
	tree->Branch("d3sigma", &d3sigma, "d3sigma/D");
	tree->Branch("elec_out", "TLorentzVector", elec_out);
	tree->Branch("Lambda_out", "TLorentzVector", Lambda_out);
	tree->Branch("Kp_out" , "TLorentzVector", Kp_out);
	tree->Branch("LambdaDecayProt_out" , "TLorentzVector", LambdaDecayProt_out);
	tree->Branch("LambdaDecayPim_out" , "TLorentzVector", LambdaDecayPim_out);
	tree->Branch("LambdaDecayNeut_out" , "TLorentzVector", LambdaDecayNeut_out);
	tree->Branch("LambdaDecayPi0_out" , "TLorentzVector", LambdaDecayPi0_out);
	tree->Branch("LambdaDecayPi0Gamma1_out" , "TLorentzVector", LambdaDecayPi0Gamma1_out);
	tree->Branch("LambdaDecayPi0Gamma2_out" , "TLorentzVector", LambdaDecayPi0Gamma2_out);
	tree->Branch("LambdaDecayVx", &LambdaDecayVx, "LambdaDecayVx/D");
	tree->Branch("LambdaDecayVy", &LambdaDecayVy, "LambdaDecayVy/D");
	tree->Branch("LambdaDecayVz", &LambdaDecayVz, "LambdaDecayVz/D");
}
void KaonExculsiveElectroproduction::SetOutputFileName(char *filename){
	strcpy(strFileName, filename);
}

void KaonExculsiveElectroproduction::SetElecBeamEnergy(double ebeamenergy){
	if(ebeamenergy<0.001){cout<<"Error: electron beam energy is too small!!!"<<endl; return;}
	if(ebeamenergy>1e6){cout<<"Error: electron beam energy is too high!!!"<<endl; return;}
	eBeamE = ebeamenergy;
	eBeam->SetXYZT(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam->SetXYZT(p_px, 0, p_pz, pBeamE);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
	eBeam->Boost(*BoostToEIC);
}
double KaonExculsiveElectroproduction::GetElecBeamEnergy(){return eBeamE;}
void KaonExculsiveElectroproduction::SetProtBeamEnergy(double pbeamenergy){
	if(pbeamenergy<1){cout<<"Error: proton beam energy is too small!!!"<<endl; return;}
	if(pbeamenergy>1e6){cout<<"Error: proton beam energy is too high!!!"<<endl; return;}
	pBeamE = pbeamenergy;
	eBeam->SetXYZT(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam->SetXYZT(p_px, 0, p_pz, pBeamE);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
	eBeam->Boost(*BoostToEIC);
}
double KaonExculsiveElectroproduction::GetProtBeamEnergy(){return pBeamE;}
//// set beam crossing angle
void KaonExculsiveElectroproduction::SetBeamCrossAngle(double _angle){
	beam_cross_angle = _angle;
	eBeam->SetXYZT(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam->SetXYZT(p_px, 0, p_pz, pBeamE);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
	eBeam->Boost(*BoostToEIC);
}
double KaonExculsiveElectroproduction::GetBeamCrossAngle(){return beam_cross_angle;}





//// the model of kaon form factor
double KaonExculsiveElectroproduction::FF_kaon(double _Q2){
	return 1.0 / (1.0 + _Q2/cutoff_K/cutoff_K);
}
//// return cross-section in the unit of nb/GeV^4.
double KaonExculsiveElectroproduction::d4sigma_dQ2dxBdtdPhi(double _Q2, double _xB, double _t, double _Phi)
{
	//// No implementation of the Phi-dependence yet
	return d3sigma_dQ2dxBdt(_Q2, _xB, _t) / 2.0 / PI;
}
//// return cross-section in the unit of nb/GeV^4.
double KaonExculsiveElectroproduction::d3sigma_dQ2dxBdt(double _Q2, double _xB, double _t)
{
	double flux = PhotonFlux(y, _xB, epsilon, _Q2); 
	/// in the natural unit, i.e., GeV^-6
	double sigma = flux * (dsigmaT() + epsilon*dsigmaL());
	/// return in the unit of nb/GeV^4
	return 3.8809e5 * sigma;
}
double KaonExculsiveElectroproduction::N_factor(double _W2, double _Q2){
	double W2_mN2 = _W2 - mN*mN;
	return 32*PI*W2_mN2 * sqrt(W2_mN2*W2_mN2 + _Q2*_Q2 + 2*_Q2*(W2+mN*mN) );
}
double KaonExculsiveElectroproduction::g_KNN(double _t){
	return 13.4*1.8/sqrt(3) * (cutoff_K*cutoff_K - mkaon*mkaon) / (cutoff_K*cutoff_K - _t);
	//return 14.688 * (cutoff_K*cutoff_K - mkaon*mkaon) / (cutoff_K*cutoff_K - _t);
}
double KaonExculsiveElectroproduction::PhotonFlux(double _y, double _xB, double _epsilon, double _Q2){
	return alpha*_y*_y*(1-_xB)/2.0/PI/_xB/(1-_epsilon)/_Q2; 
}
double KaonExculsiveElectroproduction::dsigmaT(){
	return 0; //// the transverse component can be ignored at very small |t| and high Q2
}
double KaonExculsiveElectroproduction::dsigmaL(){
	double nfactor = N_factor(W2, Q2);
	double ffkaon = FF_kaon(Q2);
	double gKNN = g_KNN(t);
	double Kpole = -t / (t-mkaon*mkaon) / (t-mkaon*mkaon);
	//// this is the Born-term contribution of kaon pole, valid at small |t|.
	return 16*PI*alpha * gKNN*gKNN * Kpole * Q2 * ffkaon*ffkaon / nfactor;
}




//// set sampling ranges
void KaonExculsiveElectroproduction::SetxBmin(double min){xBmin = min;}
void KaonExculsiveElectroproduction::SetxBmax(double max){xBmax = max;}
void KaonExculsiveElectroproduction::SetQ2min(double min){Q2min = min;}
void KaonExculsiveElectroproduction::SetQ2max(double max){Q2max = max;}
void KaonExculsiveElectroproduction::SetTmin(double min){Tmin = min;}
void KaonExculsiveElectroproduction::SetTmax(double max){Tmax = max;}
void KaonExculsiveElectroproduction::Setymin(double min){ymin = min;}
void KaonExculsiveElectroproduction::Setymax(double max){ymax = max;}


