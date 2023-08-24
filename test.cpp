#include"KineCal.h"
#include"KaonExclusiveElectroproduction.h"

#include"KineCal.cpp"
#include"KaonExclusiveElectroproduction.cpp"


int test(){

	KaonExclusiveElectroproduction demp_kaon;

	demp_kaon.SetTmax(0.5);
	demp_kaon.SetTmin(0.01);
	demp_kaon.SetQ2max(50);
	demp_kaon.SetQ2min(1);
	demp_kaon.SetxBmax(0.8);
	demp_kaon.SetxBmin(0.001);

	//char filename[50] = "DEMP-kaon-pole-at-EicC.root";
	//TString filename = "DEMP-kaon-pole-at-EicC.root";
	demp_kaon.SetOutputFileName("DEMP-kaon-pole-at-EicC.root"); 

	demp_kaon.SetElecBeamEnergy(3.5);
	demp_kaon.SetProtBeamEnergy(20);
	//demp_kaon.SetBeamCrossAngle(0.0);
	demp_kaon.SetBeamCrossAngle(0.05);

	demp_kaon.SetSamplingMode(1);
	demp_kaon.Generate(100000);

	//cout<<demp_kaon.GetElecBeamEnergy()<<"  "<<demp_kaon.GetProtBeamEnergy()<<"  "<<demp_kaon.GetBeamCrossAngle()<<endl;

	return 0;
}


