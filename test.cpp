#include"KineCal.h"
#include"KaonExculsiveElectroproduction.h"

#include"KineCal.cpp"
#include"KaonExculsiveElectroproduction.cpp"


int test(){

	KaonExculsiveElectroproduction demp_kaon;

	demp_kaon.SetTmax(0.5);
	demp_kaon.SetTmin(0.01);
	demp_kaon.SetQ2max(50);
	demp_kaon.SetQ2min(1);
	demp_kaon.SetxBmax(0.8);
	demp_kaon.SetxBmin(0.001);

	char filename[50] = "DEMP-kaon-pole-at-EicC.root";
	demp_kaon.SetOutputFileName(filename); 

	demp_kaon.SetElecBeamEnergy(3.5);
	demp_kaon.SetProtBeamEnergy(20);
	//demp_kaon.SetBeamCrossAngle(0.0);
	demp_kaon.SetBeamCrossAngle(0.05);

	demp_kaon.Generate(800000);

	//cout<<demp_kaon.GetElecBeamEnergy()<<"  "<<demp_kaon.GetProtBeamEnergy()<<"  "<<demp_kaon.GetBeamCrossAngle()<<endl;

	return 0;
}


