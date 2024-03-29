{
	/*
 	  This code is for plotting the kinematics of 
  	  the output file NOT in the sampling mode!!!
  	*/

	gStyle->SetOptStat(0);
	gStyle->SetNdivisions(505,"XY");


	TFile file("DEMP-kaon-pole-at-EicC.root");

/*
	TCanvas c1("c1","Q2 vs xB, unweighted",700,550);
	tree->Draw("Q2:xB>>hQ2xB","","colz");
	gPad->SetLogz();
	hQ2xB->SetTitle("");
	hQ2xB->GetXaxis()->SetTitle("x_{B}");
	hQ2xB->GetXaxis()->SetTitleSize(0.06);
	hQ2xB->GetXaxis()->CenterTitle();
	hQ2xB->GetXaxis()->SetTitleOffset(1.05);
	hQ2xB->GetXaxis()->SetLabelSize(0.06);
	hQ2xB->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
	hQ2xB->GetYaxis()->SetTitleSize(0.06);
	hQ2xB->GetYaxis()->CenterTitle();
	hQ2xB->GetYaxis()->SetTitleOffset(1.05);
	hQ2xB->GetYaxis()->SetLabelSize(0.06);
	hQ2xB->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c2("c2","-t vs xB, unweighted",700,550);
	tree->Draw("-t:xB>>htxB","","colz");
	gPad->SetLogz();
	htxB->SetTitle("");
	htxB->GetXaxis()->SetTitle("x_{B}");
	htxB->GetXaxis()->SetTitleSize(0.06);
	htxB->GetXaxis()->CenterTitle();
	htxB->GetXaxis()->SetTitleOffset(1.05);
	htxB->GetXaxis()->SetLabelSize(0.06);
	htxB->GetYaxis()->SetTitle("-t (GeV^{2})");
	htxB->GetYaxis()->SetTitleSize(0.06);
	htxB->GetYaxis()->CenterTitle();
	htxB->GetYaxis()->SetTitleOffset(1.05);
	htxB->GetYaxis()->SetLabelSize(0.06);
	htxB->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c3("c3","Q2 vs xB, xsec. weighted",700,550);
	tree->Draw("Q2:xB>>hQ2xB_2","d4sigma*1000","colz");
	gPad->SetLogz();
	hQ2xB_2->SetTitle("");
	hQ2xB_2->GetXaxis()->SetTitle("x_{B}");
	hQ2xB_2->GetXaxis()->SetTitleSize(0.06);
	hQ2xB_2->GetXaxis()->CenterTitle();
	hQ2xB_2->GetXaxis()->SetTitleOffset(1.05);
	hQ2xB_2->GetXaxis()->SetLabelSize(0.06);
	hQ2xB_2->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
	hQ2xB_2->GetYaxis()->SetTitleSize(0.06);
	hQ2xB_2->GetYaxis()->CenterTitle();
	hQ2xB_2->GetYaxis()->SetTitleOffset(1.05);
	hQ2xB_2->GetYaxis()->SetLabelSize(0.06);
	hQ2xB_2->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c4("c4","-t vs xB, xsec. weighted",700,550);
	tree->Draw("-t:xB>>htxB_2","d4sigma*1000","colz");
	gPad->SetLogz();
	htxB_2->SetTitle("");
	htxB_2->GetXaxis()->SetTitle("x_{B}");
	htxB_2->GetXaxis()->SetTitleSize(0.06);
	htxB_2->GetXaxis()->CenterTitle();
	htxB_2->GetXaxis()->SetTitleOffset(1.05);
	htxB_2->GetXaxis()->SetLabelSize(0.06);
	htxB_2->GetYaxis()->SetTitle("-t (GeV^{2})");
	htxB_2->GetYaxis()->SetTitleSize(0.06);
	htxB_2->GetYaxis()->CenterTitle();
	htxB_2->GetYaxis()->SetTitleOffset(1.05);
	htxB_2->GetYaxis()->SetLabelSize(0.06);
	htxB_2->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);
*/



	TCanvas c5("c5","electron P vs theta, xsec. weighted",700,550);
	tree->Draw("elec_out.P():elec_out.Theta()*TMath::RadToDeg()>>helec","d4sigma*10000","colz");
	gPad->SetLogz();
	helec->SetTitle("");
	helec->GetXaxis()->SetTitle("#theta_e (#circ)");
	helec->GetXaxis()->SetTitleSize(0.06);
	helec->GetXaxis()->CenterTitle();
	helec->GetXaxis()->SetTitleOffset(1.05);
	helec->GetXaxis()->SetLabelSize(0.06);
	helec->GetYaxis()->SetTitle("P_e (GeV/c)");
	helec->GetYaxis()->SetTitleSize(0.06);
	helec->GetYaxis()->CenterTitle();
	helec->GetYaxis()->SetTitleOffset(1.05);
	helec->GetYaxis()->SetLabelSize(0.06);
	helec->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c6("c6","Kaon^+, P vs theta, xsec. weighted",700,550);
	tree->Draw("Kp_out.P():Kp_out.Theta()*TMath::RadToDeg()>>hKp","d4sigma*10000","colz");
	gPad->SetLogz();
	hKp->SetTitle("");
	hKp->GetXaxis()->SetTitle("#theta_K+ (#circ)");
	hKp->GetXaxis()->SetTitleSize(0.06);
	hKp->GetXaxis()->CenterTitle();
	hKp->GetXaxis()->SetTitleOffset(1.05);
	hKp->GetXaxis()->SetLabelSize(0.06);
	hKp->GetYaxis()->SetTitle("P_K+ (GeV/c)");
	hKp->GetYaxis()->SetTitleSize(0.06);
	hKp->GetYaxis()->CenterTitle();
	hKp->GetYaxis()->SetTitleOffset(1.05);
	hKp->GetYaxis()->SetLabelSize(0.06);
	hKp->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c7("c7","Lambda, P vs theta, xsec. weighted",700,550);
	tree->Draw("Lambda_out.P():Lambda_out.Theta()*TMath::RadToDeg()>>hLambda","d4sigma*10000","colz");
	gPad->SetLogz();
	hLambda->SetTitle("");
	hLambda->GetXaxis()->SetTitle("#theta _ #Lambda (#circ)");
	hLambda->GetXaxis()->SetTitleSize(0.06);
	hLambda->GetXaxis()->CenterTitle();
	hLambda->GetXaxis()->SetTitleOffset(1.05);
	hLambda->GetXaxis()->SetLabelSize(0.06);
	hLambda->GetYaxis()->SetTitle("P_#Lambda (GeV/c)");
	hLambda->GetYaxis()->SetTitleSize(0.06);
	hLambda->GetYaxis()->CenterTitle();
	hLambda->GetYaxis()->SetTitleOffset(1.05);
	hLambda->GetYaxis()->SetLabelSize(0.06);
	hLambda->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c8("c8","LambdaDecayProton, P vs theta, xsec. weighted",700,550);
	tree->Draw("LambdaDecayProt_out.P():LambdaDecayProt_out.Theta()*TMath::RadToDeg()>>hLambdaDecayProt","d4sigma*10000","colz");
	gPad->SetLogz();
	hLambdaDecayProt->SetTitle("");
	hLambdaDecayProt->GetXaxis()->SetTitle("#theta_p_from_ #Lambda_decay (#circ)");
	hLambdaDecayProt->GetXaxis()->SetTitleSize(0.06);
	hLambdaDecayProt->GetXaxis()->CenterTitle();
	hLambdaDecayProt->GetXaxis()->SetTitleOffset(1.05);
	hLambdaDecayProt->GetXaxis()->SetLabelSize(0.06);
	hLambdaDecayProt->GetYaxis()->SetTitle("P_p_from_ #Lambda_decay (GeV/c)");
	hLambdaDecayProt->GetYaxis()->SetTitleSize(0.06);
	hLambdaDecayProt->GetYaxis()->CenterTitle();
	hLambdaDecayProt->GetYaxis()->SetTitleOffset(1.05);
	hLambdaDecayProt->GetYaxis()->SetLabelSize(0.06);
	hLambdaDecayProt->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c9("c9","LambdaDecayPim, P vs theta, xsec. weighted",700,550);
	tree->Draw("LambdaDecayPim_out.P():LambdaDecayPim_out.Theta()*TMath::RadToDeg()>>hLambdaDecayPim","d4sigma*10000","colz");
	gPad->SetLogz();
	hLambdaDecayPim->SetTitle("");
	hLambdaDecayPim->GetXaxis()->SetTitle("#theta_ #pi-_from_ #Lambda_decay (#circ)");
	hLambdaDecayPim->GetXaxis()->SetTitleSize(0.06);
	hLambdaDecayPim->GetXaxis()->CenterTitle();
	hLambdaDecayPim->GetXaxis()->SetTitleOffset(1.05);
	hLambdaDecayPim->GetXaxis()->SetLabelSize(0.06);
	hLambdaDecayPim->GetYaxis()->SetTitle("P_ #pi-_from_ #Lambda_decay (GeV/c)");
	hLambdaDecayPim->GetYaxis()->SetTitleSize(0.06);
	hLambdaDecayPim->GetYaxis()->CenterTitle();
	hLambdaDecayPim->GetYaxis()->SetTitleOffset(1.05);
	hLambdaDecayPim->GetYaxis()->SetLabelSize(0.06);
	hLambdaDecayPim->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);

	TCanvas c10("c10","Lambda decay vertex, Vx vs Vz, xsec. weighted",700,550);
	tree->Draw("LambdaDecayVx:LambdaDecayVz>>hvertex(50,-4,0,40,0,0.4)","d4sigma*10000","colz");
	gPad->SetLogz();
	hvertex->SetTitle("");
	hvertex->GetXaxis()->SetTitle("#Lambda decay vertex z (m)");
	hvertex->GetXaxis()->SetTitleSize(0.06);
	hvertex->GetXaxis()->CenterTitle();
	hvertex->GetXaxis()->SetTitleOffset(1.05);
	hvertex->GetXaxis()->SetLabelSize(0.06);
	hvertex->GetYaxis()->SetTitle("#Lambda decay vertex x (m)");
	hvertex->GetYaxis()->SetTitleSize(0.06);
	hvertex->GetYaxis()->CenterTitle();
	hvertex->GetYaxis()->SetTitleOffset(1.05);
	hvertex->GetYaxis()->SetLabelSize(0.06);
	hvertex->GetZaxis()->SetLabelSize(0.05);
	gPad->SetLeftMargin(0.135);
	gPad->SetBottomMargin(0.14);
	gPad->SetTopMargin(0.03);
	gPad->SetRightMargin(0.137);


}
