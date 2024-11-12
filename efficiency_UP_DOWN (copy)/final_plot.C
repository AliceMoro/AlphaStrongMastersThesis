using namespace std;

void estimate_as() {
	
	TCanvas *c = new TCanvas("c","sigma_as");
	
	// experimental dijet cross section from dijet sample with JES and JER corrections
	
	double sigma_data = 1.30698e+07;	 // pb
	double err_sigma_data = 1.07338e+06; // pb
	double err_stat = 85658;			 // pb
	double err_eff = 1.03753e+06;		 // pb
	double err_L = 261396;				 // pb
	
	TF1 *gdata = new TF1("gdata","[0]",0.10,0.14);
	TF1 *gdataUP = new TF1("gdataUP","[1]+[2]",0.10,0.14);
	TF1 *gdataDOWN = new TF1("gdataDOWN","[3]-[4]",0.10,0.14);
	gdata->SetParameter(0,sigma_data);
	gdataUP->SetParameter(1,sigma_data);
	gdataUP->SetParameter(2,err_sigma_data);
	gdataDOWN->SetParameter(3,sigma_data);
	gdataDOWN->SetParameter(4,err_sigma_data);
	
	// simulated values of dijet cross section weighted by acceptance
	
	TGraphErrors *gsim = new TGraphErrors("as_sigma.txt","%lg %lg %lg");
	TGraphErrors *gSTAR = new TGraphErrors("as_STAR.txt","%lg %lg");
	
	// linear fit of simulated cross sections
	
	TF1 *fit = new TF1("Linear fit of sigma_sim","[0]*x+[1]",0.10,0.14);
	
	fit->SetParNames("m","q");
	
	gsim->Fit(fit);
	
	double m = fit->GetParameter(0);
	double q = fit->GetParameter(1);
	double err_m = fit->GetParError(0);
	double err_q = fit->GetParError(1);
	
	TF1 *fitplot = new TF1("Linear fit of sigma_sim","[0]*x+[1]",0.10,0.14);
	fitplot->SetParameter(0,m);
	fitplot->SetParameter(1,q);
	
	// intersection between experimental cross section and fit line gives αs estimate
	
	double asSTAR = (sigma_data - q) / m; // αs estimate
	
	// intersections between fit line and limits of experimental individual error bars (vertical) gives the corresponding estimate of total error bar (horizontal) of αs
	
	double as_statUP = (sigma_data + err_stat - q) / m;
	double as_statDOWN = (sigma_data - err_stat - q) / m;
	
	double Delta_as_stat = as_statUP - as_statDOWN;
	
	cout << "\nΔαs*(stat.) = " << Delta_as_stat << endl;
	
	double as_effUP = (sigma_data + err_eff - q) / m;
	double as_effDOWN = (sigma_data - err_eff - q) / m;
	
	double Delta_as_eff = as_effUP - as_effDOWN;
	
	cout << "\nΔαs*(syst.) = " << Delta_as_eff << endl;
	
	double as_LUP = (sigma_data + err_L - q) / m;
	double as_LDOWN = (sigma_data - err_L - q) / m;
	
	double Delta_as_L = as_LUP - as_LDOWN;
	
	cout << "\nΔαs*(lumi.) = " << Delta_as_L << endl;
	
	// intersections between fit line and limits of experimental error bar (vertical) gives the corresponding estimate of total error bar (horizontal) of αs
	
	double asUP = (sigma_data + err_sigma_data - q) / m;
	double asDOWN = (sigma_data - err_sigma_data - q) / m;
	
	double Delta_as = asUP - asDOWN;
	
	cout << "\n αs* = " << asSTAR << "\nΔαs* = " << Delta_as << "\n" << endl;
	
	// plot
	
	gdata->SetTitle("#sigma_{jj}   vs #alpha_{s}(m^{2}_{Z^{0}});#alpha_{s}(m^{2}_{Z^{0}});#sigma_{jj}");
	gdata->GetYaxis()->SetRangeUser(1e7,2.2e7);
	gdata->SetLineColor(kRed);
	gdataUP->SetLineColor(kRed);
	gdataUP->SetLineStyle(7);
	gdataDOWN->SetLineColor(kRed);
	gdataDOWN->SetLineStyle(7);
	
	gsim->SetMarkerColor(kBlue+2);
	gsim->SetMarkerStyle(20);
	gsim->SetMarkerSize(1.8);
	
	gSTAR->SetMarkerColor(kGreen+2);
	gSTAR->SetMarkerStyle(29);
	gSTAR->SetMarkerSize(3.6);
	
	fitplot->SetLineColor(kBlue+2);
	
//	TLine *xline = new TLine(asSTAR,1e7,asSTAR,2.2e7);
	
//	xline->SetLineColor(kGreen+2);
//	xline->SetLineStyle(7);
	
	gStyle->SetTitleFontSize(0.038);
	gStyle->SetOptFit(111);
	
	c->cd();
	
	gdata->Draw();
	gdataUP->Draw("SAME");
	gdataDOWN->Draw("SAME");
	fitplot->Draw("SAME");
	gsim->Draw("P SAME");
	gSTAR->Draw("P SAME");
//	xline->Draw("L SAME");
	
	TLegend leg(.1,.62,.42,.9); // (x1,y1,x2,y2)
	
	leg.AddEntry(gdata,"Exp. cross section","L");
	leg.AddEntry(gsim,"Sim. cross sections","EP");
	leg.AddEntry(fitplot,"Fit line of sim. cross sections","L");
	leg.AddEntry(gSTAR,"#alpha_{s}(m^{2}_{Z^{0}}) from this work","P");
	leg.DrawClone();
	
}
