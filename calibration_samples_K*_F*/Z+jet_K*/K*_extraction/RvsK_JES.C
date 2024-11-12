#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TBranchElement.h>
#include <TBrowser.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAttAxis.h>
#include <TPad.h>
#include <TColor.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

void getRvsKbins() { // MC modified, data fixed
	
	// ******************************* //
	// * INSERT BIN LIMITS MANUALLY! * //
	// ******************************* //
	double jetpt_min = 6e4; // MeV
	double jetpt_max = 6.5e4; // MeV
	double jetpt_centr = (jetpt_min + jetpt_max) / 2;
	
	TFile *fdata = new TFile("/home/alice/Desktop/TESI/calibration_samples/Z+jet/Zmumu/Zmumu_data.root");
	TFile *fmc = new TFile("/home/alice/Desktop/TESI/calibration_samples/Z+jet/Zmumu/Zmumu_mc.root");
	
	fdata->cd("DecayTreeTuple");
	TTree *tdata = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	fmc->cd("DecayTreeTuple");
	TTree *tmc = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	int Ndata = tdata->GetEntries();
	int Nmc = tmc->GetEntries();
		
	cout << "\nNumber of entries of data tuple = " << Ndata << "\nNumber of entries of MC tuple   = " << Nmc << "\n" << endl;
	
	double jetpt_data = 0.0;
	double jeteta_data = 0.0;
	
	double Zpt_data = 0.0;
	
	double muppt_data = 0.0;
	double muppx_data = 0.0;
	double muppy_data = 0.0;
	double muppz_data = 0.0;
	
	double mumpt_data = 0.0;
	double mumpx_data = 0.0;
	double mumpy_data = 0.0;
	double mumpz_data = 0.0;
	
	double jetpt_mc = 0.0;
	double jeteta_mc = 0.0;
	
	double Zpt_mc = 0.0;
	
	double muppt_mc = 0.0;
	double muppx_mc = 0.0;
	double muppy_mc = 0.0;
	double muppz_mc = 0.0;
	
	double mumpt_mc = 0.0;
	double mumpx_mc = 0.0;
	double mumpy_mc = 0.0;
	double mumpz_mc = 0.0;
	
	tdata->SetBranchAddress("Jet_PT",&jetpt_data);
	tdata->SetBranchAddress("Jet_Eta",&jeteta_data);
	
	tdata->SetBranchAddress("Z0_PT",&Zpt_data);
	
	tdata->SetBranchAddress("mup_PT",&muppt_data);
	tdata->SetBranchAddress("mup_PX",&muppx_data);
	tdata->SetBranchAddress("mup_PY",&muppy_data);
	tdata->SetBranchAddress("mup_PZ",&muppz_data);
	
	tdata->SetBranchAddress("mum_PT",&mumpt_data);
	tdata->SetBranchAddress("mum_PX",&mumpx_data);
	tdata->SetBranchAddress("mum_PY",&mumpy_data);
	tdata->SetBranchAddress("mum_PZ",&mumpz_data);
	
	tmc->SetBranchAddress("Jet_PT",&jetpt_mc);
	tmc->SetBranchAddress("Jet_Eta",&jeteta_mc);
	
	tmc->SetBranchAddress("Z0_PT",&Zpt_mc);
	
	tmc->SetBranchAddress("mup_PT",&muppt_mc);
	tmc->SetBranchAddress("mup_PX",&muppx_mc);
	tmc->SetBranchAddress("mup_PY",&muppy_mc);
	tmc->SetBranchAddress("mup_PZ",&muppz_mc);
	
	tmc->SetBranchAddress("mum_PT",&mumpt_mc);
	tmc->SetBranchAddress("mum_PX",&mumpx_mc);
	tmc->SetBranchAddress("mum_PY",&mumpy_mc);
	tmc->SetBranchAddress("mum_PZ",&mumpz_mc);
	
	// trigger lines
	bool mupL0MuonEW_data = 0.0;
	bool mupHLT1MuonHighPt_data = 0.0;
	bool mupHLT2_data = 0.0;
	
	bool mumL0MuonEW_data = 0.0;
	bool mumHLT1MuonHighPt_data = 0.0;
	bool mumHLT2_data = 0.0;
	
	bool mupL0MuonEW_mc = 0.0;
	bool mupHLT1MuonHighPt_mc = 0.0;
	bool mupHLT2_mc = 0.0;
	
	bool mumL0MuonEW_mc = 0.0;
	bool mumHLT1MuonHighPt_mc = 0.0;
	bool mumHLT2_mc = 0.0;
	
	tdata->SetBranchAddress("mup_L0MuonEWDecision_TOS",&mupL0MuonEW_data);
	tdata->SetBranchAddress("mup_Hlt1SingleMuonHighPTDecision_TOS",&mupHLT1MuonHighPt_data);
	tdata->SetBranchAddress("mup_Hlt2Global_TOS",&mupHLT2_data);
	
	tdata->SetBranchAddress("mum_L0MuonEWDecision_TOS",&mumL0MuonEW_data);
	tdata->SetBranchAddress("mum_Hlt1SingleMuonHighPTDecision_TOS",&mumHLT1MuonHighPt_data);
	tdata->SetBranchAddress("mum_Hlt2Global_TOS",&mumHLT2_data);
	
	tmc->SetBranchAddress("mup_L0MuonEWDecision_TOS",&mupL0MuonEW_mc);
	tmc->SetBranchAddress("mup_Hlt1SingleMuonHighPTDecision_TOS",&mupHLT1MuonHighPt_mc);
	tmc->SetBranchAddress("mup_Hlt2Global_TOS",&mupHLT2_mc);
	
	tmc->SetBranchAddress("mum_L0MuonEWDecision_TOS",&mumL0MuonEW_mc);
	tmc->SetBranchAddress("mum_Hlt1SingleMuonHighPTDecision_TOS",&mumHLT1MuonHighPt_mc);
	tmc->SetBranchAddress("mum_Hlt2Global_TOS",&mumHLT2_mc);
	
	// to compute DeltaPhi
	double jetpx_data = 0.0;
	double jetpy_data = 0.0;
	double Zpx_data = 0.0;
	double Zpy_data = 0.0;
	
	double jetpx_mc = 0.0;
	double jetpy_mc = 0.0;
	double Zpx_mc = 0.0;
	double Zpy_mc = 0.0;
	
	tdata->SetBranchAddress("Jet_PX",&jetpx_data);
	tdata->SetBranchAddress("Jet_PY",&jetpy_data);
	tdata->SetBranchAddress("Z0_PX",&Zpx_data);
	tdata->SetBranchAddress("Z0_PY",&Zpy_data);
	tmc->SetBranchAddress("Jet_PX",&jetpx_mc);
	tmc->SetBranchAddress("Jet_PY",&jetpy_mc);
	tmc->SetBranchAddress("Z0_PX",&Zpx_mc);
	tmc->SetBranchAddress("Z0_PY",&Zpy_mc);
	
	// to select only events with one jet
	ULong64_t totCand_data = 0;
	ULong64_t totCand_mc = 0;
	
	tdata->SetBranchAddress("totCandidates",&totCand_data);
	tmc->SetBranchAddress("totCandidates",&totCand_mc);
	
	// arrays to store R values
	double Rdata[190849] = {0.0};
	double Rmc[374981] = {0.0};
	
	TH1D *hRdata = new TH1D("hRdata","Calibration sample Z^{0}+jet - R = p_{T}(jet) / p_{T}(Z^{0}) - data;R;A.U.",1000,0,10);
	TH1D *hRmc = new TH1D("hRmc","Calibration sample Z^{0}+jet - R = p_{T}(jet) / p_{T}(Z^{0}) - MC;R;A.U.",1000,0,10);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double phi_data = acos(((Zpx_data*jetpx_data) + (Zpy_data*jetpy_data))/(pow((((Zpx_data*Zpx_data) + (Zpy_data*Zpy_data))*((jetpx_data*jetpx_data) + (jetpy_data*jetpy_data))),0.5)));
		
		double mupeta_data = -log(tan(acos(muppz_data/sqrt(muppx_data*muppx_data+muppy_data*muppy_data+muppz_data*muppz_data))/2));
		
		double mumeta_data = -log(tan(acos(mumpz_data/sqrt(mumpx_data*mumpx_data+mumpy_data*mumpy_data+mumpz_data*mumpz_data))/2));
		
		if (totCand_data != 1) continue; // select events with only one jet
		
		if (jetpt_data < jetpt_min || jetpt_data > jetpt_max) continue; // cut in pT of jet
		
		if (jeteta_data < 2.2 || jeteta_data > 4.2) continue; // cut in eta of jet
		
		if (phi_data < 2.8) continue; // cut in DeltaPhi
		
		if (muppt_data < 2e4 || mumpt_data < 2e4) continue; // cut in pT of muons
		
		if (mupeta_data < 2 || mupeta_data > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_data < 2 || mumeta_data > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_data == 1 && mupHLT1MuonHighPt_data == 1 && mupHLT2_data == 1) || (mumL0MuonEW_data == 1 && mumHLT1MuonHighPt_data == 1 && mumHLT2_data == 1)) { // cuts in trigger lines on muons
			
			Rdata[d] = jetpt_data / Zpt_data;
			
			hRdata->Fill(Rdata[d]);
			
		}
		
	}
	
	cout << endl;
	
	// array to store K values
	double K[10] = {0.0};
	K[0] = 0.97;
	
	for (int i=0; i<10; i++) {
		
		K[i+1] = K[i] + 0.01;
		
		cout << "K[" << i << "] = " << K[i] << endl;
		
	}
	
	cout << endl;
	
	// arrays to store Rmean(K) and sigmaRmean(K) values
	double RmeanK[10] = {0.0};
	double sigmaRmeanK[10] = {0.0};
	
	for (int k=0; k<10; k++) {
		
		for (int m=0; m<Nmc+1; m++) {
			
			if (m % 1000 == 0) cout << "MC[" << k << "] event " << m << endl;
			
			tmc->GetEntry(m);
			
			double phi_mc = acos(((Zpx_mc*jetpx_mc) + (Zpy_mc*jetpy_mc))/(pow((((Zpx_mc*Zpx_mc) + (Zpy_mc*Zpy_mc))*((jetpx_mc*jetpx_mc) + (jetpy_mc*jetpy_mc))),0.5)));
			
			double mupeta_mc = -log(tan(acos(muppz_mc/sqrt(muppx_mc*muppx_mc+muppy_mc*muppy_mc+muppz_mc*muppz_mc))/2));
			
			double mumeta_mc = -log(tan(acos(mumpz_mc/sqrt(mumpx_mc*mumpx_mc+mumpy_mc*mumpy_mc+mumpz_mc*mumpz_mc))/2));
			
			if (totCand_mc != 1) continue; // select events with only one jet
			
			if (jetpt_mc < jetpt_min || jetpt_mc > jetpt_max) continue; // cut in pT of jet
	
			if (jeteta_mc < 2.2 || jeteta_mc > 4.2) continue; // cut in eta of jet
			
			if (phi_mc < 2.8) continue; // cut in DeltaPhi
			
			if (muppt_mc < 2e4 || mumpt_mc < 2e4) continue; // cut in pT of muons
			
			if (mupeta_mc < 2 || mupeta_mc > 4.5) continue; // cut in eta of mu+
			
			if (mumeta_mc < 2 || mumeta_mc > 4.5) continue; // cut in eta of mu-
			
			if ((mupL0MuonEW_mc == 1 && mupHLT1MuonHighPt_mc == 1 && mupHLT2_mc == 1) || (mumL0MuonEW_mc == 1 && mumHLT1MuonHighPt_mc == 1 && mumHLT2_mc == 1)) { // cuts in trigger lines on muons
				
				Rmc[m] = (K[k] * jetpt_mc) / Zpt_mc;
				
				hRmc->Fill(Rmc[m]);
				
			}
			
		}
		
		double Nsel_mc = hRmc->Integral();
		
		cout << "\nNsel_mc = " << Nsel_mc << "\n" << endl;
		
		hRmc->Scale(1/(hRmc->Integral()));
		
		RmeanK[k] = hRmc->GetMean();
		sigmaRmeanK[k] = hRmc->GetMeanError();
		
	}
	
	double Nsel_data = hRdata->Integral();
	
	cout << "Nsel_data = " << Nsel_data << endl;
	
	hRdata->Scale(1/(hRdata->Integral()));
			
	double Rmean_data = hRdata->GetMean();
	double sigmaRmean_data = hRdata->GetMeanError();
	
	cout << "\n<R>data      = " << Rmean_data << endl;
	cout << "sigma<R>data = " << sigmaRmean_data << endl;
	
	// ********************************** //
	// * INSERT NUMBER OF BIN MANUALLY! * //
	// ********************************** //
	ofstream outfile ("RvsKbin11.txt");
	
	for (int f=0; f<10; f++) {
		
		outfile << K[f] << "\t" << RmeanK[f] << "\t" << sigmaRmeanK[f] << endl;
		
		cout << "\n<R>(K[" << f << "])      = " << RmeanK[f] << endl;
		cout << "sigma<R>(K[" << f << "]) = " << sigmaRmeanK[f] << endl;
		
	}
	
	cout << endl;
	
	TCanvas *c = new TCanvas("c","RvsK");
	
	TF1 *RmeanData = new TF1("RmeanData","[0]",0.94,1.1);
	TF1 *RmeanUp = new TF1("RmeanUp","[1]+[2]",0.94,1.1);
	TF1 *RmeanDown = new TF1("RmeanDown","[3]-[4]",0.94,1.1);
	RmeanData->SetParameter(0,Rmean_data);
	RmeanUp->SetParameter(1,Rmean_data);
	RmeanUp->SetParameter(2,sigmaRmean_data);
	RmeanDown->SetParameter(3,Rmean_data);
	RmeanDown->SetParameter(4,sigmaRmean_data);
	
	// ********************************** //
	// * INSERT NUMBER OF BIN MANUALLY! * //
	// ********************************** //
	TGraphErrors *R_mc = new TGraphErrors("RvsKbin11.txt","%lg %lg"); // NO ERRORS!!!
	
	TF1 *fit = new TF1("Linear fit of Rmc","[0]*x+[1]",0.94,1.1);
	
	fit->SetParNames("m","q");
	
	R_mc->Fit(fit);
	
	double m = fit->GetParameter(0);
	double q = fit->GetParameter(1);
	double sigma_m = fit->GetParError(0);
	double sigma_q = fit->GetParError(1);
	
	TF1 *fitplot = new TF1("Linear fit of Rmc","[0]*x+[1]",0.94,1.1);
	fitplot->SetParameter(0,m);
	fitplot->SetParameter(1,q);
	
	RmeanData->SetTitle("<R_{Z}> vs K - 60 GeV < p_{T}(jet) < 65 GeV;K;<R_{Z}>");
	RmeanData->GetXaxis()->SetRangeUser(0.94,1.1);
	RmeanData->GetYaxis()->SetRangeUser(1,1.18);
	RmeanUp->SetLineStyle(2);
	RmeanUp->SetLineColor(kRed);
	RmeanDown->SetLineStyle(2);
	RmeanDown->SetLineColor(kRed);
	
	R_mc->SetMarkerColor(kBlue+2);
	R_mc->SetMarkerStyle(20);
	R_mc->SetMarkerSize(1);
	
	fitplot->SetLineColor(kBlue+2);
	
	double Kl = (Rmean_data - sigmaRmean_data - q) / m;
	double Kh = (Rmean_data + sigmaRmean_data - q) / m;
	double Kstar = (Rmean_data - q) / m; // JES
	
	double DeltaJES = Kh - Kl;
	
	double ErrJES = (1 / m) * sqrt( (sigmaRmean_data*sigmaRmean_data) + (sigma_q*sigma_q) + ( ((Rmean_data-q)/m)*((Rmean_data-q)/m)*sigma_m*sigma_m ) );
	
	double ErrDeltaJES = (1 / m) * sqrt( (2*sigmaRmean_data*sigmaRmean_data) + (2*sigma_q*sigma_q) + ( ((Kh*Kh)+(Kl*Kl))*((sigma_m/m)*(sigma_m/m))) );
	
	double ErrRelErrJES = (1 / Kstar) * sqrt( (ErrDeltaJES*ErrDeltaJES) + ( ((DeltaJES/Kstar)*(DeltaJES/Kstar))*(ErrJES*ErrJES) ) );
	
	gStyle->SetOptFit(111);
	
	c->cd();
	
	RmeanData->Draw();
	RmeanUp->Draw("SAME");
	RmeanDown->Draw("SAME");
	fitplot->Draw("SAME");
	R_mc->Draw("P SAME");
	
	TLegend leg(.1,.65,.4,.9); // (x1,y1,x2,y2)
	
	leg.AddEntry(RmeanData,"<R_{Z}>_{data}","L");
	leg.AddEntry(R_mc,"<R_{Z}>_{MC}(K)","P");
	leg.AddEntry(fitplot,"Fit line of <R_{Z}>_{MC}(K)","L");
	leg.DrawClone("SAME");
	
	cout << "\nKl       = " << Kl << "\nKh       = " << Kh << "\nKstar    = " << Kstar << "\n\nDeltaJES = " << DeltaJES << "\n" << endl;
	
	// add row to files with jetpT (central value of bin, computed at the beginning of this function) and corresponding DeltaJES, JES and DeltaJES/JES
	
//	ofstream outfileDJES("DeltaJESvsJetpT.txt",ios::app);
//	outfileDJES << jetpt_centr/1000 << "\t" << DeltaJES << "\t" << ErrDeltaJES << endl;
	
//	ofstream outfileJES("JESvsJetpT.txt",ios::app);
//	outfileJES << jetpt_centr/1000 << "\t" << Kstar << "\t" << ErrJES << endl;
	
//	ofstream outfileRelErrJES("RelErrJESvsJetpT.txt",ios::app);
//	outfileRelErrJES << jetpt_centr/1000 << "\t" << DeltaJES/Kstar << "\t" << ErrRelErrJES << endl;
	
	// for each bin add Rmean_data, sigmaRmean_data, m, sigma_m, q, sigma_q, Kl, Kh, Kstar=JES, ErrJES
	
	// ********************************** //
	// * INSERT NUMBER OF BIN MANUALLY! * //
	// ********************************** //
	ofstream outfilebin("RvsKbin11.txt",ios::app);
	outfilebin << Rmean_data << "\n" << sigmaRmean_data << "\n" << m << "\n" << sigma_m << "\n" << q << "\n" << sigma_q << "\n" << Kl << "\n" << Kh << "\n" << Kstar << "\n" << ErrJES << endl;
	
}








void DeltaJESvsJetpT() {
	
	TCanvas *c1 = new TCanvas("c1","DeltaJES vs jetpt");
	
	TGraphErrors *gDJES = new TGraphErrors("DeltaJESvsJetpT.txt","%lg %lg %lg");

	gDJES->SetTitle("#DeltaK* vs p_{T}(jet);p_{T}(jet) [GeV];#DeltaK*");
	gDJES->SetMarkerColor(kViolet+6);
	gDJES->SetMarkerStyle(20);
	gDJES->SetMarkerSize(1);
	gDJES->GetXaxis()->SetRangeUser(20,100);

	gDJES->Draw("AP");
	
}








void JESvsJetpT() {
	
	TCanvas *c2 = new TCanvas("c2","JES vs jetpt");
	
	TGraphErrors *gJES = new TGraphErrors("JESvsJetpT.txt","%lg %lg %lg");

	gJES->SetTitle("K* vs p_{T}(jet);p_{T}(jet) [GeV];K*");
	gJES->SetMarkerColor(kViolet+6);
	gJES->SetMarkerStyle(20);
	gJES->SetMarkerSize(1.2);
	gJES->GetXaxis()->SetRangeUser(20,100);

	gJES->Draw("AP");
	
}








void RelErrJESvsJetpT() {
	
	TCanvas *c3 = new TCanvas("c3","RelErrJES vs jetpt");
	
	TGraphErrors *gRelErr = new TGraphErrors("RelErrJESvsJetpT.txt","%lg %lg %lg");
	
	gRelErr->SetTitle("#DeltaK*/K* vs p_{T}(jet);p_{T}(jet) [GeV];#DeltaK*/K*");
	gRelErr->SetMarkerColor(kViolet+6);
	gRelErr->SetMarkerStyle(20);
	gRelErr->SetMarkerSize(1);
	gRelErr->GetXaxis()->SetRangeUser(20,100);

	gRelErr->Draw("AP");
	
}








void ErrJES() {
	
	ifstream infile("RvsKbin18.txt");
	ofstream outfile("RvsKbin18.txt",ios::app);
	
	for (int i=0; i<10; i++) { // skip the first 10 lines
		
		string line;
		getline(infile,line);
		
	}
	
	double R = 0.0;
	double sR = 0.0;
	double m = 0.0;
	double sm = 0.0;
	double q = 0.0;
	double sq = 0.0;
	double ErrJES = 0.0;
	
	while (infile >> R >> sR >> m >> sm >> q >> sq) {
		
		cout << "\n<R> = " << R << "\n" << "sigma_<R> = " << sR << "\n\n" << "m = " << m << "\n" << "sigma_m = " << sm << "\n\n" << "q = " << q << "\n" << "sigma_q = " << sq << endl;
		
		ErrJES = (1 / m) * sqrt( (sR*sR) + (sq*sq) + ( ((R-q)/m)*((R-q)/m)*sm*sm ) );
		
		cout << "\nErrJES = " << ErrJES << "\n" << endl;
		
		outfile << ErrJES << endl;
		
	}
	
}
