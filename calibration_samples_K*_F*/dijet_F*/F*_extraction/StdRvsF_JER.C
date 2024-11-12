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
#include <TRandom3.h>
#include <TUUID.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

void getStdRvsFbins() { // MC smeared, data fixed
	
	// ******************************* //
	// * INSERT BIN LIMITS MANUALLY! * //
	// ******************************* //
	double jetpt_min = 2e4; // MeV - jetpt_data >= 20 GeV
	double jetpt_max = 10e4; // MeV
	double jetpt_centr = (jetpt_min + jetpt_max) / 2;
	
	TFile *fdata = new TFile("dijet_data_2016.root");
	TFile *fmc = new TFile("qq.root");
	
	TTree *tdata = (TTree*)fdata->Get("Tuple");
	TTree *tmc = (TTree*)fmc->Get("Tuple");
	
	int Ndata = tdata->GetEntries();
	int Nmc = tmc->GetEntries();
	
	cout << "\nNumber of entries of data tuple = " << Ndata << "\nNumber of entries of MC tuple   = " << Nmc << "\n" << endl;
	
	double jet1pt_data = 0.0;
	double jet2pt_data = 0.0;
	
	double jet1pt_mc = 0.0;
	double jet2pt_mc = 0.0;
	double w = 0.0;
	
	tdata->SetBranchAddress("jet0_pt",&jet1pt_data);
	tdata->SetBranchAddress("jet1_pt",&jet2pt_data);
	
	tmc->SetBranchAddress("jet0_pt",&jet1pt_mc);
	tmc->SetBranchAddress("jet1_pt",&jet2pt_mc);
	tmc->SetBranchAddress("weight",&w);
	
	// to compute DeltaPhi
	double jet1px_data = 0.0;
	double jet1py_data = 0.0;
	double jet2px_data = 0.0;
	double jet2py_data = 0.0;
	
	double jet1px_mc = 0.0;
	double jet1py_mc = 0.0;
	double jet2px_mc = 0.0;
	double jet2py_mc = 0.0;
	
	tdata->SetBranchAddress("jet0_px",&jet1px_data);
	tdata->SetBranchAddress("jet0_py",&jet1py_data);
	tdata->SetBranchAddress("jet1_px",&jet2px_data);
	tdata->SetBranchAddress("jet1_py",&jet2py_data);
	tmc->SetBranchAddress("jet0_px",&jet1px_mc);
	tmc->SetBranchAddress("jet0_py",&jet1py_mc);
	tmc->SetBranchAddress("jet1_px",&jet2px_mc);
	tmc->SetBranchAddress("jet1_py",&jet2py_mc);
	
	// arrays to store R values
	double Rdata[24970] = {0.0};
	double Rmc[513657] = {0.0};
	
	TH1D *hRdata = new TH1D("hRdata","Calibration sample dijet - R = (p_{T}(jet_{1})-p_{T}(jet_{2}))/(p_{T}(jet_{1})+p_{T}(jet_{2})) - data;R;A.U.",1000,-1,1);
	TH1D *hRmc = new TH1D("hRmc","Calibration sample dijet - R = (p_{T}(jet_{1})-p_{T}(jet_{2}))/(p_{T}(jet_{1})+p_{T}(jet_{2})) - MC;R;A.U.",1000,-1,1);
	
	// random generator for mixing order of jet1,jet2 and smearing
	TRandom3 *RanGen = new TRandom3(0);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double Dphi_data = acos(((jet1px_data*jet2px_data) + (jet1py_data*jet2py_data))/(pow((((jet1px_data*jet1px_data) + (jet1py_data*jet1py_data))*((jet2px_data*jet2px_data) + (jet2py_data*jet2py_data))),0.5)));
		
		if (jet1pt_data < jetpt_min || jet1pt_data > jetpt_max || jet2pt_data < jetpt_min || jet2pt_data > jetpt_max) continue; // cut in pT of jets
		
		if (Dphi_data < 2.8) continue; // cut in DeltaPhi
		
		double rd = RanGen->Uniform();
		
		if (rd <= 0.5) {
			Rdata[d] = (jet1pt_data-jet2pt_data) / (jet1pt_data+jet2pt_data);
		}
		
		else if (rd > 0.5) {
			Rdata[d] = (jet2pt_data-jet1pt_data) / (jet2pt_data+jet1pt_data);
		}
		
		hRdata->Fill(Rdata[d]);
		
	}
	
	cout << endl;
	
	// array to store F values
	double F[100] = {0.0};
	F[0] = 0.005;
	
	for (int i=0; i<100; i++) {
		
		F[i+1] = F[i] + 0.001;
		
		cout << "F[" << i << "] = " << F[i] << endl;
		
	}
	
	cout << endl;
	
	// arrays to store stdR(F) and sigma_stdR(F) values
	double stdR_F[100] = {0.0};
	double sigma_stdR_F[100] = {0.0};
	
	for (int k=0; k<100; k++) {
		
		for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC[" << k << "] event " << m << endl;
		
		tmc->GetEntry(m);
		
		double Dphi_mc = acos(((jet1px_mc*jet2px_mc) + (jet1py_mc*jet2py_mc))/(pow((((jet1px_mc*jet1px_mc) + (jet1py_mc*jet1py_mc))*((jet2px_mc*jet2px_mc) + (jet2py_mc*jet2py_mc))),0.5)));
		
		if (jet1pt_mc < jetpt_min || jet1pt_mc > jetpt_max || jet2pt_mc < jetpt_min || jet2pt_mc > jetpt_max) continue; // cut in pT of jets
		
		if (Dphi_mc < 2.8) continue; // cut in DeltaPhi
		
		double rm = RanGen->Uniform();
		
		if (rm <= 0.5) {
			Rmc[m] = (RanGen->Gaus(jet1pt_mc,F[k]*jet1pt_mc)-RanGen->Gaus(jet2pt_mc,F[k]*jet2pt_mc)) / (RanGen->Gaus(jet1pt_mc,F[k]*jet1pt_mc)+RanGen->Gaus(jet2pt_mc,F[k]*jet2pt_mc));
		}
		
		if (rm > 0.5) {
			Rmc[m] = (RanGen->Gaus(jet2pt_mc,F[k]*jet2pt_mc)-RanGen->Gaus(jet1pt_mc,F[k]*jet1pt_mc)) / (RanGen->Gaus(jet2pt_mc,F[k]*jet2pt_mc)+RanGen->Gaus(jet1pt_mc,F[k]*jet1pt_mc));
		}
		
		hRmc->Fill(Rmc[m],w);
		
		}
		
		double Nsel_mc = hRmc->Integral();
		
		cout << "\nNsel_mc = " << Nsel_mc << "\n" << endl;
		
		hRmc->Scale(1/(hRmc->Integral()));
		
		stdR_F[k] = hRmc->GetStdDev();
		sigma_stdR_F[k] = hRmc->GetStdDevError();
		
	}
	
	double Nsel_data = hRdata->Integral();
	
	cout << "Nsel_data = " << Nsel_data << endl;
	
	hRdata->Scale(1/(hRdata->Integral()));
			
	double stdR_data = hRdata->GetStdDev();
	double sigma_stdR_data = hRdata->GetStdDevError();
	
	cout << "\nstdRdata      = " << stdR_data << endl;
	cout << "sigma_stdRdata = " << sigma_stdR_data << endl;
	
	// ********************************** //
	// * INSERT NUMBER OF BIN MANUALLY! * //
	// ********************************** //
	ofstream outfile ("stdRvsFbinTOT.txt");
	
	for (int f=0; f<100; f++) {
		
		outfile << F[f] << "\t" << stdR_F[f] << "\t" << sigma_stdR_F[f] << endl;
		
		cout << "\nstdR(F[" << f << "])      = " << stdR_F[f] << endl;
		cout << "sigma_stdR(F[" << f << "]) = " << sigma_stdR_F[f] << endl;
		
	}
	
	cout << endl;
	
	TCanvas *c = new TCanvas("c","stdRvsF");
	
	TF1 *stdRData = new TF1("stdRData","[0]",0,2);
	TF1 *stdRUp = new TF1("stdRUp","[1]+[2]",0,2);
	TF1 *stdRDown = new TF1("stdRDown","[3]-[4]",0,2);
	stdRData->SetParameter(0,stdR_data);
	stdRUp->SetParameter(1,stdR_data);
	stdRUp->SetParameter(2,sigma_stdR_data);
	stdRDown->SetParameter(3,stdR_data);
	stdRDown->SetParameter(4,sigma_stdR_data);
	
	// ********************************** //
	// * INSERT NUMBER OF BIN MANUALLY! * //
	// ********************************** //
	TGraphErrors *stdR_mc = new TGraphErrors("stdRvsFbinTOT.txt","%lg %lg");
	
	TF1 *fit = new TF1("Linear fit of stdRmc","[0]*x+[1]",0,2);
	
	fit->SetParNames("m","q");
	
	stdR_mc->Fit(fit);
	
	double m = fit->GetParameter(0);
	double q = fit->GetParameter(1);
	double sigma_m = fit->GetParError(0);
	double sigma_q = fit->GetParError(1);
	
	stdRData->SetTitle("<#sigma_{R}> vs F;F;<#sigma_{R}>");
	stdRData->GetYaxis()->SetRangeUser(0,2);
	stdRUp->SetLineWidth(1);
	stdRUp->SetLineColor(kBlack);
	stdRDown->SetLineWidth(1);
	stdRDown->SetLineColor(kBlack);
	
	stdR_mc->SetMarkerColor(kBlue+2);
	stdR_mc->SetMarkerStyle(20);
	stdR_mc->SetMarkerSize(0.8);
	
	double Fl = (stdR_data - sigma_stdR_data - q) / m;
	double Fh = (stdR_data + sigma_stdR_data - q) / m;
	double Fstar = (stdR_data - q) / m;
	
	double DeltaJER = Fh - Fl;
	
	double ErrJER = (1 / m) * sqrt( (sigma_stdR_data*sigma_stdR_data) + (sigma_q*sigma_q) + ( ((stdR_data-q)/m)*((stdR_data-q)/m)*sigma_m*sigma_m ) );
	
	double ErrDeltaJER = (1 / m) * sqrt( (2*sigma_stdR_data*sigma_stdR_data) + (2*sigma_q*sigma_q) + ( ((Fh*Fh)+(Fl*Fl))*((sigma_m/m)*(sigma_m/m))) );
	
	double ErrRelErrJER = (1 / Fstar) * sqrt( (ErrDeltaJER*ErrDeltaJER) + ( ((DeltaJER/Fstar)*(DeltaJER/Fstar))*(ErrJER*ErrJER) ) );
	
	stdRData->GetXaxis()->SetRangeUser(0,0.11);
	
	gStyle->SetOptFit(111);
	
	c->cd();
	
	stdRData->Draw();
	stdRUp->Draw("SAME");
	stdRDown->Draw("SAME");
	stdR_mc->Draw("P SAME");
	
	TLegend leg(.1,.75,.35,.9); // (x1,y1,x2,y2)
	
	leg.AddEntry("stdRData","<#sigma_{R}>_{data}");
	leg.AddEntry("stdR_mc","<#sigma_{R}>_{MC}(F)","PL");
	leg.DrawClone("SAME");
	
	cout << "\nFl       = " << Fl << "\nFh       = " << Fh << "\nFstar    = " << Fstar << "\n\nDeltaJER = " << DeltaJER << "\n" << endl;
	
	// add row to files with jetpT (central value of bin, computed at the beginning of this function) and corresponding DeltaJER, JER and DeltaJER/JER
	
//	ofstream outfileDJER("DeltaJERvsJetpT.txt",ios::app);
//	outfileDJER << jetpt_centr/1000 << "\t" << DeltaJER << "\t" << ErrDeltaJER << endl;
	
//	ofstream outfileJER("JERvsJetpT.txt",ios::app);
//	outfileJER << jetpt_centr/1000 << "\t" << Fstar << "\t" << ErrJER << endl;
	
//	ofstream outfileRelErrJER("RelErrJERvsJetpT.txt",ios::app);
//	outfileRelErrJER << jetpt_centr/1000 << "\t" << DeltaJER/Fstar << "\t" << ErrRelErrJER << endl;
	
	// for each bin add stdR_data, sigma_stdR_data, m, sigma_m, q, sigma_q, Fl, Fh, Fstar=JER, ErrJER, DeltaJER, ErrDeltaJER, RelErrJER, ErrRelErrJER
	
	// ********************************** //
	// * INSERT NUMBER OF BIN MANUALLY! * //
	// ********************************** //
	ofstream outfilebin("stdRvsFbinTOT.txt",ios::app);
	outfilebin << stdR_data << "\n" << sigma_stdR_data << "\n" << m << "\n" << sigma_m << "\n" << q << "\n" << sigma_q << "\n" << Fl << "\n" << Fh << "\n" << Fstar << "\n" << ErrJER << "\n" << DeltaJER << "\n" << ErrDeltaJER << "\n" << DeltaJER/Fstar << "\n" << ErrRelErrJER << endl;
	
}








void DeltaJERvsJetpT() {
	
	TCanvas *c1 = new TCanvas("c1","DeltaJER vs jetpt");
	
	TGraphErrors *gDJER = new TGraphErrors("DeltaJERvsJetpT.txt","%lg %lg %lg");

	gDJER->SetTitle("#DeltaF* vs p_{T}(jet_{1,2});p_{T}(jet_{1,2}) [GeV];#DeltaF*");
	gDJER->SetMarkerColor(kViolet+6);
	gDJER->SetMarkerStyle(20);
	gDJER->SetMarkerSize(1);
	gDJER->GetXaxis()->SetRangeUser(20,100);

	gDJER->Draw("AP");
	
}








void JERvsJetpT() {
	
	TCanvas *c2 = new TCanvas("c2","JER vs jetpt");
	
	TGraphErrors *gJER = new TGraphErrors("JERvsJetpT.txt","%lg %lg %lg");

	gJER->SetTitle("F* vs p_{T}(jet_{1,2});p_{T}(jet_{1,2}) [GeV];F*");
	gJER->SetMarkerColor(kViolet+6);
	gJER->SetMarkerStyle(20);
	gJER->SetMarkerSize(1.2);
	gJER->GetXaxis()->SetRangeUser(20,100);

	gJER->Draw("AP");
	
}








void RelErrJERvsJetpT() {
	
	TCanvas *c3 = new TCanvas("c3","RelErrJER vs jetpt");
	
	TGraphErrors *gRelErr = new TGraphErrors("RelErrJERvsJetpT.txt","%lg %lg %lg");
	
	gRelErr->SetTitle("#DeltaF*/F* vs p_{T}(jet_{1,2});p_{T}(jet_{1,2}) [GeV];#DeltaF*/F*");
	gRelErr->SetMarkerColor(kViolet+6);
	gRelErr->SetMarkerStyle(20);
	gRelErr->SetMarkerSize(1);
	gRelErr->GetXaxis()->SetRangeUser(20,100);

	gRelErr->Draw("AP");
	
}
