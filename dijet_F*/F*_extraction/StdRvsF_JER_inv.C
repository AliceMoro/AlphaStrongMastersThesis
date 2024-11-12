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

void getStdRvsFbins() { // MC fixed, data smeared
	
	// ******************************* //
	// * INSERT BIN LIMITS MANUALLY! * //
	// ******************************* //
	double jetpt_min = 6e4; // MeV - jetpt_data >= 20 GeV
	double jetpt_max = 7e4; // MeV
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
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
		
		double Dphi_mc = acos(((jet1px_mc*jet2px_mc) + (jet1py_mc*jet2py_mc))/(pow((((jet1px_mc*jet1px_mc) + (jet1py_mc*jet1py_mc))*((jet2px_mc*jet2px_mc) + (jet2py_mc*jet2py_mc))),0.5)));
		
		if (jet1pt_mc < jetpt_min || jet1pt_mc > jetpt_max || jet2pt_mc < jetpt_min || jet2pt_mc > jetpt_max) continue; // cut in pT of jets
		
		if (Dphi_mc < 2.8) continue; // cut in DeltaPhi
		
		double rm = RanGen->Uniform();
		
		if (rm <= 0.5) {
			Rmc[m] = (jet1pt_mc-jet2pt_mc) / (jet1pt_mc+jet2pt_mc);
		}
		
		else if (rm > 0.5) {
			Rmc[m] = (jet2pt_mc-jet1pt_mc) / (jet2pt_mc+jet1pt_mc);
		}
		
		hRmc->Fill(Rmc[m],w);
		
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
		
		for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data[" << k << "] event " << d << endl;
		
		tdata->GetEntry(d);
		
		double Dphi_data = acos(((jet1px_data*jet2px_data) + (jet1py_data*jet2py_data))/(pow((((jet1px_data*jet1px_data) + (jet1py_data*jet1py_data))*((jet2px_data*jet2px_data) + (jet2py_data*jet2py_data))),0.5)));
		
		if (jet1pt_data < jetpt_min || jet1pt_data > jetpt_max || jet2pt_data < jetpt_min || jet2pt_data > jetpt_max) continue; // cut in pT of jets
		
		if (Dphi_data < 2.8) continue; // cut in DeltaPhi
		
		double rd = RanGen->Uniform();
		
		if (rd <= 0.5) {
			Rdata[d] = (RanGen->Gaus(jet1pt_data,F[k]*jet1pt_data)-RanGen->Gaus(jet2pt_data,F[k]*jet2pt_data)) / (RanGen->Gaus(jet1pt_data,F[k]*jet1pt_data)+RanGen->Gaus(jet2pt_data,F[k]*jet2pt_data));
		}
		
		if (rd > 0.5) {
			Rdata[d] = (RanGen->Gaus(jet2pt_data,F[k]*jet2pt_data)-RanGen->Gaus(jet1pt_data,F[k]*jet1pt_data)) / (RanGen->Gaus(jet2pt_data,F[k]*jet2pt_data)+RanGen->Gaus(jet1pt_data,F[k]*jet1pt_data));
		}
		
		hRdata->Fill(Rdata[d]);
		
		}
		
		double Nsel_data = hRdata->Integral();
		
		cout << "\nNsel_data = " << Nsel_data << "\n" << endl;
		
		hRdata->Scale(1/(hRdata->Integral()));
		
		stdR_F[k] = hRdata->GetStdDev();
		sigma_stdR_F[k] = hRdata->GetStdDevError();
		
	}
	
	double Nsel_mc = hRmc->Integral();
	
	cout << "Nsel_mc = " << Nsel_mc << endl;
	
	hRmc->Scale(1/(hRmc->Integral()));
			
	double stdR_mc = hRmc->GetStdDev();
	double sigma_stdR_mc = hRmc->GetStdDevError();
	
	cout << "\nstdRmc      = " << stdR_mc << endl;
	cout << "sigma_stdRmc = " << sigma_stdR_mc << endl;
	
	// ********************************** //
	// * INSERT NUMBER OF BIN MANUALLY! * //
	// ********************************** //
	ofstream outfile ("stdRvsFbin05.txt");
	
	for (int f=0; f<100; f++) {
		
		outfile << F[f] << "\t" << stdR_F[f] << "\t" << sigma_stdR_F[f] << endl;
		
		cout << "\nstdR(F[" << f << "])      = " << stdR_F[f] << endl;
		cout << "sigma_stdR(F[" << f << "]) = " << sigma_stdR_F[f] << endl;
		
	}
	
	cout << endl;
	
	TCanvas *c = new TCanvas("c","stdRvsF");
	
	TF1 *stdRMC = new TF1("stdRMC","[0]",0,0.12);
	TF1 *stdRUp = new TF1("stdRUp","[1]+[2]",0,0.12);
	TF1 *stdRDown = new TF1("stdRDown","[3]-[4]",0,0.12);
	stdRMC->SetParameter(0,stdR_mc);
	stdRUp->SetParameter(1,stdR_mc);
	stdRUp->SetParameter(2,sigma_stdR_mc);
	stdRDown->SetParameter(3,stdR_mc);
	stdRDown->SetParameter(4,sigma_stdR_mc);
	
	// ********************************** //
	// * INSERT NUMBER OF BIN MANUALLY! * //
	// ********************************** //
	TGraphErrors *stdR_data = new TGraphErrors("stdRvsFbin05.txt","%lg %lg");
	
	ifstream IF("stdRvsFbin05.txt");
	ofstream OF("stdRvsFbin05FIT.txt");
	
	string line;
	int curr_line = 1;
	
	while (getline(IF,line)) {
		
		if (curr_line >= 6 && curr_line <= 26) {
			
			OF << line << endl;
			
		}
		
		if (curr_line > 26) break;
		
		curr_line++;
		
	}
	
	TGraphErrors *stdR_dataFIT = new TGraphErrors("stdRvsFbin05FIT.txt","%lg %lg");
	
	TF1 *fit = new TF1("Linear fit of stdRdata","[0]*x+[1]",0,0.12);
	
	fit->SetParNames("m","q");
	
	stdR_dataFIT->Fit(fit);
	
	double m = fit->GetParameter(0);
	double q = fit->GetParameter(1);
	double sigma_m = fit->GetParError(0);
	double sigma_q = fit->GetParError(1);
	
	TF1 *fitplot = new TF1("Linear fit of stdRdata","[0]*x+[1]",0,0.12);
	fitplot->SetParameter(0,m);
	fitplot->SetParameter(1,q);
	
	stdRMC->SetTitle("<#sigma_{R_{12}}> vs F - 60 GeV < p_{T}(jet_{1,2}) < 70 GeV;F;<#sigma_{R_{12}}>");
	stdRMC->GetXaxis()->SetRangeUser(0,0.06);
	stdRMC->GetYaxis()->SetRangeUser(0.02,0.08);
	stdRMC->SetLineColor(kBlue+2);
	stdRUp->SetLineStyle(2);
	stdRUp->SetLineColor(kBlue+2);
	stdRDown->SetLineStyle(2);
	stdRDown->SetLineColor(kBlue+2);
	
	stdR_data->SetMarkerColor(kRed);
	stdR_data->SetMarkerStyle(20);
	stdR_data->SetMarkerSize(1);
	
	stdR_dataFIT->SetMarkerColor(kRed);
	stdR_dataFIT->SetMarkerStyle(20);
	stdR_dataFIT->SetMarkerSize(1);
	
	fitplot->SetLineColor(kRed);
	
	double Fl = (stdR_mc - sigma_stdR_mc - q) / m;
	double Fh = (stdR_mc + sigma_stdR_mc - q) / m;
	double Fstar = (stdR_mc - q) / m;
	
	double DeltaJER = Fh - Fl;
	
	double ErrJER = (1 / m) * sqrt( (sigma_stdR_mc*sigma_stdR_mc) + (sigma_q*sigma_q) + ( ((stdR_mc-q)/m)*((stdR_mc-q)/m)*sigma_m*sigma_m ) );
	
	double ErrDeltaJER = (1 / m) * sqrt( (2*sigma_stdR_mc*sigma_stdR_mc) + (2*sigma_q*sigma_q) + ( ((Fh*Fh)+(Fl*Fl))*((sigma_m/m)*(sigma_m/m))) );
	
	double ErrRelErrJER = (1 / Fstar) * sqrt( (ErrDeltaJER*ErrDeltaJER) + ( ((DeltaJER/Fstar)*(DeltaJER/Fstar))*(ErrJER*ErrJER) ) );
	
	gStyle->SetOptFit(111);
	
	c->cd();
	
	stdRMC->Draw();
	stdRUp->Draw("SAME");
	stdRDown->Draw("SAME");
	fitplot->Draw("SAME");
	stdR_data->Draw("P SAME");
	stdR_dataFIT->Draw("P SAME");
	
	TLegend leg(.1,.65,.4,.9); // (x1,y1,x2,y2)
	
	leg.AddEntry(stdRMC,"<#sigma_{R_{12}}>_{MC}","L");
	leg.AddEntry(stdR_data,"<#sigma_{R_{12}}>_{data}(F)","P");
	leg.AddEntry(fitplot,"Fit line of <#sigma_{R_{12}}>_{data}(F)","L");
	leg.DrawClone("SAME");
	
	cout << "\nFl       = " << Fl << "\nFh       = " << Fh << "\nFstar    = " << Fstar << "\n\nDeltaJER = " << DeltaJER << "\n" << endl;
	
	// add row to files with jetpT (central value of bin, computed at the beginning of this function) and corresponding DeltaJER, JER and DeltaJER/JER
	
	ofstream outfileDJER("DeltaJERvsJetpT.txt",ios::app);
	outfileDJER << jetpt_centr/1000 << "\t" << DeltaJER << "\t" << ErrDeltaJER << endl;
	
	ofstream outfileJER("JERvsJetpT.txt",ios::app);
	outfileJER << jetpt_centr/1000 << "\t" << Fstar << "\t" << ErrJER << endl;
	
	ofstream outfileRelErrJER("RelErrJERvsJetpT.txt",ios::app);
	outfileRelErrJER << jetpt_centr/1000 << "\t" << DeltaJER/Fstar << "\t" << ErrRelErrJER << endl;
	
	// for each bin add stdR_mc, sigma_stdR_mc, m, sigma_m, q, sigma_q, Fl, Fh, Fstar=JER, ErrJER
	
	// ********************************** //
	// * INSERT NUMBER OF BIN MANUALLY! * //
	// ********************************** //
	ofstream outfilebin("stdRvsFbin05.txt",ios::app);
	outfilebin << stdR_mc << "\n" << sigma_stdR_mc << "\n" << m << "\n" << sigma_m << "\n" << q << "\n" << sigma_q << "\n" << Fl << "\n" << Fh << "\n" << Fstar << "\n" << ErrJER << endl;
	
}
