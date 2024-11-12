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
#include <TRandom3.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

void Rdijet() { // R = (pt(jet1)-pt(jet2))/(pt(jet1)+pt(jet2))
	
	double jetpt_min = 2e4; // MeV - jetpt_data >= 20 GeV
	double jetpt_max = 10e4; // MeV
	
	TFile *fdata = new TFile("/home/alice/Desktop/TESI/calibration_samples/dijet/more_F/dijet_data_2016.root");
	TFile *fmc = new TFile("/home/alice/Desktop/TESI/calibration_samples/dijet/more_F/qq.root");
	
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
	
	TH1D *hRdata = new TH1D("hRdata","Calibration sample dijet - R_{12} = (p_{T}(jet_{1})-p_{T}(jet_{2}))/(p_{T}(jet_{1})+p_{T}(jet_{2})) - 20 GeV < p_{T}(jet_{1,2}) < 100 GeV;R_{12};A.U.",40,-1,1);
	TH1D *hRmc = new TH1D("hRmc","Calibration sample dijet - R_{12} = (p_{T}(jet_{1})-p_{T}(jet_{2}))/(p_{T}(jet_{1})+p_{T}(jet_{2})) - 20 GeV < p_{T}(jet_{1,2}) < 100 GeV;R_{12};A.U.",40,-1,1);
	TH2D *hpt12data = new TH2D("hpt12data","Calibration sample dijet - p_{T}(jet_{1}) vs p_{T}(jet_{2}) - data;p_{T}(jet_{1}) [GeV];p_{T}(jet_{2}) [GeV]",120,0,120,120,0,120);
	TH2D *hpt12mc = new TH2D("hpt12mc","Calibration sample dijet - p_{T}(jet_{1}) vs p_{T}(jet_{2}) - MC;p_{T}(jet_{1}) [GeV];p_{T}(jet_{2}) [GeV]",120,0,120,120,0,120);
	
	// random generator for mixing order of jet1,jet2
	TRandom3 *RanGen = new TRandom3(0);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double Dphi_data = acos(((jet1px_data*jet2px_data) + (jet1py_data*jet2py_data))/(pow((((jet1px_data*jet1px_data) + (jet1py_data*jet1py_data))*((jet2px_data*jet2px_data) + (jet2py_data*jet2py_data))),0.5)));
		
		if (jet1pt_data < jetpt_min || jet1pt_data > jetpt_max || jet2pt_data < jetpt_min || jet2pt_data > jetpt_max) continue; // cut in pT of jets
		
		if (Dphi_data < 2.8) continue; // cut in DeltaPhi
		
		hpt12data->Fill(jet1pt_data/1000,jet2pt_data/1000);
		
		double rd = RanGen->Uniform();
		
		if (rd <= 0.5) {
			Rdata[d] = (jet1pt_data-jet2pt_data) / (jet1pt_data+jet2pt_data);
		}
		
		else if (rd > 0.5) {
			Rdata[d] = (jet2pt_data-jet1pt_data) / (jet2pt_data+jet1pt_data);
		}
		
		hRdata->Fill(Rdata[d]);
		
	}
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
		
		double Dphi_mc = acos(((jet1px_mc*jet2px_mc) + (jet1py_mc*jet2py_mc))/(pow((((jet1px_mc*jet1px_mc) + (jet1py_mc*jet1py_mc))*((jet2px_mc*jet2px_mc) + (jet2py_mc*jet2py_mc))),0.5)));
		
		if (jet1pt_mc < jetpt_min || jet1pt_mc > jetpt_max || jet2pt_mc < jetpt_min || jet2pt_mc > jetpt_max) continue; // cut in pT of jets
		
		if (Dphi_mc < 2.8) continue; // cut in DeltaPhi
		
		hpt12mc->Fill(jet1pt_mc/1000,jet2pt_mc/1000,w);
		
		double rm = RanGen->Uniform();
		
		if (rm <= 0.5) {
			Rmc[m] = (jet1pt_mc-jet2pt_mc) / (jet1pt_mc+jet2pt_mc);
		}
		
		if (rm > 0.5) {
			Rmc[m] = (jet2pt_mc-jet1pt_mc) / (jet2pt_mc+jet1pt_mc);
		}
		
		hRmc->Fill(Rmc[m],w);
		
	}
	
	double Nsel_data = hRdata->Integral();
	double Nsel_mc = hRmc->Integral();
		
	cout << "\nNsel data = " << Nsel_data << endl;
	cout << "\nNsel MC = " << Nsel_mc << "\n" << endl;
	
	hRdata->Scale(1/(hRdata->Integral()));
	hRmc->Scale(1/(hRmc->Integral()));
	
	TCanvas *c1 = new TCanvas("c1","Rdijet");
	
	hRdata->GetXaxis()->SetRangeUser(-1,1);
	hRdata->SetLineColor(kRed);
	hRdata->SetMarkerStyle(20);
	hRdata->SetMarkerSize(1.2);
	hRdata->SetMarkerColor(kRed);
	hRmc->GetXaxis()->SetRangeUser(-1,1);
	hRmc->SetLineColor(kBlue+2);
	hRmc->SetMarkerStyle(20);
	hRmc->SetMarkerSize(1.2);
	hRmc->SetMarkerColor(kBlue+2);
	gStyle->SetOptStat("");
	
	c1->cd();
	
	hRmc->Draw();
	hRdata->Draw("SAME");
	
	TLegend leg(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg.AddEntry(hRdata,"Data","PL");
	leg.AddEntry(hRmc,"MC","PL");
	leg.DrawClone("SAME");
	
	TCanvas *c9 = new TCanvas("c9","hpt12data",1000,1000);
	
	hpt12data->SetTitle("Calibration sample dijet - p_{T}(jet_{1}) vs p_{T}(jet_{2}) - data;p_{T}(jet_{1}) [GeV];p_{T}(jet_{2}) [GeV]");
	hpt12data->SetMarkerColor(kRed);
	hpt12data->SetMarkerStyle(1);
	hpt12data->SetMarkerSize(1.2);
	gStyle->SetOptStat("");
	
	c9->cd();
	
	hpt12data->Draw();
	
	TCanvas *c10 = new TCanvas("c10","hpt12mc",1000,1000);
	
	hpt12mc->SetTitle("Calibration sample dijet - p_{T}(jet_{1}) vs p_{T}(jet_{2}) - MC;p_{T}(jet_{1}) [GeV];p_{T}(jet_{2}) [GeV]");
	hpt12mc->SetMarkerColor(kBlue+2);
	hpt12mc->SetMarkerStyle(1);
	hpt12mc->SetMarkerSize(1.2);
	gStyle->SetOptStat("");
	
	c10->cd();
	
	hpt12mc->Draw();
	
}








void Dphidijet() {
	
	double jetpt_min = 2e4; // MeV - jetpt_data >= 20 GeV
	double jetpt_max = 10e4; // MeV
	
	TFile *fdata = new TFile("/home/alice/Desktop/TESI/calibration_samples/dijet/more_F/dijet_data_2016.root");
	TFile *fmc = new TFile("/home/alice/Desktop/TESI/calibration_samples/dijet/more_F/qq.root");
	
	TTree *tdata = (TTree*)fdata->Get("Tuple");
	TTree *tmc = (TTree*)fmc->Get("Tuple");
	
	int Ndata = tdata->GetEntries();
	int Nmc = tmc->GetEntries();
	
	cout << "\nNumber of entries of data tuple = " << Ndata << "\nNumber of entries of MC tuple   = " << Nmc << "\n" << endl;
	
	double jet1pt_data = 0.0;
	double jet2pt_data = 0.0;
	
	double jet1pt_mc = 0.0;
	double jet2pt_mc = 0.0;
	
	tdata->SetBranchAddress("jet0_pt",&jet1pt_data);
	tdata->SetBranchAddress("jet1_pt",&jet2pt_data);
	
	tmc->SetBranchAddress("jet0_pt",&jet1pt_mc);
	tmc->SetBranchAddress("jet1_pt",&jet2pt_mc);
	
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
	
	double Dphidata[24970] = {0.0};
	double Dphimc[513657] = {0.0};
	
	TH1D *hDphidata = new TH1D("hDphidata","Calibration sample dijet - #Delta#phi_{jet_{1,2}} - data;#Delta#phi_{jet_{1,2}};A.U.",35,0,3.5);
	TH1D *hDphimc = new TH1D("hDphimc","Calibration sample dijet - #Delta#phi_{jet_{1,2}} - MC;#Delta#phi_{jet_{1,2}};A.U.",35,0,3.5);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		if (jet1pt_data < jetpt_min || jet1pt_data > jetpt_max || jet2pt_data < jetpt_min || jet2pt_data > jetpt_max) continue; // cut in pT of jets
		
		Dphidata[d] = acos(((jet1px_data*jet2px_data) + (jet1py_data*jet2py_data))/(pow((((jet1px_data*jet1px_data) + (jet1py_data*jet1py_data))*((jet2px_data*jet2px_data) + (jet2py_data*jet2py_data))),0.5)));
		
		hDphidata->Fill(Dphidata[d]);
		
	}
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
		
		if (jet1pt_mc < jetpt_min || jet1pt_mc > jetpt_max || jet2pt_mc < jetpt_min || jet2pt_mc > jetpt_max) continue; // cut in pT of jets
		
		Dphimc[m] = acos(((jet1px_mc*jet2px_mc) + (jet1py_mc*jet2py_mc))/(pow((((jet1px_mc*jet1px_mc) + (jet1py_mc*jet1py_mc))*((jet2px_mc*jet2px_mc) + (jet2py_mc*jet2py_mc))),0.5)));
		
		hDphimc->Fill(Dphimc[m]);
		
	}
	
	double Nsel_data = hDphidata->Integral();
	double Nsel_mc = hDphimc->Integral();
	
	cout << "\nNsel data = " << Nsel_data << endl;
	cout << "\nNsel MC = " << Nsel_mc << "\n" << endl;
	
	hDphidata->Scale(1/(hDphidata->Integral()));
	hDphimc->Scale(1/(hDphimc->Integral()));
	
	TCanvas *c2 = new TCanvas("c2","Dphidijet");
	
	hDphidata->GetXaxis()->SetRangeUser(0,3.4);
	hDphimc->GetXaxis()->SetRangeUser(0,3.4);
	hDphidata->GetYaxis()->SetRangeUser(0,0.25);
	hDphimc->GetYaxis()->SetRangeUser(0,0.25);
	hDphidata->SetTitle("Calibration sample dijet - #Delta#phi_{jet_{1,2}};#Delta#phi_{jet_{1,2}};A.U.");
	hDphimc->SetTitle("Calibration sample dijet - #Delta#phi_{jet_{1,2}};#Delta#phi_{jet_{1,2}};A.U.");
	hDphidata->SetLineColor(kRed);
	hDphidata->SetMarkerStyle(20);
	hDphidata->SetMarkerSize(1);
	hDphidata->SetMarkerColor(kRed);
	hDphimc->SetLineColor(kBlue+2);
	hDphimc->SetMarkerStyle(20);
	hDphimc->SetMarkerSize(1);
	hDphimc->SetMarkerColor(kBlue+2);
	
	TLine *xline = new TLine(2.8,0,2.8,0.25);
	
	xline->SetLineColor(kBlack);
	xline->SetLineWidth(3);
	
	gStyle->SetOptStat("");
	
	c2->cd();
	
	hDphidata->Draw();
	hDphimc->Draw("SAME");
	xline->Draw("SAME");
	
	TLegend leg2(.1,.65,.25,.9); // (x1,y1,x2,y2)
	
	leg2.AddEntry(hDphidata,"Data","PL");
	leg2.AddEntry(hDphimc,"MC","PL");
	leg2.AddEntry(xline,"Cut","L");
	leg2.DrawClone("SAME");
	
}








void GetJetpT() {
	
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
	
	tdata->SetBranchAddress("jet0_pt",&jet1pt_data);
	tdata->SetBranchAddress("jet1_pt",&jet2pt_data);
	
	tmc->SetBranchAddress("jet0_pt",&jet1pt_mc);
	tmc->SetBranchAddress("jet1_pt",&jet2pt_mc);
	
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
	
	TH1D *hpT1data = new TH1D("hpT1data","Calibration sample dijet - p_{T}(jet_{1}) - data;p_{T}(jet_{1}) [GeV];A.U.",60,0,300);
	TH1D *hpT2data = new TH1D("hpT2data","Calibration sample dijet - p_{T}(jet_{2}) - data;p_{T}(jet_{2}) [GeV];A.U.",60,0,300);
	TH1D *hpT1mc = new TH1D("hpT1mc","Calibration sample dijet - p_{T}(jet_{1}) - MC;p_{T}(jet_{1}) [GeV];A.U.",60,0,300);
	TH1D *hpT2mc = new TH1D("hpT2mc","Calibration sample dijet - p_{T}(jet_{2}) - MC;p_{T}(jet_{2}) [GeV];A.U.",60,0,300);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double Dphi_data = acos(((jet1px_data*jet2px_data) + (jet1py_data*jet2py_data))/(pow((((jet1px_data*jet1px_data) + (jet1py_data*jet1py_data))*((jet2px_data*jet2px_data) + (jet2py_data*jet2py_data))),0.5)));
		
		if (Dphi_data < 2.8) continue; // cut in DeltaPhi
		
		hpT1data->Fill(jet1pt_data/1000);
		hpT2data->Fill(jet2pt_data/1000);
		
	}
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
		
		double Dphi_mc = acos(((jet1px_mc*jet2px_mc) + (jet1py_mc*jet2py_mc))/(pow((((jet1px_mc*jet1px_mc) + (jet1py_mc*jet1py_mc))*((jet2px_mc*jet2px_mc) + (jet2py_mc*jet2py_mc))),0.5)));
		
		if (Dphi_mc < 2.8) continue; // cut in DeltaPhi
			
		hpT1mc->Fill(jet1pt_mc/1000);
		hpT2mc->Fill(jet2pt_mc/1000);
		
	}
	
	hpT1data->Scale(1/(hpT1data->Integral()));
	hpT2data->Scale(1/(hpT2data->Integral()));
	hpT1mc->Scale(1/(hpT1mc->Integral()));
	hpT2mc->Scale(1/(hpT2mc->Integral()));
	
	TCanvas *c3 = new TCanvas("c3","pTjet1");
	
	hpT1data->SetLineColor(kRed);
	hpT1data->SetTitle("Calibration sample dijet - p_{T}(jet_{1});p_{T}(jet_{1}) [GeV];A.U.");
	gStyle->SetOptStat("");
	
	c3->cd();
	
	hpT1data->Draw();
	hpT1mc->Draw("SAME");
	
	TLegend leg3(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg3.AddEntry("hpT1data","Data");
	leg3.AddEntry("hpT1mc","MC");
	leg3.DrawClone("SAME");
	
	TCanvas *c4 = new TCanvas("c4","pTjet2");
	
	hpT2data->SetLineColor(kRed);
	hpT2data->SetTitle("Calibration sample dijet - p_{T}(jet_{2});p_{T}(jet_{2}) [GeV];A.U.");
	gStyle->SetOptStat("");
	
	c4->cd();
	
	hpT2data->Draw();
	hpT2mc->Draw("SAME");
	
	TLegend leg4(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg4.AddEntry("hpT2data","Data");
	leg4.AddEntry("hpT2mc","MC");
	leg4.DrawClone("SAME");
	
}








void GetJeteta() {
	
	double jetpt_min = 2e4; // MeV - jetpt_data >= 20 GeV
	double jetpt_max = 10e4; // MeV
	
	TFile *fdata = new TFile("dijet_data_2016.root");
	TFile *fmc = new TFile("qq.root");
	
	TTree *tdata = (TTree*)fdata->Get("Tuple");
	TTree *tmc = (TTree*)fmc->Get("Tuple");
	
	int Ndata = tdata->GetEntries();
	int Nmc = tmc->GetEntries();
	
	cout << "\nNumber of entries of data tuple = " << Ndata << "\nNumber of entries of MC tuple   = " << Nmc << "\n" << endl;
	
	double jet1pt_data = 0.0;
	double jet2pt_data = 0.0;
	double jet1eta_data = 0.0;
	double jet2eta_data = 0.0;
	
	double jet1pt_mc = 0.0;
	double jet2pt_mc = 0.0;
	double jet1eta_mc = 0.0;
	double jet2eta_mc = 0.0;
	
	tdata->SetBranchAddress("jet0_pt",&jet1pt_data);
	tdata->SetBranchAddress("jet1_pt",&jet2pt_data);
	tdata->SetBranchAddress("jet0_eta",&jet1eta_data);
	tdata->SetBranchAddress("jet1_eta",&jet2eta_data);
	
	tmc->SetBranchAddress("jet0_pt",&jet1pt_mc);
	tmc->SetBranchAddress("jet1_pt",&jet2pt_mc);
	tmc->SetBranchAddress("jet0_eta",&jet1eta_mc);
	tmc->SetBranchAddress("jet1_eta",&jet2eta_mc);
	
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
	
	TH1D *heta1data = new TH1D("heta1data","Calibration sample dijet - #eta(jet_{1}) - data;#eta(jet_{1});A.U.",24,2,4.4);
	TH1D *heta2data = new TH1D("heta2data","Calibration sample dijet - #eta(jet_{2}) - data;#eta(jet_{2});A.U.",24,2,4.4);
	TH1D *heta1mc = new TH1D("heta1mc","Calibration sample dijet - #eta(jet_{1}) - MC;#eta(jet_{1});A.U.",24,2,4.4);
	TH1D *heta2mc = new TH1D("heta2mc","Calibration sample dijet - #eta(jet_{2}) - MC;#eta(jet_{2});A.U.",24,2,4.4);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double Dphi_data = acos(((jet1px_data*jet2px_data) + (jet1py_data*jet2py_data))/(pow((((jet1px_data*jet1px_data) + (jet1py_data*jet1py_data))*((jet2px_data*jet2px_data) + (jet2py_data*jet2py_data))),0.5)));
		
		if (jet1pt_data < jetpt_min || jet1pt_data > jetpt_max || jet2pt_data < jetpt_min || jet2pt_data > jetpt_max) continue; // cut in pT of jets
		
		if (Dphi_data < 2.8) continue; // cut in DeltaPhi
		
		heta1data->Fill(jet1eta_data);
		heta2data->Fill(jet2eta_data);
		
	}
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
			
		double Dphi_mc = acos(((jet1px_mc*jet2px_mc) + (jet1py_mc*jet2py_mc))/(pow((((jet1px_mc*jet1px_mc) + (jet1py_mc*jet1py_mc))*((jet2px_mc*jet2px_mc) + (jet2py_mc*jet2py_mc))),0.5)));
		
		if (jet1pt_mc < jetpt_min || jet1pt_mc > jetpt_max || jet2pt_mc < jetpt_min || jet2pt_mc > jetpt_max) continue; // cut in pT of jets

		if (Dphi_mc < 2.8) continue; // cut in DeltaPhi
			
		heta1mc->Fill(jet1eta_mc);	
		heta2mc->Fill(jet2eta_mc);
		
	}
	
	heta1data->Scale(1/(heta1data->Integral()));
	heta2data->Scale(1/(heta2data->Integral()));
	heta1mc->Scale(1/(heta1mc->Integral()));
	heta2mc->Scale(1/(heta2mc->Integral()));
	
	TCanvas *c5 = new TCanvas("c5","etajet1");
	
	heta1data->SetLineColor(kRed);
	heta1mc->SetTitle("Calibration sample dijet - #eta(jet_{1});#eta(jet_{1});A.U.");
	gStyle->SetOptStat("");
	
	c5->cd();
	
	heta1mc->Draw();
	heta1data->Draw("SAME");
	
	TLegend leg5(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg5.AddEntry("heta1data","Data");
	leg5.AddEntry("heta1mc","MC");
	leg5.DrawClone("SAME");
	
	TCanvas *c6 = new TCanvas("c6","etajet2");
	
	heta2data->SetLineColor(kRed);
	heta2mc->SetTitle("Calibration sample dijet - #eta(jet_{2});#eta(jet_{2});A.U.");
	gStyle->SetOptStat("");
	
	c6->cd();
	
	heta2mc->Draw();
	heta2data->Draw("SAME");
	
	TLegend leg6(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg6.AddEntry("heta2data","Data");
	leg6.AddEntry("heta2mc","MC");
	leg6.DrawClone("SAME");
	
}








void GetJetphi() {
	
	double jetpt_min = 2e4; // MeV - jetpt_data >= 20 GeV
	double jetpt_max = 10e4; // MeV
	
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
	
	tdata->SetBranchAddress("jet0_pt",&jet1pt_data);
	tdata->SetBranchAddress("jet1_pt",&jet2pt_data);
	
	tmc->SetBranchAddress("jet0_pt",&jet1pt_mc);
	tmc->SetBranchAddress("jet1_pt",&jet2pt_mc);
	
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
	
	TH1D *hphi1data = new TH1D("hphi1data","Calibration sample dijet - #phi(jet_{1}) - data;#phi(jet_{1});A.U.",32,0,3.2);
	TH1D *hphi2data = new TH1D("hphi2data","Calibration sample dijet - #phi(jet_{2}) - data;#phi(jet_{2});A.U.",32,0,3.2);
	TH1D *hphi1mc = new TH1D("hphi1mc","Calibration sample dijet - #phi(jet_{1}) - MC;#phi(jet_{1});A.U.",32,0,3.2);
	TH1D *hphi2mc = new TH1D("hphi2mc","Calibration sample dijet - #phi(jet_{2}) - MC;#phi(jet_{2});A.U.",32,0,3.2);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double Dphi_data = acos(((jet1px_data*jet2px_data) + (jet1py_data*jet2py_data))/(pow((((jet1px_data*jet1px_data) + (jet1py_data*jet1py_data))*((jet2px_data*jet2px_data) + (jet2py_data*jet2py_data))),0.5)));
		
		if (jet1pt_data < jetpt_min || jet1pt_data > jetpt_max || jet2pt_data < jetpt_min || jet2pt_data > jetpt_max) continue; // cut in pT of jets
		
		if (Dphi_data < 2.8) continue; // cut in DeltaPhi
		
		hphi1data->Fill(acos(jet1px_data/(pow((jet1px_data*jet1px_data+jet1py_data*jet1py_data),0.5))));
		hphi2data->Fill(acos(jet2px_data/(pow((jet2px_data*jet2px_data+jet2py_data*jet2py_data),0.5))));
		
	}
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
			
		double Dphi_mc = acos(((jet1px_mc*jet2px_mc) + (jet1py_mc*jet2py_mc))/(pow((((jet1px_mc*jet1px_mc) + (jet1py_mc*jet1py_mc))*((jet2px_mc*jet2px_mc) + (jet2py_mc*jet2py_mc))),0.5)));
		
		if (jet1pt_mc < jetpt_min || jet1pt_mc > jetpt_max || jet2pt_mc < jetpt_min || jet2pt_mc > jetpt_max) continue; // cut in pT of jets

		if (Dphi_mc < 2.8) continue; // cut in DeltaPhi
			
		hphi1mc->Fill(acos(jet1px_mc/(pow((jet1px_mc*jet1px_mc+jet1py_mc*jet1py_mc),0.5))));
		hphi2mc->Fill(acos(jet2px_mc/(pow((jet2px_mc*jet2px_mc+jet2py_mc*jet2py_mc),0.5))));
		
	}
	
	hphi1data->Scale(1/(hphi1data->Integral()));
	hphi2data->Scale(1/(hphi2data->Integral()));
	hphi1mc->Scale(1/(hphi1mc->Integral()));
	hphi2mc->Scale(1/(hphi2mc->Integral()));
	
	TCanvas *c7 = new TCanvas("c7","phijet1");
	
	hphi1data->SetLineColor(kRed);
	hphi1data->SetTitle("Calibration sample dijet - #phi(jet_{1});#phi(jet_{1});A.U.");
	gStyle->SetOptStat("");
	
	c7->cd();
	
	hphi1data->Draw();
	hphi1mc->Draw("SAME");
	
	TLegend leg7(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg7.AddEntry("hphi1data","Data");
	leg7.AddEntry("hphi1mc","MC");
	leg7.DrawClone("SAME");
	
	TCanvas *c8 = new TCanvas("c8","phijet2");
	
	hphi2data->SetLineColor(kRed);
	hphi2data->SetTitle("Calibration sample dijet - #phi(jet_{2});#phi(jet_{2});A.U.");
	gStyle->SetOptStat("");
	
	c8->cd();
	
	hphi2data->Draw();
	hphi2mc->Draw("SAME");
	
	TLegend leg8(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg8.AddEntry("hphi2data","Data");
	leg8.AddEntry("hphi2mc","MC");
	leg8.DrawClone("SAME");
	
}
