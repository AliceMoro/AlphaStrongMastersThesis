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
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

double L = 1.6;		// integrated luminosity of Run 2 [fb^-1]
double Ngen = 1e6;		// number of generated events

double sigma3 = 4.26194328e12;	// cross section [fb] for pt > 10 GeV, as = 0.110
double sigma4 = 5.3888668e12;	// cross section [fb] for pt > 10 GeV, as = 0.130
double sigma5 = 4.6643252e12;	// cross section [fb] for pt > 10 GeV, as = 0.118

double F3 = ((L*sigma3)/Ngen); // scaling factor for pt > 10 GeV, as = 0.110
double F4 = ((L*sigma4)/Ngen); // scaling factor for pt > 10 GeV, as = 0.130
double F5 = ((L*sigma5)/Ngen); // scaling factor for pt > 10 GeV, as = 0.118

void pt1Norm() {
	
	TCanvas *c1 = new TCanvas("c1","pt1");
		
	cout << "\nScaling factors:\tF3 (as = 0.110) = " << F3 << "\n\t\t\tF4 (as = 0.130) = " << F4 << "\n\t\t\tF5 (as = 0.180) = " << F5 << endl;
	
	TFile *f3 = new TFile("jets_ct10_10gev_0110.root","READ");
	TFile *f4 = new TFile("jets_ct10_10gev_0130.root","READ");
	TFile *f5 = new TFile("jets_ct10_10gev_0118.root","READ");
	
	TTree *t3 = (TTree*)f3->Get("jets_ct10_10gev_0110");
	TTree *t4 = (TTree*)f4->Get("jets_ct10_10gev_0130");
	TTree *t5 = (TTree*)f5->Get("jets_ct10_10gev_0118");
	
	int N3 = t3->GetEntries();
	int N4 = t4->GetEntries();
	int N5 = t5->GetEntries();
	
	cout << "\nNumber of entries:\tN3 = " << N3 << "\n\t\t\tN4 = " << N4 << "\n\t\t\tN5 = " << N5 << endl;
	
	// variables to store branch content of tree t3
	double jet1pt_3;
	double jet2pt_3;
	double jet3pt_3;
	
	t3->SetBranchAddress("jet1pt",&jet1pt_3);
	t3->SetBranchAddress("jet2pt",&jet2pt_3);
	t3->SetBranchAddress("jet3pt",&jet3pt_3);
	
	// variables to store branch content of tree t4
	double jet1pt_4;
	double jet2pt_4;
	double jet3pt_4;
	
	t4->SetBranchAddress("jet1pt",&jet1pt_4);
	t4->SetBranchAddress("jet2pt",&jet2pt_4);
	t4->SetBranchAddress("jet3pt",&jet3pt_4);
	
	// variables to store branch content of tree t5
	double jet1pt_5;
	double jet2pt_5;
	double jet3pt_5;
	
	t5->SetBranchAddress("jet1pt",&jet1pt_5);
	t5->SetBranchAddress("jet2pt",&jet2pt_5);
	t5->SetBranchAddress("jet3pt",&jet3pt_5);
	
	// histograms of distributions
	TH1D *hpt1_3 = new TH1D("hpt1_3","p_{T}(jet_{1}) - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;p_{T}(jet_{1}) [GeV];A.U.",50,10,110);
	TH1D *hpt1_4 = new TH1D("hpt1_4","p_{T}(jet_{1}) - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;p_{T}(jet_{1}) [GeV];A.U.",50,10,110);
	TH1D *hpt1_5 = new TH1D("hpt1_5","p_{T}(jet_{1}) - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;p_{T}(jet_{1}) [GeV];A.U.",50,10,110);
	
	for (int n=0; n<N3; n++) {

		t3->GetEntry(n);
		
		if (jet1pt_3 < 20 || jet2pt_3 < 20 || jet3pt_3 < 20) continue;
					
		hpt1_3->Fill(jet1pt_3);
		
	}
	
	for (int m=0; m<N4; m++) {

		t4->GetEntry(m);
		
		if (jet1pt_4 < 20 || jet2pt_4 < 20 || jet3pt_4 < 20) continue;
					
		hpt1_4->Fill(jet1pt_4);
		
	}
	
	for (int p=0; p<N5; p++) {

		t5->GetEntry(p);
		
		if (jet1pt_5 < 20 || jet2pt_5 < 20 || jet3pt_5 < 20) continue;
					
		hpt1_5->Fill(jet1pt_5);
		
	}
	
	cout << "\nNumber of entries of the histogram 3: " << hpt1_3->GetEntries() << endl;
	cout << "Number of entries of the histogram 4: " << hpt1_4->GetEntries() << endl;
	cout << "Number of entries of the histogram 5: " << hpt1_5->GetEntries() << endl;
	cout << endl;
	
	hpt1_3->Scale(F3);
	hpt1_4->Scale(F4);
	hpt1_5->Scale(F5);
	
	hpt1_3->Scale(1/(hpt1_3->Integral()));
	hpt1_4->Scale(1/(hpt1_4->Integral()));
	hpt1_5->Scale(1/(hpt1_5->Integral()));
	
	hpt1_3->GetXaxis()->SetRangeUser(10,110);
	hpt1_3->GetXaxis()->SetNdivisions(511);
	hpt1_3->SetLineWidth(1);
	hpt1_3->SetMarkerStyle(20);
	hpt1_3->SetMarkerSize(0.8);
	hpt1_3->SetMarkerColor(kBlue+2);
	
	hpt1_4->GetXaxis()->SetRangeUser(10,110);
	hpt1_4->GetXaxis()->SetNdivisions(511);
	hpt1_4->SetLineWidth(1);
	hpt1_4->SetMarkerStyle(20);
	hpt1_4->SetMarkerSize(0.8);
	hpt1_4->SetLineColor(kRed);
	hpt1_4->SetMarkerColor(kRed);
	
	hpt1_5->GetXaxis()->SetRangeUser(10,110);
	hpt1_5->GetXaxis()->SetNdivisions(511);
	hpt1_5->SetLineWidth(1);
	hpt1_5->SetMarkerStyle(20);
	hpt1_5->SetMarkerSize(0.8);
	hpt1_5->SetLineColor(kGreen+1);
	hpt1_5->SetMarkerColor(kGreen+1);
	
	gStyle->SetOptStat("");
	
	c1->cd();
	
	hpt1_5->Draw();
	hpt1_3->Draw("SAME");
	hpt1_4->Draw("SAME");
	
	TLegend leg2(.5,.7,.9,.9); // (x1,y1,x2,y2)
	
	leg2.SetFillColor(0);
	leg2.AddEntry("hpt1_3","CT10nnlo #alpha_{s}=0.1100");
	leg2.AddEntry("hpt1_5","CT10nnlo #alpha_{s}=0.1180");
	leg2.AddEntry("hpt1_4","CT10nnlo #alpha_{s}=0.1300");
	leg2.DrawClone("SAME");
	
}

void pt2Norm() {
	
	TCanvas *c2 = new TCanvas("c2","pt2");
	
	cout << "\nScaling factors:\tF3 (as = 0.110) = " << F3 << "\n\t\t\tF4 (as = 0.130) = " << F4 << "\n\t\t\tF5 (as = 0.180) = " << F5 << endl;
	
	TFile *f3 = new TFile("jets_ct10_10gev_0110.root","READ");
	TFile *f4 = new TFile("jets_ct10_10gev_0130.root","READ");
	TFile *f5 = new TFile("jets_ct10_10gev_0118.root","READ");
	
	TTree *t3 = (TTree*)f3->Get("jets_ct10_10gev_0110");
	TTree *t4 = (TTree*)f4->Get("jets_ct10_10gev_0130");
	TTree *t5 = (TTree*)f5->Get("jets_ct10_10gev_0118");
	
	int N3 = t3->GetEntries();
	int N4 = t4->GetEntries();
	int N5 = t5->GetEntries();
	
	cout << "\nNumber of entries:\tN3 = " << N3 << "\n\t\t\tN4 = " << N4 << "\n\t\t\tN5 = " << N5 << endl;
	
	// variables to store branch content of tree t3
	double jet1pt_3;
	double jet2pt_3;
	double jet3pt_3;
	
	t3->SetBranchAddress("jet1pt",&jet1pt_3);
	t3->SetBranchAddress("jet2pt",&jet2pt_3);
	t3->SetBranchAddress("jet3pt",&jet3pt_3);
	
	// variables to store branch content of tree t4
	double jet1pt_4;
	double jet2pt_4;
	double jet3pt_4;
	
	t4->SetBranchAddress("jet1pt",&jet1pt_4);
	t4->SetBranchAddress("jet2pt",&jet2pt_4);
	t4->SetBranchAddress("jet3pt",&jet3pt_4);
	
	// variables to store branch content of tree t5
	double jet1pt_5;
	double jet2pt_5;
	double jet3pt_5;
	
	t5->SetBranchAddress("jet1pt",&jet1pt_5);
	t5->SetBranchAddress("jet2pt",&jet2pt_5);
	t5->SetBranchAddress("jet3pt",&jet3pt_5);
	
	// histograms of distributions
	TH1D *hpt2_3 = new TH1D("hpt2_3","p_{T}(jet_{2}) - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;p_{T}(jet_{2}) [GeV];A.U.",50,10,110);
	TH1D *hpt2_4 = new TH1D("hpt2_4","p_{T}(jet_{2}) - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;p_{T}(jet_{2}) [GeV];A.U.",50,10,110);
	TH1D *hpt2_5 = new TH1D("hpt2_5","p_{T}(jet_{2}) - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;p_{T}(jet_{2}) [GeV];A.U.",50,10,110);
	
	for (int n=0; n<N3; n++) {

		t3->GetEntry(n);
		
		if (jet1pt_3 < 20 || jet2pt_3 < 20 || jet3pt_3 < 20) continue;
					
		hpt2_3->Fill(jet2pt_3);
		
	}
	
	for (int m=0; m<N4; m++) {

		t4->GetEntry(m);
		
		if (jet1pt_4 < 20 || jet2pt_4 < 20 || jet3pt_4 < 20) continue;
					
		hpt2_4->Fill(jet2pt_4);
		
	}
	
	for (int p=0; p<N5; p++) {

		t5->GetEntry(p);
		
		if (jet1pt_5 < 20 || jet2pt_5 < 20 || jet3pt_5 < 20) continue;
					
		hpt2_5->Fill(jet2pt_5);
		
	}
	
	cout << "\nNumber of entries of the histogram 3: " << hpt2_3->GetEntries() << endl;
	cout << "Number of entries of the histogram 4: " << hpt2_4->GetEntries() << endl;
	cout << "Number of entries of the histogram 5: " << hpt2_5->GetEntries() << endl;
	cout << endl;
	
	hpt2_3->Scale(F3);
	hpt2_4->Scale(F4);
	hpt2_5->Scale(F5);
	
	hpt2_3->Scale(1/(hpt2_3->Integral()));
	hpt2_4->Scale(1/(hpt2_4->Integral()));
	hpt2_5->Scale(1/(hpt2_5->Integral()));
	
	hpt2_3->GetXaxis()->SetRangeUser(10,110);
	hpt2_3->GetXaxis()->SetNdivisions(511);
	hpt2_3->SetLineWidth(1);
	hpt2_3->SetMarkerStyle(20);
	hpt2_3->SetMarkerSize(0.8);
	hpt2_3->SetMarkerColor(kBlue+2);
	
	hpt2_4->GetXaxis()->SetRangeUser(10,110);
	hpt2_4->GetXaxis()->SetNdivisions(511);
	hpt2_4->SetLineWidth(1);
	hpt2_4->SetMarkerStyle(20);
	hpt2_4->SetMarkerSize(0.8);
	hpt2_4->SetLineColor(kRed);
	hpt2_4->SetMarkerColor(kRed);
	
	hpt2_5->GetXaxis()->SetRangeUser(10,110);
	hpt2_5->GetXaxis()->SetNdivisions(511);
	hpt2_5->SetLineWidth(1);
	hpt2_5->SetMarkerStyle(20);
	hpt2_5->SetMarkerSize(0.8);
	hpt2_5->SetLineColor(kGreen+1);
	hpt2_5->SetMarkerColor(kGreen+1);
	
	gStyle->SetOptStat("");
	
	c2->cd();
	
	hpt2_4->Draw();
	hpt2_3->Draw("SAME");
	hpt2_5->Draw("SAME");
	
	TLegend leg2(.5,.7,.9,.9); // (x1,y1,x2,y2)
	
	leg2.SetFillColor(0);
	leg2.AddEntry("hpt2_3","CT10nnlo #alpha_{s}=0.1100");
	leg2.AddEntry("hpt2_5","CT10nnlo #alpha_{s}=0.1180");
	leg2.AddEntry("hpt2_4","CT10nnlo #alpha_{s}=0.1300");
	leg2.DrawClone("SAME");
	
}

void QNorm() {
	
	TCanvas *c3 = new TCanvas("c3","Q");
	
	cout << "\nScaling factors:\tF3 = " << F3 << "\n\t\t\tF4 = " << F4 << "\n\t\t\tF5 = " << F5 << endl;
	
	TFile *f3 = new TFile("jets_ct10_10gev_0110.root","READ");
	TFile *f4 = new TFile("jets_ct10_10gev_0130.root","READ");
	TFile *f5 = new TFile("jets_ct10_10gev_0118.root","READ");
	
	TTree *t3 = (TTree*)f3->Get("jets_ct10_10gev_0110");
	TTree *t4 = (TTree*)f4->Get("jets_ct10_10gev_0130");
	TTree *t5 = (TTree*)f5->Get("jets_ct10_10gev_0118");
	
	int N3 = t3->GetEntries();
	int N4 = t4->GetEntries();
	int N5 = t5->GetEntries();
	
	cout << "\nNumber of entries:\tN3 = " << N3 << "\n\t\t\tN4 = " << N4 << "\n\t\t\tN5 = " << N5 << endl;
	
	// variables to store branch content of tree t3
	double jet1pt_3;
	double jet2pt_3;
	double jet3pt_3;
	double Q_3;
	
	t3->SetBranchAddress("jet1pt",&jet1pt_3);
	t3->SetBranchAddress("jet2pt",&jet2pt_3);
	t3->SetBranchAddress("jet3pt",&jet3pt_3);
	t3->SetBranchAddress("Q",&Q_3);
	
	// variables to store branch content of tree t4
	double jet1pt_4;
	double jet2pt_4;
	double jet3pt_4;
	double Q_4;
	
	t4->SetBranchAddress("jet1pt",&jet1pt_4);
	t4->SetBranchAddress("jet2pt",&jet2pt_4);
	t4->SetBranchAddress("jet3pt",&jet3pt_4);
	t4->SetBranchAddress("Q",&Q_4);
	t3->SetBranchAddress("Q",&Q_3);
	
	// variables to store branch content of tree t5
	double jet1pt_5;
	double jet2pt_5;
	double jet3pt_5;
	double Q_5;
	
	t5->SetBranchAddress("jet1pt",&jet1pt_5);
	t5->SetBranchAddress("jet2pt",&jet2pt_5);
	t5->SetBranchAddress("jet3pt",&jet3pt_5);
	t5->SetBranchAddress("Q",&Q_5);
	
	// histograms of distributions
	TH1D *hQ_3 = new TH1D("hQ_3","Q - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;Q [GeV];A.U.",50,10,110);
	TH1D *hQ_4 = new TH1D("hQ_4","Q - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;Q [GeV];A.U.",50,10,110);
	TH1D *hQ_5 = new TH1D("hQ_5","Q - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;Q [GeV];A.U.",50,10,110);
	
	for (int n=0; n<N3; n++) {

		t3->GetEntry(n);
		
		if (jet1pt_3 < 20 || jet2pt_3 < 20 || jet3pt_3 < 20) continue;
					
		hQ_3->Fill(Q_3);
		
	}
	
	for (int m=0; m<N4; m++) {

		t4->GetEntry(m);
		
		if (jet1pt_4 < 20 || jet2pt_4 < 20 || jet3pt_4 < 20) continue;
					
		hQ_4->Fill(Q_4);
		
	}
	
	for (int p=0; p<N5; p++) {

		t5->GetEntry(p);
		
		if (jet1pt_5 < 20 || jet2pt_5 < 20 || jet3pt_5 < 20) continue;
					
		hQ_5->Fill(Q_5);
		
	}
	
	cout << "\nNumber of entries of the histogram 3: " << hQ_3->GetEntries() << endl;
	cout << "Number of entries of the histogram 4: " << hQ_4->GetEntries() << endl;
	cout << "Number of entries of the histogram 5: " << hQ_5->GetEntries() << endl;
	cout << endl;
	
	hQ_3->Scale(F3);
	hQ_4->Scale(F4);
	hQ_5->Scale(F5);
	
	hQ_3->Scale(1/(hQ_3->Integral()));
	hQ_4->Scale(1/(hQ_4->Integral()));
	hQ_5->Scale(1/(hQ_5->Integral()));
	
	hQ_3->GetXaxis()->SetRangeUser(10,110);
	hQ_3->GetXaxis()->SetNdivisions(511);
	hQ_3->SetLineWidth(1);
	hQ_3->SetMarkerStyle(20);
	hQ_3->SetMarkerSize(0.8);
	hQ_3->SetMarkerColor(kBlue+2);
	
	hQ_4->GetXaxis()->SetRangeUser(10,110);
	hQ_4->GetXaxis()->SetNdivisions(511);
	hQ_4->SetLineWidth(1);
	hQ_4->SetMarkerStyle(20);
	hQ_4->SetMarkerSize(0.8);
	hQ_4->SetMarkerColor(kRed);
	hQ_4->SetLineColor(kRed);
	
	hQ_5->GetXaxis()->SetRangeUser(10,110);
	hQ_5->GetXaxis()->SetNdivisions(511);
	hQ_5->SetLineWidth(1);
	hQ_5->SetMarkerStyle(20);
	hQ_5->SetMarkerSize(0.8);
	hQ_5->SetMarkerColor(kGreen+1);
	hQ_5->SetLineColor(kGreen+1);
	
	gStyle->SetOptStat("");
	
	c3->cd();
	
	hQ_3->Draw();
	hQ_4->Draw("SAME");
	hQ_5->Draw("SAME");
	
	TLegend leg2(.5,.7,.9,.9); // (x1,y1,x2,y2)
	
	leg2.SetFillColor(0);
	leg2.AddEntry("hQ_3","CT10nnlo #alpha_{s}=0.1100");
	leg2.AddEntry("hQ_5","CT10nnlo #alpha_{s}=0.1180");
	leg2.AddEntry("hQ_4","CT10nnlo #alpha_{s}=0.1300");
	leg2.DrawClone("SAME");
	
}

void m12Norm() {
	
	TCanvas *c4 = new TCanvas("c4","m12");
	
	cout << "\nScaling factors:\tF3 = " << F3 << "\n\t\t\tF4 = " << F4 << "\n\t\t\tF5 = " << F5 << endl;
	
	TFile *f3 = new TFile("jets_ct10_10gev_0110.root","READ");
	TFile *f4 = new TFile("jets_ct10_10gev_0130.root","READ");
	TFile *f5 = new TFile("jets_ct10_10gev_0118.root","READ");
	
	TTree *t3 = (TTree*)f3->Get("jets_ct10_10gev_0110");
	TTree *t4 = (TTree*)f4->Get("jets_ct10_10gev_0130");
	TTree *t5 = (TTree*)f5->Get("jets_ct10_10gev_0118");
	
	int N3 = t3->GetEntries();
	int N4 = t4->GetEntries();
	int N5 = t5->GetEntries();
	
	cout << "\nNumber of entries:\tN3 = " << N3 << "\n\t\t\tN4 = " << N4 << "\n\t\t\tN5 = " << N5 << endl;
	
	// variables to store branch content of tree t3
	double jet1pt_3;
	double jet2pt_3;
	double jet3pt_3;
	double jet1mass_3;
	double jet2mass_3;
	double jet1px_3;
	double jet1py_3;
	double jet1pz_3;
	double jet2px_3;
	double jet2py_3;
	double jet2pz_3;
	
	t3->SetBranchAddress("jet1pt",&jet1pt_3);
	t3->SetBranchAddress("jet2pt",&jet2pt_3);
	t3->SetBranchAddress("jet3pt",&jet3pt_3);
	t3->SetBranchAddress("jet1mass",&jet1mass_3);
	t3->SetBranchAddress("jet2mass",&jet2mass_3);
	t3->SetBranchAddress("jet1px",&jet1px_3);
	t3->SetBranchAddress("jet1py",&jet1py_3);
	t3->SetBranchAddress("jet1pz",&jet1pz_3);
	t3->SetBranchAddress("jet2px",&jet2px_3);
	t3->SetBranchAddress("jet2py",&jet2py_3);
	t3->SetBranchAddress("jet2pz",&jet2pz_3);
	
	// array to store m12 values
	double m12_3[619972] = {0.0};
	
	// variables to store branch content of tree t4
	double jet1pt_4;
	double jet2pt_4;
	double jet3pt_4;
	double jet1mass_4;
	double jet2mass_4;
	double jet1px_4;
	double jet1py_4;
	double jet1pz_4;
	double jet2px_4;
	double jet2py_4;
	double jet2pz_4;
	
	t4->SetBranchAddress("jet1pt",&jet1pt_4);
	t4->SetBranchAddress("jet2pt",&jet2pt_4);
	t4->SetBranchAddress("jet3pt",&jet3pt_4);
	t4->SetBranchAddress("jet1mass",&jet1mass_4);
	t4->SetBranchAddress("jet2mass",&jet2mass_4);
	t4->SetBranchAddress("jet1px",&jet1px_4);
	t4->SetBranchAddress("jet1py",&jet1py_4);
	t4->SetBranchAddress("jet1pz",&jet1pz_4);
	t4->SetBranchAddress("jet2px",&jet2px_4);
	t4->SetBranchAddress("jet2py",&jet2py_4);
	t4->SetBranchAddress("jet2pz",&jet2pz_4);
	
	// array to store m12 values
	double m12_4[622422] = {0.0};
	
	// variables to store branch content of tree t5
	double jet1pt_5;
	double jet2pt_5;
	double jet3pt_5;
	double jet1mass_5;
	double jet2mass_5;
	double jet1px_5;
	double jet1py_5;
	double jet1pz_5;
	double jet2px_5;
	double jet2py_5;
	double jet2pz_5;
	
	t5->SetBranchAddress("jet1pt",&jet1pt_5);
	t5->SetBranchAddress("jet2pt",&jet2pt_5);
	t5->SetBranchAddress("jet3pt",&jet3pt_5);
	t5->SetBranchAddress("jet1mass",&jet1mass_5);
	t5->SetBranchAddress("jet2mass",&jet2mass_5);
	t5->SetBranchAddress("jet1px",&jet1px_5);
	t5->SetBranchAddress("jet1py",&jet1py_5);
	t5->SetBranchAddress("jet1pz",&jet1pz_5);
	t5->SetBranchAddress("jet2px",&jet2px_5);
	t5->SetBranchAddress("jet2py",&jet2py_5);
	t5->SetBranchAddress("jet2pz",&jet2pz_5);
	
	// array to store m12 values
	double m12_5[621619] = {0.0};
	
	// histograms of distributions
	TH1D *hm12_3 = new TH1D("hm12_3","m_{12} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;m_{12} [GeV];A.U.",40,0,240);
	TH1D *hm12_4 = new TH1D("hm12_4","m_{12} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;m_{12} [GeV];A.U.",40,0,240);
	TH1D *hm12_5 = new TH1D("hm12_5","m_{12} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;m_{12} [GeV];A.U.",40,0,240);
	
	for (int n=0; n<N3; n++) {

		t3->GetEntry(n);
		
		m12_3[n] = pow(((jet1mass_3*jet1mass_3) + (jet2mass_3*jet2mass_3) + (2*pow((((jet1mass_3*jet1mass_3) + (jet1px_3*jet1px_3) + (jet1py_3*jet1py_3) + (jet1pz_3*jet1pz_3))*((jet2mass_3*jet2mass_3) + (jet2px_3*jet2px_3) + (jet2py_3*jet2py_3) + (jet2pz_3*jet2pz_3))),0.5)) - (2*((jet1px_3*jet2px_3) + (jet1py_3*jet2py_3) + (jet1pz_3*jet2pz_3)))),0.5);
		
		if (jet1pt_3 < 20 || jet2pt_3 < 20 || jet3pt_3 < 20) continue;
					
		hm12_3->Fill(m12_3[n]);
		
	}
	
	for (int m=0; m<N4; m++) {

		t4->GetEntry(m);
		
		m12_4[m] = pow(((jet1mass_4*jet1mass_4) + (jet2mass_4*jet2mass_4) + (2*pow((((jet1mass_4*jet1mass_4) + (jet1px_4*jet1px_4) + (jet1py_4*jet1py_4) + (jet1pz_4*jet1pz_4))*((jet2mass_4*jet2mass_4) + (jet2px_4*jet2px_4) + (jet2py_4*jet2py_4) + (jet2pz_4*jet2pz_4))),0.5)) - (2*((jet1px_4*jet2px_4) + (jet1py_4*jet2py_4) + (jet1pz_4*jet2pz_4)))),0.5);
		
		if (jet1pt_4 < 20 || jet2pt_4 < 20 || jet3pt_4 < 20) continue;
					
		hm12_4->Fill(m12_4[m]);
		
	}
	
	for (int p=0; p<N5; p++) {

		t5->GetEntry(p);
		
		m12_5[p] = pow(((jet1mass_5*jet1mass_5) + (jet2mass_5*jet2mass_5) + (2*pow((((jet1mass_5*jet1mass_5) + (jet1px_5*jet1px_5) + (jet1py_5*jet1py_5) + (jet1pz_5*jet1pz_5))*((jet2mass_5*jet2mass_5) + (jet2px_5*jet2px_5) + (jet2py_5*jet2py_5) + (jet2pz_5*jet2pz_5))),0.5)) - (2*((jet1px_5*jet2px_5) + (jet1py_5*jet2py_5) + (jet1pz_5*jet2pz_5)))),0.5);
		
		if (jet1pt_5 < 20 || jet2pt_5 < 20 || jet3pt_5 < 20) continue;
					
		hm12_5->Fill(m12_5[p]);
		
	}
	
	cout << "\nNumber of entries of the histogram 3: " << hm12_3->GetEntries() << endl;
	cout << "Number of entries of the histogram 4: " << hm12_4->GetEntries() << endl;
	cout << "Number of entries of the histogram 5: " << hm12_5->GetEntries() << endl;
	cout << endl;
	
	hm12_3->Scale(F3);
	hm12_4->Scale(F4);
	hm12_5->Scale(F5);
	
	hm12_3->Scale(1/(hm12_3->Integral()));
	hm12_4->Scale(1/(hm12_4->Integral()));
	hm12_5->Scale(1/(hm12_5->Integral()));
	
	hm12_3->GetXaxis()->SetRangeUser(0,240);
	hm12_3->GetXaxis()->SetNdivisions(511);
	hm12_3->SetLineWidth(1);
	hm12_4->SetMarkerStyle(20);
	hm12_4->SetMarkerSize(0.8);
	hm12_4->SetMarkerColor(kBlue+2);
	
	hm12_4->GetXaxis()->SetRangeUser(0,240);
	hm12_4->GetXaxis()->SetNdivisions(511);
	hm12_4->SetLineWidth(1);
	hm12_4->SetLineColor(kRed);
	hm12_4->SetMarkerStyle(20);
	hm12_4->SetMarkerSize(0.8);
	hm12_4->SetMarkerColor(kRed);
	
	hm12_5->GetXaxis()->SetRangeUser(0,240);
	hm12_5->GetXaxis()->SetNdivisions(511);
	hm12_5->SetLineWidth(1);
	hm12_5->SetLineColor(kGreen+1);
	hm12_5->SetMarkerStyle(20);
	hm12_5->SetMarkerSize(0.8);
	hm12_5->SetMarkerColor(kGreen+1);
	
	gStyle->SetOptStat("");
	
	c4->cd();
	
	hm12_4->Draw();
	hm12_5->Draw("SAME");
	hm12_3->Draw("SAME");
	
	TLegend leg2(.5,.7,.9,.9); // (x1,y1,x2,y2)
	
	leg2.SetFillColor(0);
	leg2.AddEntry("hm12_3","CT10nnlo #alpha_{s}=0.1100");
	leg2.AddEntry("hm12_5","CT10nnlo #alpha_{s}=0.1180");
	leg2.AddEntry("hm12_4","CT10nnlo #alpha_{s}=0.1300");
	leg2.DrawClone("SAME");
	
}

void mtotNorm() {
	
	TCanvas *c5 = new TCanvas("c5","mtot");
	
	cout << "\nScaling factors:\tF3 = " << F3 << "\n\t\t\tF4 = " << F4 << "\n\t\t\tF5 = " << F5 << endl;
	
	TFile *f3 = new TFile("jets_ct10_10gev_0110.root","READ");
	TFile *f4 = new TFile("jets_ct10_10gev_0130.root","READ");
	TFile *f5 = new TFile("jets_ct10_10gev_0118.root","READ");
	
	TTree *t3 = (TTree*)f3->Get("jets_ct10_10gev_0110");
	TTree *t4 = (TTree*)f4->Get("jets_ct10_10gev_0130");
	TTree *t5 = (TTree*)f5->Get("jets_ct10_10gev_0118");
	
	int N3 = t3->GetEntries();
	int N4 = t4->GetEntries();
	int N5 = t5->GetEntries();
	
	cout << "\nNumber of entries:\tN3 = " << N3 << "\n\t\t\tN4 = " << N4 << "\n\t\t\tN5 = " << N5 << endl;
	
	// variables to store branch content of tree t3
	double jet1pt_3;
	double jet2pt_3;
	double jet3pt_3;
	double jet1mass_3;
	double jet2mass_3;
	double jet3mass_3;
	double jet1px_3;
	double jet1py_3;
	double jet1pz_3;
	double jet2px_3;
	double jet2py_3;
	double jet2pz_3;
	double jet3px_3;
	double jet3py_3;
	double jet3pz_3;
	
	t3->SetBranchAddress("jet1pt",&jet1pt_3);
	t3->SetBranchAddress("jet2pt",&jet2pt_3);
	t3->SetBranchAddress("jet3pt",&jet3pt_3);
	t3->SetBranchAddress("jet1mass",&jet1mass_3);
	t3->SetBranchAddress("jet2mass",&jet2mass_3);
	t3->SetBranchAddress("jet3mass",&jet3mass_3);
	t3->SetBranchAddress("jet1px",&jet1px_3);
	t3->SetBranchAddress("jet1py",&jet1py_3);
	t3->SetBranchAddress("jet1pz",&jet1pz_3);
	t3->SetBranchAddress("jet2px",&jet2px_3);
	t3->SetBranchAddress("jet2py",&jet2py_3);
	t3->SetBranchAddress("jet2pz",&jet2pz_3);
	t3->SetBranchAddress("jet3px",&jet3px_3);
	t3->SetBranchAddress("jet3py",&jet3py_3);
	t3->SetBranchAddress("jet3pz",&jet3pz_3);
	
	// array to store mtot values
	double mtot_3[619972] = {0.0};
	
	// variables to store branch content of tree t4
	double jet1pt_4;
	double jet2pt_4;
	double jet3pt_4;
	double jet1mass_4;
	double jet2mass_4;
	double jet3mass_4;
	double jet1px_4;
	double jet1py_4;
	double jet1pz_4;
	double jet2px_4;
	double jet2py_4;
	double jet2pz_4;
	double jet3px_4;
	double jet3py_4;
	double jet3pz_4;
	
	t4->SetBranchAddress("jet1pt",&jet1pt_4);
	t4->SetBranchAddress("jet2pt",&jet2pt_4);
	t4->SetBranchAddress("jet3pt",&jet3pt_4);
	t4->SetBranchAddress("jet1mass",&jet1mass_4);
	t4->SetBranchAddress("jet2mass",&jet2mass_4);
	t4->SetBranchAddress("jet3mass",&jet3mass_4);
	t4->SetBranchAddress("jet1px",&jet1px_4);
	t4->SetBranchAddress("jet1py",&jet1py_4);
	t4->SetBranchAddress("jet1pz",&jet1pz_4);
	t4->SetBranchAddress("jet2px",&jet2px_4);
	t4->SetBranchAddress("jet2py",&jet2py_4);
	t4->SetBranchAddress("jet2pz",&jet2pz_4);
	t4->SetBranchAddress("jet3px",&jet3px_4);
	t4->SetBranchAddress("jet3py",&jet3py_4);
	t4->SetBranchAddress("jet3pz",&jet3pz_4);
	
	// array to store mtot values
	double mtot_4[622422] = {0.0};
	
	// variables to store branch content of tree t5
	double jet1pt_5;
	double jet2pt_5;
	double jet3pt_5;
	double jet1mass_5;
	double jet2mass_5;
	double jet3mass_5;
	double jet1px_5;
	double jet1py_5;
	double jet1pz_5;
	double jet2px_5;
	double jet2py_5;
	double jet2pz_5;
	double jet3px_5;
	double jet3py_5;
	double jet3pz_5;
	
	t5->SetBranchAddress("jet1pt",&jet1pt_5);
	t5->SetBranchAddress("jet2pt",&jet2pt_5);
	t5->SetBranchAddress("jet3pt",&jet3pt_5);
	t5->SetBranchAddress("jet1mass",&jet1mass_5);
	t5->SetBranchAddress("jet2mass",&jet2mass_5);
	t5->SetBranchAddress("jet3mass",&jet3mass_5);
	t5->SetBranchAddress("jet1px",&jet1px_5);
	t5->SetBranchAddress("jet1py",&jet1py_5);
	t5->SetBranchAddress("jet1pz",&jet1pz_5);
	t5->SetBranchAddress("jet2px",&jet2px_5);
	t5->SetBranchAddress("jet2py",&jet2py_5);
	t5->SetBranchAddress("jet2pz",&jet2pz_5);
	t5->SetBranchAddress("jet3px",&jet3px_5);
	t5->SetBranchAddress("jet3py",&jet3py_5);
	t5->SetBranchAddress("jet3pz",&jet3pz_5);
	
	// array to store mtot values
	double mtot_5[621619] = {0.0};
	
	// histograms of distributions
	TH1D *hmtot_3 = new TH1D("hmtot_3","m_{TOT} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;m_{TOT} [GeV];A.U.",48,0,240);
	TH1D *hmtot_4 = new TH1D("hmtot_4","m_{TOT} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;m_{TOT} [GeV];A.U.",48,0,240);
	TH1D *hmtot_5 = new TH1D("hmtot_5","m_{TOT} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;m_{TOT} [GeV];A.U.",48,0,240);
	
	for (int n=0; n<N3; n++) {

		t3->GetEntry(n);
		
		mtot_3[n] = pow(((jet1mass_3*jet1mass_3) + (jet2mass_3*jet2mass_3) + (jet3mass_3*jet3mass_3) + (2*pow((((jet1mass_3*jet1mass_3) + (jet1px_3*jet1px_3) + (jet1py_3*jet1py_3) + (jet1pz_3*jet1pz_3))*((jet2mass_3*jet2mass_3) + (jet2px_3*jet2px_3) + (jet2py_3*jet2py_3) + (jet2pz_3*jet2pz_3))),0.5)) + (2*pow((((jet1mass_3*jet1mass_3) + (jet1px_3*jet1px_3) + (jet1py_3*jet1py_3) + (jet1pz_3*jet1pz_3))*((jet3mass_3*jet3mass_3) + (jet3px_3*jet3px_3) + (jet3py_3*jet3py_3) + (jet3pz_3*jet3pz_3))),0.5)) + (2*pow((((jet2mass_3*jet2mass_3) + (jet2px_3*jet2px_3) + (jet2py_3*jet2py_3) + (jet2pz_3*jet2pz_3))*((jet3mass_3*jet3mass_3) + (jet3px_3*jet3px_3) + (jet3py_3*jet3py_3) + (jet3pz_3*jet3pz_3))),0.5)) - (2*((jet1px_3*jet2px_3) + (jet1py_3*jet2py_3) + (jet1pz_3*jet2pz_3))) - (2*((jet1px_3*jet3px_3) + (jet1py_3*jet3py_3) + (jet1pz_3*jet3pz_3))) - (2*((jet2px_3*jet3px_3) + (jet2py_3*jet3py_3) + (jet2pz_3*jet3pz_3)))),0.5);
		
		if (jet1pt_3 < 20 || jet2pt_3 < 20 || jet3pt_3 < 20) continue;
					
		hmtot_3->Fill(mtot_3[n]);
		
	}
	
	for (int m=0; m<N4; m++) {

		t4->GetEntry(m);
		
		mtot_4[m] = pow(((jet1mass_4*jet1mass_4) + (jet2mass_4*jet2mass_4) + (jet3mass_4*jet3mass_4) + (2*pow((((jet1mass_4*jet1mass_4) + (jet1px_4*jet1px_4) + (jet1py_4*jet1py_4) + (jet1pz_4*jet1pz_4))*((jet2mass_4*jet2mass_4) + (jet2px_4*jet2px_4) + (jet2py_4*jet2py_4) + (jet2pz_4*jet2pz_4))),0.5)) + (2*pow((((jet1mass_4*jet1mass_4) + (jet1px_4*jet1px_4) + (jet1py_4*jet1py_4) + (jet1pz_4*jet1pz_4))*((jet3mass_4*jet3mass_4) + (jet3px_4*jet3px_4) + (jet3py_4*jet3py_4) + (jet3pz_4*jet3pz_4))),0.5)) + (2*pow((((jet2mass_4*jet2mass_4) + (jet2px_4*jet2px_4) + (jet2py_4*jet2py_4) + (jet2pz_4*jet2pz_4))*((jet3mass_4*jet3mass_4) + (jet3px_4*jet3px_4) + (jet3py_4*jet3py_4) + (jet3pz_4*jet3pz_4))),0.5)) - (2*((jet1px_4*jet2px_4) + (jet1py_4*jet2py_4) + (jet1pz_4*jet2pz_4))) - (2*((jet1px_4*jet3px_4) + (jet1py_4*jet3py_4) + (jet1pz_4*jet3pz_4))) - (2*((jet2px_4*jet3px_4) + (jet2py_4*jet3py_4) + (jet2pz_4*jet3pz_4)))),0.5);
		
		if (jet1pt_4 < 20 || jet2pt_4 < 20 || jet3pt_4 < 20) continue;
					
		hmtot_4->Fill(mtot_4[m]);
		
	}
	
	for (int p=0; p<N5; p++) {

		t5->GetEntry(p);
		
		mtot_5[p] = pow(((jet1mass_5*jet1mass_5) + (jet2mass_5*jet2mass_5) + (jet3mass_5*jet3mass_5) + (2*pow((((jet1mass_5*jet1mass_5) + (jet1px_5*jet1px_5) + (jet1py_5*jet1py_5) + (jet1pz_5*jet1pz_5))*((jet2mass_5*jet2mass_5) + (jet2px_5*jet2px_5) + (jet2py_5*jet2py_5) + (jet2pz_5*jet2pz_5))),0.5)) + (2*pow((((jet1mass_5*jet1mass_5) + (jet1px_5*jet1px_5) + (jet1py_5*jet1py_5) + (jet1pz_5*jet1pz_5))*((jet3mass_5*jet3mass_5) + (jet3px_5*jet3px_5) + (jet3py_5*jet3py_5) + (jet3pz_5*jet3pz_5))),0.5)) + (2*pow((((jet2mass_5*jet2mass_5) + (jet2px_5*jet2px_5) + (jet2py_5*jet2py_5) + (jet2pz_5*jet2pz_5))*((jet3mass_5*jet3mass_5) + (jet3px_5*jet3px_5) + (jet3py_5*jet3py_5) + (jet3pz_5*jet3pz_5))),0.5)) - (2*((jet1px_5*jet2px_5) + (jet1py_5*jet2py_5) + (jet1pz_5*jet2pz_5))) - (2*((jet1px_5*jet3px_5) + (jet1py_5*jet3py_5) + (jet1pz_5*jet3pz_5))) - (2*((jet2px_5*jet3px_5) + (jet2py_5*jet3py_5) + (jet2pz_5*jet3pz_5)))),0.5);
		
		if (jet1pt_5 < 20 || jet2pt_5 < 20 || jet3pt_5 < 20) continue;
					
		hmtot_5->Fill(mtot_5[p]);
		
	}
	
	cout << "\nNumber of entries of the histogram 3: " << hmtot_3->GetEntries() << endl;
	cout << "Number of entries of the histogram 4: " << hmtot_4->GetEntries() << endl;
	cout << "Number of entries of the histogram 5: " << hmtot_5->GetEntries() << endl;
	cout << endl;
	
	hmtot_3->Scale(F3);
	hmtot_4->Scale(F4);
	hmtot_5->Scale(F5);
	
	hmtot_3->Scale(1/(hmtot_3->Integral()));
	hmtot_4->Scale(1/(hmtot_4->Integral()));
	hmtot_5->Scale(1/(hmtot_5->Integral()));
	
	hmtot_3->GetXaxis()->SetRangeUser(0,240);
	hmtot_3->GetXaxis()->SetNdivisions(511);
	hmtot_3->SetLineWidth(1);
	hmtot_3->SetMarkerStyle(20);
	hmtot_3->SetMarkerSize(0.8);
	hmtot_3->SetMarkerColor(kBlue+2);
	
	hmtot_4->GetXaxis()->SetRangeUser(0,240);
	hmtot_4->GetXaxis()->SetNdivisions(511);
	hmtot_4->SetLineWidth(1);
	hmtot_4->SetLineColor(kRed);
	hmtot_4->SetMarkerStyle(20);
	hmtot_4->SetMarkerSize(0.8);
	hmtot_4->SetMarkerColor(kRed);
	
	hmtot_5->GetXaxis()->SetRangeUser(0,240);
	hmtot_5->GetXaxis()->SetNdivisions(511);
	hmtot_5->SetLineWidth(1);
	hmtot_5->SetLineColor(kGreen+1);
	hmtot_5->SetMarkerStyle(20);
	hmtot_5->SetMarkerSize(0.8);
	hmtot_5->SetMarkerColor(kGreen+1);
	
	gStyle->SetOptStat("");
	
	c5->cd();
	
	hmtot_4->Draw();
	hmtot_3->Draw("SAME");
	hmtot_5->Draw("SAME");
	
	TLegend leg2(.5,.7,.9,.9); // (x1,y1,x2,y2)
	
	leg2.SetFillColor(0);
	leg2.AddEntry("hmtot_3","CT10nnlo #alpha_{s}=0.1100");
	leg2.AddEntry("hmtot_5","CT10nnlo #alpha_{s}=0.1180");
	leg2.AddEntry("hmtot_4","CT10nnlo #alpha_{s}=0.1300");
	leg2.DrawClone("SAME");
	
}

void eta12Norm() {
	
	TCanvas *c6 = new TCanvas("c6","eta12");
		
	cout << "\nScaling factors:\tF3 = " << F3 << "\n\t\t\tF4 = " << F4 << "\n\t\t\tF5 = " << F5 << endl;
	
	TFile *f3 = new TFile("jets_ct10_10gev_0110.root","READ");
	TFile *f4 = new TFile("jets_ct10_10gev_0130.root","READ");
	TFile *f5 = new TFile("jets_ct10_10gev_0118.root","READ");
	
	TTree *t3 = (TTree*)f3->Get("jets_ct10_10gev_0110");
	TTree *t4 = (TTree*)f4->Get("jets_ct10_10gev_0130");
	TTree *t5 = (TTree*)f5->Get("jets_ct10_10gev_0118");
	
	int N3 = t3->GetEntries();
	int N4 = t4->GetEntries();
	int N5 = t5->GetEntries();
	
	cout << "\nNumber of entries:\tN3 = " << N3 << "\n\t\t\tN4 = " << N4 << "\n\t\t\tN5 = " << N5 << endl;
	
	// variables to store branch content of tree t3
	double jet1pt_3;
	double jet2pt_3;
	double jet3pt_3;
	double jet1eta_3;
	double jet2eta_3;
	
	t3->SetBranchAddress("jet1pt",&jet1pt_3);
	t3->SetBranchAddress("jet2pt",&jet2pt_3);
	t3->SetBranchAddress("jet3pt",&jet3pt_3);
	t3->SetBranchAddress("jet1eta",&jet1eta_3);
	t3->SetBranchAddress("jet2eta",&jet2eta_3);
	
	// array to store eta12 values
	double eta12_3[619972] = {0.0};
	
	// variables to store branch content of tree t4
	double jet1pt_4;
	double jet2pt_4;
	double jet3pt_4;
	double jet1eta_4;
	double jet2eta_4;
	
	t4->SetBranchAddress("jet1pt",&jet1pt_4);
	t4->SetBranchAddress("jet2pt",&jet2pt_4);
	t4->SetBranchAddress("jet3pt",&jet3pt_4);
	t4->SetBranchAddress("jet1eta",&jet1eta_4);
	t4->SetBranchAddress("jet2eta",&jet2eta_4);
	
	// array to store eta12 values
	double eta12_4[622422] = {0.0};
	
	// variables to store branch content of tree t5
	double jet1pt_5;
	double jet2pt_5;
	double jet3pt_5;
	double jet1eta_5;
	double jet2eta_5;
	
	t5->SetBranchAddress("jet1pt",&jet1pt_5);
	t5->SetBranchAddress("jet2pt",&jet2pt_5);
	t5->SetBranchAddress("jet3pt",&jet3pt_5);
	t5->SetBranchAddress("jet1eta",&jet1eta_5);
	t5->SetBranchAddress("jet2eta",&jet2eta_5);
	
	// array to store eta12 values
	double eta12_5[621619] = {0.0};
	
	// histograms of distributions
	TH1D *heta12_3 = new TH1D("heta12_3","#Delta#eta_{12} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;#Delta#eta_{12};A.U.",25,-2.5,2.5);
	TH1D *heta12_4 = new TH1D("heta12_4","#Delta#eta_{12} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;#Delta#eta_{12};A.U.",25,-2.5,2.5);
	TH1D *heta12_5 = new TH1D("heta12_5","#Delta#eta_{12} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;#Delta#eta_{12};A.U.",25,-2.5,2.5);
	
	for (int n=0; n<N3; n++) {

		t3->GetEntry(n);
		
		eta12_3[n] = jet1eta_3 - jet2eta_3;
		
		if (jet1pt_3 < 20 || jet2pt_3 < 20 || jet3pt_3 < 20) continue;
					
		heta12_3->Fill(eta12_3[n]);
		
	}
	
	for (int m=0; m<N4; m++) {

		t4->GetEntry(m);
		
		eta12_4[m] = jet1eta_4 - jet2eta_4;
		
		if (jet1pt_4 < 20 || jet2pt_4 < 20 || jet3pt_4 < 20) continue;
					
		heta12_4->Fill(eta12_4[m]);
		
	}
	
	for (int p=0; p<N4; p++) {

		t5->GetEntry(p);
		
		eta12_5[p] = jet1eta_5 - jet2eta_5;
		
		if (jet1pt_5 < 20 || jet2pt_5 < 20 || jet3pt_5 < 20) continue;
					
		heta12_5->Fill(eta12_5[p]);
		
	}
	
	cout << "\nNumber of entries of the histogram 3: " << heta12_3->GetEntries() << endl;
	cout << "Number of entries of the histogram 4: " << heta12_4->GetEntries() << endl;
	cout << "Number of entries of the histogram 5: " << heta12_5->GetEntries() << endl;
	cout << endl;
	
	heta12_3->Scale(F3);
	heta12_4->Scale(F4);
	heta12_5->Scale(F5);
	
	heta12_3->Scale(1/(heta12_3->Integral()));
	heta12_4->Scale(1/(heta12_4->Integral()));
	heta12_5->Scale(1/(heta12_5->Integral()));
	
	heta12_3->GetXaxis()->SetRangeUser(-2.5,2.5);
	heta12_3->GetXaxis()->SetNdivisions(511);
	heta12_3->SetLineWidth(1);
	heta12_3->SetMarkerStyle(20);
	heta12_3->SetMarkerSize(0.8);
	heta12_3->SetMarkerColor(kBlue+2);
	
	heta12_4->GetXaxis()->SetRangeUser(-2.5,2.5);
	heta12_4->GetXaxis()->SetNdivisions(511);
	heta12_4->SetLineWidth(1);
	heta12_4->SetLineColor(kRed);
	heta12_4->SetMarkerStyle(20);
	heta12_4->SetMarkerSize(0.8);
	heta12_4->SetMarkerColor(kRed);
	
	heta12_5->GetXaxis()->SetRangeUser(-2.5,2.5);
	heta12_5->GetXaxis()->SetNdivisions(511);
	heta12_5->SetLineWidth(1);
	heta12_5->SetLineColor(kGreen+1);
	heta12_5->SetMarkerStyle(20);
	heta12_5->SetMarkerSize(0.8);
	heta12_5->SetMarkerColor(kGreen+1);
	
	gStyle->SetOptStat("");
	
	c6->cd();
	
	heta12_4->Draw();
	heta12_3->Draw("SAME");
	heta12_5->Draw("SAME");
	
	TLegend leg2(.5,.7,.9,.9); // (x1,y1,x2,y2)
	
	leg2.SetFillColor(0);
	leg2.AddEntry("heta12_3","CT10nnlo #alpha_{s}=0.1100");
	leg2.AddEntry("heta12_5","CT10nnlo #alpha_{s}=0.1180");
	leg2.AddEntry("heta12_4","CT10nnlo #alpha_{s}=0.1300");
	leg2.DrawClone("SAME");
	
}

void phi12Norm() {
	
	TCanvas *c7 = new TCanvas("c7","phi12");
	
	cout << "\nScaling factors:\tF3 = " << F3 << "\n\t\t\tF4 = " << F4 << "\n\t\t\tF5 = " << F5 << endl;
	
	TFile *f3 = new TFile("jets_ct10_10gev_0110.root","READ");
	TFile *f4 = new TFile("jets_ct10_10gev_0130.root","READ");
	TFile *f5 = new TFile("jets_ct10_10gev_0118.root","READ");
	
	TTree *t3 = (TTree*)f3->Get("jets_ct10_10gev_0110");
	TTree *t4 = (TTree*)f4->Get("jets_ct10_10gev_0130");
	TTree *t5 = (TTree*)f5->Get("jets_ct10_10gev_0118");
	
	int N3 = t3->GetEntries();
	int N4 = t4->GetEntries();
	int N5 = t5->GetEntries();
	
	cout << "\nNumber of entries:\tN3 = " << N3 << "\n\t\t\tN4 = " << N4 << "\n\t\t\tN5 = " << N5 << endl;
	
	// variables to store branch content of tree t3
	double jet1pt_3;
	double jet2pt_3;
	double jet3pt_3;
	double jet1px_3;
	double jet1py_3;
	double jet2px_3;
	double jet2py_3;
	
	t3->SetBranchAddress("jet1pt",&jet1pt_3);
	t3->SetBranchAddress("jet2pt",&jet2pt_3);
	t3->SetBranchAddress("jet3pt",&jet3pt_3);
	t3->SetBranchAddress("jet1px",&jet1px_3);
	t3->SetBranchAddress("jet1py",&jet1py_3);
	t3->SetBranchAddress("jet2px",&jet2px_3);
	t3->SetBranchAddress("jet2py",&jet2py_3);
	
	// array to store phi12 values
	double phi12_3[619972] = {0.0};
	
	// variables to store branch content of tree t4
	double jet1pt_4;
	double jet2pt_4;
	double jet3pt_4;
	double jet1px_4;
	double jet1py_4;
	double jet2px_4;
	double jet2py_4;
	
	t4->SetBranchAddress("jet1pt",&jet1pt_4);
	t4->SetBranchAddress("jet2pt",&jet2pt_4);
	t4->SetBranchAddress("jet3pt",&jet3pt_4);
	t4->SetBranchAddress("jet1px",&jet1px_4);
	t4->SetBranchAddress("jet1py",&jet1py_4);
	t4->SetBranchAddress("jet2px",&jet2px_4);
	t4->SetBranchAddress("jet2py",&jet2py_4);
	
	// array to store phi12 values
	double phi12_4[622422] = {0.0};
	
	// variables to store branch content of tree t5
	double jet1pt_5;
	double jet2pt_5;
	double jet3pt_5;
	double jet1px_5;
	double jet1py_5;
	double jet2px_5;
	double jet2py_5;
	
	t5->SetBranchAddress("jet1pt",&jet1pt_5);
	t5->SetBranchAddress("jet2pt",&jet2pt_5);
	t5->SetBranchAddress("jet3pt",&jet3pt_5);
	t5->SetBranchAddress("jet1px",&jet1px_5);
	t5->SetBranchAddress("jet1py",&jet1py_5);
	t5->SetBranchAddress("jet2px",&jet2px_5);
	t5->SetBranchAddress("jet2py",&jet2py_5);
	
	// array to store phi12 values
	double phi12_5[621619] = {0.0};
	
	// histograms of distributions
	TH1D *hphi12_3 = new TH1D("hphi12_3","#Delta#phi_{12} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;#Delta#phi_{12};A.U.",25,0,3.5);
	TH1D *hphi12_4 = new TH1D("hphi12_4","#Delta#phi_{12} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;#Delta#phi_{12};A.U.",25,0,3.5);
	TH1D *hphi12_5 = new TH1D("hphi12_5","#Delta#phi_{12} - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;#Delta#phi_{12};A.U.",25,0,3.5);
	
	for (int n=0; n<N3; n++) {

		t3->GetEntry(n);
		
		phi12_3[n] = acos(((jet1px_3*jet2px_3) + (jet1py_3*jet2py_3))/(pow((((jet1px_3*jet1px_3) + (jet1py_3*jet1py_3))*((jet2px_3*jet2px_3) + (jet2py_3*jet2py_3))),0.5)));
		
		if (jet1pt_3 < 20 || jet2pt_3 < 20 || jet3pt_3 < 20) continue;
					
		hphi12_3->Fill(phi12_3[n]);
		
	}
	
	for (int m=0; m<N4; m++) {

		t4->GetEntry(m);
		
		phi12_4[m] = acos(((jet1px_4*jet2px_4) + (jet1py_4*jet2py_4))/(pow((((jet1px_4*jet1px_4) + (jet1py_4*jet1py_4))*((jet2px_4*jet2px_4) + (jet2py_4*jet2py_4))),0.5)));
		
		if (jet1pt_4 < 20 || jet2pt_4 < 20 || jet3pt_4 < 20) continue;
					
		hphi12_4->Fill(phi12_4[m]);
		
	}
	
	for (int p=0; p<N5; p++) {

		t5->GetEntry(p);
		
		phi12_5[p] = acos(((jet1px_5*jet2px_5) + (jet1py_5*jet2py_5))/(pow((((jet1px_5*jet1px_5) + (jet1py_5*jet1py_5))*((jet2px_5*jet2px_5) + (jet2py_5*jet2py_5))),0.5)));
		
		if (jet1pt_5 < 20 || jet2pt_5 < 20 || jet3pt_5 < 20) continue;
					
		hphi12_5->Fill(phi12_5[p]);
		
	}
	
	cout << "\nNumber of entries of the histogram 3: " << hphi12_3->GetEntries() << endl;
	cout << "Number of entries of the histogram 4: " << hphi12_4->GetEntries() << endl;
	cout << "Number of entries of the histogram 5: " << hphi12_5->GetEntries() << endl;
	cout << endl;
	
	hphi12_3->Scale(F3);
	hphi12_4->Scale(F4);
	hphi12_5->Scale(F5);
	
	hphi12_3->Scale(1/(hphi12_3->Integral()));
	hphi12_4->Scale(1/(hphi12_4->Integral()));
	hphi12_5->Scale(1/(hphi12_5->Integral()));
	
	hphi12_3->GetXaxis()->SetRangeUser(0,3.5);
	hphi12_3->GetXaxis()->SetNdivisions(511);
	hphi12_3->SetLineWidth(1);
	hphi12_3->SetMarkerStyle(20);
	hphi12_3->SetMarkerSize(0.8);
	hphi12_3->SetMarkerColor(kBlue+2);
	
	hphi12_4->GetXaxis()->SetRangeUser(0,3.5);
	hphi12_4->GetXaxis()->SetNdivisions(511);
	hphi12_4->SetLineWidth(1);
	hphi12_4->SetLineColor(kRed);
	hphi12_4->SetMarkerStyle(20);
	hphi12_4->SetMarkerSize(0.8);
	hphi12_4->SetMarkerColor(kRed);
	
	hphi12_5->GetXaxis()->SetRangeUser(0,3.5);
	hphi12_5->GetXaxis()->SetNdivisions(511);
	hphi12_5->SetLineWidth(1);
	hphi12_5->SetLineColor(kGreen+1);
	hphi12_5->SetMarkerStyle(20);
	hphi12_5->SetMarkerSize(0.8);
	hphi12_5->SetMarkerColor(kGreen+1);
	
	gStyle->SetOptStat("");
	
	c7->cd();
	
	hphi12_4->Draw();
	hphi12_3->Draw("SAME");
	hphi12_5->Draw("SAME");
	
	TLegend leg2(.1,.7,.5,.9); // (x1,y1,x2,y2)
	
	leg2.SetFillColor(0);
	leg2.AddEntry("hphi12_3","CT10nnlo #alpha_{s}=0.1100");
	leg2.AddEntry("hphi12_5","CT10nnlo #alpha_{s}=0.1180");
	leg2.AddEntry("hphi12_4","CT10nnlo #alpha_{s}=0.1300");
	leg2.DrawClone("SAME");
	
}
