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
#include <vector>

using namespace std;

double L = 1.6;		// integrated luminosity of Run 2 [fb^-1]
double Ngen = 1e6;		// number of generated events

double sigma3 = 4.26194328e12; // cross section [fb] for pt > 10 GeV, as = 0.110
double sigma4 = 5.3888668e12;	// cross section [fb] for pt > 10 GeV, as = 0.130
double sigma5 = 4.6643252e12;	// cross section [fb] for pt > 10 GeV, as = 0.118

double err_sigma3 = 9.891907e8; // error of cross section [fb] for as = 0.110
double err_sigma4 = 1.379683e9; // error of cross section [fb] for as = 0.130
double err_sigma5 = 1.063324e9; // error of cross section [fb] for as = 0.118

double F3 = ((L*sigma3)/Ngen); // scaling factor for pt > 10 GeV, as = 0.110
double F4 = ((L*sigma4)/Ngen); // scaling factor for pt > 10 GeV, as = 0.130
double F5 = ((L*sigma5)/Ngen); // scaling factor for pt > 10 GeV, as = 0.118

double err_F3 = L*1e-8*sqrt((10000*err_sigma3*err_sigma3)+(4*sigma3*sigma3));
double err_F4 = L*1e-8*sqrt((10000*err_sigma4*err_sigma4)+(4*sigma4*sigma4));
double err_F5 = L*1e-8*sqrt((10000*err_sigma5*err_sigma5)+(4*sigma5*sigma5));

void ErrF() {
	
	cout << "Scale factors:\tF3 = " << F3 << " ± " << err_F3 << "\n\t\tF4 = " << F4 << " ± " << err_F4 << "\n\t\tF5 = " << F5 << " ± " << err_F5 << endl;
	
}

void NumJets10_110() {
	
	TFile *f = new TFile("pp_jj_ct10_10gev_0110.root");
	
	TTree *t = (TTree*)f->Get("hepmc3_tree");
	
	TCanvas *c1 = new TCanvas("c1","NumJetsPerEvent0110");
	
	TH1D *hnj = new TH1D("hnj","Inclusive number of jets per event - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;NumJetsPerEvent;A.U.",10,0,10);
	
	int N = t->GetEntries();
	
	cout << "\nNumber of entries = " << N << endl;
	
	// arrays to store momentum components of jets
	double px[50] = {0.0};
	double py[50] = {0.0};
	double pz[50] = {0.0};
	
	TBranch *b_px;
	TBranch *b_py;
	TBranch *b_pz;
	
	t->SetBranchAddress("particles.momentum.m_v1",px,&b_px);
	t->SetBranchAddress("particles.momentum.m_v2",py,&b_py);
	t->SetBranchAddress("particles.momentum.m_v3",pz,&b_pz);
	
	// array to store number of jets per event
	double Njets[999999] = {0.0};
	
	// loop over events
	for (int e=0; e<N+1; e++) {
		
		cout << "\n++++++++++ EVENT NUMBER " << e << " ++++++++++" << endl;
		
		// reset arrays for each event to avoid data carryover
		memset(px,0.0,sizeof(px));
		memset(py,0.0,sizeof(py));
		memset(pz,0.0,sizeof(pz));
		
		t->GetEntry(e);
		
		// counter for number of jets in each event
		int nj = 0;
		
		// loop over jets per event
		for (int j=0; j<50; j++) {
			
			double eta = -log(tan(acos(pz[j]/sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]))/2));
			
			double pt = sqrt((px[j]*px[j])+(py[j]*py[j]));
			
			if (isnan(eta)) continue;
			
			// reject jets outside LHCb acceptance 2.2 < eta < 4.2
			if (eta < 2.2 || eta > 4.2) continue;
			
			// select only jets with pt >= 20 GeV
			if (pt < 20) continue;
			
			nj++;
			
			Njets[e] = nj;
			
			cout << "Jet inside LHCb acceptance number " << nj << endl;
			cout << "eta = " << eta << endl;
			cout << "pT = " << pt << endl;
			
		}
		
		hnj->Fill(Njets[e]);
		
	}
	
	int nBins = hnj->GetNbinsX();
	
	cout << endl;
	
	for (int i=1; i<=nBins; i++) {
		
		double cont = hnj->GetBinContent(i);
		cout << "Number of events with n >= " << i-1 << " jets: " << cont << endl;
		
	}
	
	cout << "\nScaling factors:\tF3 = " << F3 << "\n\t\t\tF4 = " << F4 << "\n\t\t\tF5 = " << F5 << "\n" << endl;
	
	hnj->Sumw2(); // to correctly compute error bars when scaling
	
	hnj->Scale(F3);
	
	hnj->Scale(1/(hnj->Integral()));
	
	hnj->GetXaxis()->SetNdivisions(10);
//	hnj->GetXaxis()->SetRangeUser(1,10);
	gStyle->SetOptStat("emr");
	
	c1->cd();
	
	hnj->Draw();
	
	TFile *outfile = new TFile("NumJetsSameArea10_0110.root","RECREATE");
	
	hnj->Write();
	
//	cout << endl;
	
//	int nBins = hnj->GetNbinsX();
	
//	double NF[10] = {0.0};
//	double err_NF[10] = {0.0};
	
//	for (int i=1; i<=nBins; i++) {
		
//		err_NF[i-1] = sqrt((F3*F3*(hnj->GetBinContent(i)))+((hnj->GetBinContent(i))*(hnj->GetBinContent(i))*err_F3*err_F3));
		
//	}
	
//	hnj->Scale(F3);
	
//	ofstream outfile ("NumJetsArea0110.txt");
	
//	for (int i=1; i<=nBins; i++) {
		
//		NF[i-1] = hnj->GetBinContent(i);
		
//		outfile << i-1 << "\t" << NF[i-1] << "\t" << err_NF[i-1] << endl;
		
//	}
	
}

void NumJets10_130() {
	
	TFile *f = new TFile("pp_jj_ct10_10gev_0130.root");
	
	TTree *t = (TTree*)f->Get("hepmc3_tree");
	
	TCanvas *c2 = new TCanvas("c2","NumJetsPerEvent0130");
	
	TH1D *hnj = new TH1D("hnj","Inclusive number of jets per event - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;NumJetsPerEvent;A.U.",10,0,10);
	
	int N = t->GetEntries();
	
	cout << "\nNumber of entries = " << N << endl;
	
	// arrays to store momentum components of jets
	double px[50] = {0.0};
	double py[50] = {0.0};
	double pz[50] = {0.0};
	
	TBranch *b_px;
	TBranch *b_py;
	TBranch *b_pz;
	
	t->SetBranchAddress("particles.momentum.m_v1",px,&b_px);
	t->SetBranchAddress("particles.momentum.m_v2",py,&b_py);
	t->SetBranchAddress("particles.momentum.m_v3",pz,&b_pz);
	
	// array to store number of jets per event
	double Njets[999999] = {0.0};
	
	// loop over events
	for (int e=0; e<N+1; e++) {
		
		cout << "\n++++++++++ EVENT NUMBER " << e << " ++++++++++" << endl;
		
		// reset arrays for each event to avoid data carryover
		memset(px,0.0,sizeof(px));
		memset(py,0.0,sizeof(py));
		memset(pz,0.0,sizeof(pz));
		
		t->GetEntry(e);
		
		// counter for number of jets in each event
		int nj = 0;
		
		// loop over jets per event
		for (int j=0; j<50; j++) {
			
			double eta = -log(tan(acos(pz[j]/sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]))/2));
			
			double pt = sqrt((px[j]*px[j])+(py[j]*py[j]));
			
			if (isnan(eta)) continue;
			
			// reject jets outside LHCb acceptance 2.2 < eta < 4.2
			if (eta < 2.2 || eta > 4.2) continue;
			
			// select only jets with pt >= 20 GeV
			if (pt < 20) continue;
			
			nj++;
			
			Njets[e] = nj;
			
			cout << "Jet inside LHCb acceptance number " << nj << endl;
			cout << "eta = " << eta << endl;
			cout << "pT = " << pt << endl;
			
		}
		
		hnj->Fill(Njets[e]);
		
	}
	
	int nBins = hnj->GetNbinsX();
	
	cout << endl;
	
	for (int i=1; i<=nBins; i++) {
		
		double cont = hnj->GetBinContent(i);
		cout << "Number of events with n >= " << i-1 << " jets: " << cont << endl;
		
	}
	
	cout << "\nScaling factors:\tF3 = " << F3 << "\n\t\t\tF4 = " << F4 << "\n\t\t\tF5 = " << F5 << "\n" << endl;
	
	hnj->Sumw2(); // to correctly compute error bars when scaling
	
	hnj->Scale(F4);
	
	hnj->Scale(1/(hnj->Integral()));
	
	hnj->GetXaxis()->SetNdivisions(10);
//	hnj->GetXaxis()->SetRangeUser(1,10);
	gStyle->SetOptStat("emr");
	
	c2->cd();
	
	hnj->Draw();
	
	TFile *outfile = new TFile("NumJetsSameArea10_0130.root","RECREATE");
	
	hnj->Write();
	
}

void NumJets10_118() {
	
	TFile *f = new TFile("pp_jj_ct10_10gev_0118.root");
	
	TTree *t = (TTree*)f->Get("hepmc3_tree");
	
	TCanvas *c3 = new TCanvas("c3","NumJetsPerEvent0118");
	
	TH1D *hnj = new TH1D("hnj","Inclusive number of jets per event - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;NumJetsPerEvent;A.U.",10,0,10);
	
	int N = t->GetEntries();
	
	cout << "\nNumber of entries = " << N << endl;
	
	// arrays to store momentum components of jets
	double px[50] = {0.0};
	double py[50] = {0.0};
	double pz[50] = {0.0};
	
	TBranch *b_px;
	TBranch *b_py;
	TBranch *b_pz;
	
	t->SetBranchAddress("particles.momentum.m_v1",px,&b_px);
	t->SetBranchAddress("particles.momentum.m_v2",py,&b_py);
	t->SetBranchAddress("particles.momentum.m_v3",pz,&b_pz);
	
	// array to store number of jets per event
	double Njets[999999] = {0.0};
	
	// loop over events
	for (int e=0; e<N+1; e++) {
		
		cout << "\n++++++++++ EVENT NUMBER " << e << " ++++++++++" << endl;
		
		// reset arrays for each event to avoid data carryover
		memset(px,0.0,sizeof(px));
		memset(py,0.0,sizeof(py));
		memset(pz,0.0,sizeof(pz));
		
		t->GetEntry(e);
		
		// counter for number of jets in each event
		int nj = 0;
		
		// loop over jets per event
		for (int j=0; j<50; j++) {
			
			double eta = -log(tan(acos(pz[j]/sqrt(px[j]*px[j]+py[j]*py[j]+pz[j]*pz[j]))/2));
			
			double pt = sqrt((px[j]*px[j])+(py[j]*py[j]));
			
			if (isnan(eta)) continue;
			
			// reject jets outside LHCb acceptance 2.2 < eta < 4.2
			if (eta < 2.2 || eta > 4.2) continue;
			
			// select only jets with pt >= 20 GeV
			if (pt < 20) continue;
			
			nj++;
			
			Njets[e] = nj;
			
			cout << "Jet inside LHCb acceptance number " << nj << endl;
			cout << "eta = " << eta << endl;
			cout << "pT = " << pt << endl;
			
		}
		
		hnj->Fill(Njets[e]);
		
	}
	
	int nBins = hnj->GetNbinsX();
	
	cout << endl;
	
	for (int i=1; i<=nBins; i++) {
		
		double cont = hnj->GetBinContent(i);
		cout << "Number of events with n >= " << i-1 << " jets: " << cont << endl;
		
	}
	
	cout << "\nScaling factors:\tF3 = " << F3 << "\n\t\t\tF4 = " << F4 << "\n\t\t\tF5 = " << F5 << "\n" << endl;
	
	hnj->Sumw2(); // to correctly compute error bars when scaling
	
	hnj->Scale(F5);
	
	hnj->Scale(1/(hnj->Integral()));
	
	hnj->GetXaxis()->SetNdivisions(10);
//	hnj->GetXaxis()->SetRangeUser(1,10);
	gStyle->SetOptStat("emr");
	
	c3->cd();
	
	hnj->Draw();
	
	TFile *outfile = new TFile("NumJetsSameArea10_0118.root","RECREATE");
	
	hnj->Write();
	
}

void NumJetsPlots() {
	
	TFile *f3 = TFile::Open("NumJetsSameArea10_0110.root");
	TFile *f4 = TFile::Open("NumJetsSameArea10_0130.root");
	TFile *f5 = TFile::Open("NumJetsSameArea10_0118.root");
	
	TH1D *h3;
	f3->GetObject("hnj",h3);
	
	TH1D *h4;
	f4->GetObject("hnj",h4);
	
	TH1D *h5;
	f5->GetObject("hnj",h5);
	
	TCanvas *c = new TCanvas("c","NumJetsPerEventSameArea");
	
	h3->SetName("h3");
	h4->SetName("h4");
	h5->SetName("h5");
	
//	h3->GetXaxis()->SetRangeUser(1,8);
	h3->GetXaxis()->SetNdivisions(10);
	h3->GetXaxis()->CenterLabels();
	h3->SetLineColor(kBlue+2);
	h3->SetMarkerStyle(20);
	h3->SetMarkerSize(1);
	h3->SetMarkerColor(kBlue+2);
	
//	h4->GetXaxis()->SetRangeUser(1,8);
	h4->GetXaxis()->SetNdivisions(10);
	h4->GetXaxis()->CenterLabels();
	h4->SetLineColor(kRed);
	h4->SetMarkerStyle(20);
	h4->SetMarkerSize(1);
	h4->SetMarkerColor(kRed);
	
//	h5->GetXaxis()->SetRangeUser(1,8);
	h5->GetXaxis()->SetNdivisions(10);
	h5->GetXaxis()->CenterLabels();
	h5->SetLineColor(kGreen+2);
	h5->SetMarkerStyle(20);
	h5->SetMarkerSize(1);
	h5->SetMarkerColor(kGreen+2);

	gStyle->SetOptStat("");
	
	c->cd();
	
	h4->Draw();
	h3->Draw("SAME");
	h5->Draw("SAME");
	TLegend leg(.5,.7,.9,.9); // (x1,y1,x2,y2)
	
	leg.AddEntry("h3","CT10nnlo #alpha_{s}=0.1100","PL");
	leg.AddEntry("h5","CT10nnlo #alpha_{s}=0.1180","PL");
	leg.AddEntry("h4","CT10nnlo #alpha_{s}=0.1300","PL");
	leg.DrawClone("SAME");
	
}
