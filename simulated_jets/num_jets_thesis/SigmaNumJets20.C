// Plot sezione d'urto in funzione del numero di jet: sezione d'urto calcolata come sigma*A, dove sigma è la sezione d'urto di MadGraph e A è l'accettanza (cioè la frazione di eventi che passano i tagli, Nsel/Ngen); A calcolata richiedendo N>=1, N>=2, N>=3 etc. e poi plottata

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

double Ngen = 1e6;		   // number of generated events

double sigma110 = 4.26194328e9;   // cross section [pb] for as = 0.110
double sigma130 = 5.3888668e9;    // cross section [pb] for as = 0.130
double sigma118 = 4.6643252e9;    // cross section [pb] for as = 0.118

double err_sigma110 = 9.891907e5; // error of cross section [pb] for as = 0.110
double err_sigma130 = 1.379683e6; // error of cross section [pb] for as = 0.130
double err_sigma118 = 1.063324e6; // error of cross section [pb] for as = 0.118

void SigmaNumJets10_0110() {
	
	TFile *f = new TFile("pp_jj_ct10_10gev_0110.root");
	
	TTree *t = (TTree*)f->Get("hepmc3_tree");
	
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
	
	// array to store number of jets per event in loop
	double Njets[999999] = {0.0};
	
	// array to store number of jets per event
	double Nsel[10] = {0.0};
	
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
		
		Nsel[i-1] = hnj->GetBinContent(i);
		
		cout << "Number of events with n >= " << i-1 << " jets: " << Nsel[i-1] << endl;
		
	}
	
	double A0 = Nsel[0]/Ngen;
	double A1 = Nsel[1]/Ngen;
	double A2 = Nsel[2]/Ngen;
	double A3 = Nsel[3]/Ngen;
	double A4 = Nsel[4]/Ngen;
	double A5 = Nsel[5]/Ngen;
	double A6 = Nsel[6]/Ngen;
	
	cout << "\nAcceptance for NumJets>=0: " << A0 << "\nAcceptance for NumJets>=1: " << A1 << "\nAcceptance for NumJets>=2: " << A2 << "\nAcceptance for NumJets>=3: " << A3 << "\nAcceptance for NumJets>=4: " << A4 << "\nAcceptance for NumJets>=5: " << A5 << "\nAcceptance for NumJets>=6: " << A6 << endl;
	
	double s0 = sigma110*A0;
	double s1 = sigma110*A1;
	double s2 = sigma110*A2;
	double s3 = sigma110*A3;
	double s4 = sigma110*A4;
	double s5 = sigma110*A5;
	double s6 = sigma110*A6;
	
	cout << "\nCross section for NumJets>=0: " << s0 << " pb\nCross section for NumJets>=1: " << s1 << " pb\nCross section for NumJets>=2: " << s2 << " pb\nCross section for NumJets>=3: " << s3 << " pb\nCross section for NumJets>=4: " << s4 << " pb\nCross section for NumJets>=5: " << s5 << " pb\nCross section for NumJets>=6: " << s6 << " pb\n" << endl;
	
	double err_A0 = (sqrt(Nsel[0]))/Ngen;
	double err_A1 = (sqrt(Nsel[1]))/Ngen;
	double err_A2 = (sqrt(Nsel[2]))/Ngen;
	double err_A3 = (sqrt(Nsel[3]))/Ngen;
	double err_A4 = (sqrt(Nsel[4]))/Ngen;
	double err_A5 = (sqrt(Nsel[5]))/Ngen;
	double err_A6 = (sqrt(Nsel[6]))/Ngen;
	
	double err_s0 = sqrt((sigma110*sigma110*err_A0*err_A0)+(A0*A0*err_sigma110*err_sigma110));
	double err_s1 = sqrt((sigma110*sigma110*err_A1*err_A1)+(A1*A1*err_sigma110*err_sigma110));
	double err_s2 = sqrt((sigma110*sigma110*err_A2*err_A2)+(A2*A2*err_sigma110*err_sigma110));
	double err_s3 = sqrt((sigma110*sigma110*err_A3*err_A3)+(A3*A3*err_sigma110*err_sigma110));
	double err_s4 = sqrt((sigma110*sigma110*err_A4*err_A4)+(A4*A4*err_sigma110*err_sigma110));
	double err_s5 = sqrt((sigma110*sigma110*err_A5*err_A5)+(A5*A5*err_sigma110*err_sigma110));
	double err_s6 = sqrt((sigma110*sigma110*err_A6*err_A6)+(A6*A6*err_sigma110*err_sigma110));
	
	ofstream outfile ("SigmaNumJets10_0110.txt");
	
	outfile << 0 << "\t" << s0 << "\t" << err_s0 << "\n" << 1 << "\t" << s1 << "\t" << err_s1 << "\n" << 2 << "\t" << s2 << "\t" << err_s2 << "\n" << 3 << "\t" << s3 << "\t" << err_s3 << "\n" << 4 << "\t" << s4 << "\t" << err_s4 << "\n" << 5 << "\t" << s5 << "\t" << err_s5 << "\n" << 6 << "\t" << s6 << "\t" << err_s6 << endl;
	
}

void SigmaNumJets10_0130() {
	
	TFile *f = new TFile("pp_jj_ct10_10gev_0130.root");
	
	TTree *t = (TTree*)f->Get("hepmc3_tree");
	
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
	
	// array to store number of jets per event in loop
	double Njets[999999] = {0.0};
	
	// array to store number of jets per event
	double Nsel[10] = {0.0};
	
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
		
		Nsel[i-1] = hnj->GetBinContent(i);
		
		cout << "Number of events with n >= " << i-1 << " jets: " << Nsel[i-1] << endl;
		
	}
	
	double A0 = Nsel[0]/Ngen;
	double A1 = Nsel[1]/Ngen;
	double A2 = Nsel[2]/Ngen;
	double A3 = Nsel[3]/Ngen;
	double A4 = Nsel[4]/Ngen;
	double A5 = Nsel[5]/Ngen;
	double A6 = Nsel[6]/Ngen;
	
	cout << "\nAcceptance for NumJets>=0: " << A0 << "\nAcceptance for NumJets>=1: " << A1 << "\nAcceptance for NumJets>=2: " << A2 << "\nAcceptance for NumJets>=3: " << A3 << "\nAcceptance for NumJets>=4: " << A4 << "\nAcceptance for NumJets>=5: " << A5 << "\nAcceptance for NumJets>=6: " << A6 << endl;
	
	double s0 = sigma130*A0;
	double s1 = sigma130*A1;
	double s2 = sigma130*A2;
	double s3 = sigma130*A3;
	double s4 = sigma130*A4;
	double s5 = sigma130*A5;
	double s6 = sigma130*A6;
	
	cout << "\nCross section for NumJets>=0: " << s0 << " pb\nCross section for NumJets>=1: " << s1 << " pb\nCross section for NumJets>=2: " << s2 << " pb\nCross section for NumJets>=3: " << s3 << " pb\nCross section for NumJets>=4: " << s4 << " pb\nCross section for NumJets>=5: " << s5 << " pb\nCross section for NumJets>=6: " << s6 << " pb\n" << endl;
	
	double err_A0 = (sqrt(Nsel[0]))/Ngen;
	double err_A1 = (sqrt(Nsel[1]))/Ngen;
	double err_A2 = (sqrt(Nsel[2]))/Ngen;
	double err_A3 = (sqrt(Nsel[3]))/Ngen;
	double err_A4 = (sqrt(Nsel[4]))/Ngen;
	double err_A5 = (sqrt(Nsel[5]))/Ngen;
	double err_A6 = (sqrt(Nsel[6]))/Ngen;
	
	double err_s0 = sqrt((sigma130*sigma130*err_A0*err_A0)+(A0*A0*err_sigma130*err_sigma130));
	double err_s1 = sqrt((sigma130*sigma130*err_A1*err_A1)+(A1*A1*err_sigma130*err_sigma130));
	double err_s2 = sqrt((sigma130*sigma130*err_A2*err_A2)+(A2*A2*err_sigma130*err_sigma130));
	double err_s3 = sqrt((sigma130*sigma130*err_A3*err_A3)+(A3*A3*err_sigma130*err_sigma130));
	double err_s4 = sqrt((sigma130*sigma130*err_A4*err_A4)+(A4*A4*err_sigma130*err_sigma130));
	double err_s5 = sqrt((sigma130*sigma130*err_A5*err_A5)+(A5*A5*err_sigma130*err_sigma130));
	double err_s6 = sqrt((sigma130*sigma130*err_A6*err_A6)+(A6*A6*err_sigma130*err_sigma130));
	
	ofstream outfile ("SigmaNumJets10_0130.txt");
	
	outfile << 0 << "\t" << s0 << "\t" << err_s0 << "\n" << 1 << "\t" << s1 << "\t" << err_s1 << "\n" << 2 << "\t" << s2 << "\t" << err_s2 << "\n" << 3 << "\t" << s3 << "\t" << err_s3 << "\n" << 4 << "\t" << s4 << "\t" << err_s4 << "\n" << 5 << "\t" << s5 << "\t" << err_s5 << "\n" << 6 << "\t" << s6 << "\t" << err_s6 << endl;
	
}

void SigmaNumJets10_0118() {
	
	TFile *f = new TFile("pp_jj_ct10_10gev_0118.root");
	
	TTree *t = (TTree*)f->Get("hepmc3_tree");
	
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
	
	// array to store number of jets per event in loop
	double Njets[999999] = {0.0};
	
	// array to store number of jets per event
	double Nsel[10] = {0.0};
	
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
		
		Nsel[i-1] = hnj->GetBinContent(i);
		
		cout << "Number of events with n >= " << i-1 << " jets: " << Nsel[i-1] << endl;
		
	}
	
	double A0 = Nsel[0]/Ngen;
	double A1 = Nsel[1]/Ngen;
	double A2 = Nsel[2]/Ngen;
	double A3 = Nsel[3]/Ngen;
	double A4 = Nsel[4]/Ngen;
	double A5 = Nsel[5]/Ngen;
	double A6 = Nsel[6]/Ngen;
	
	cout << "\nAcceptance for NumJets>=0: " << A0 << "\nAcceptance for NumJets>=1: " << A1 << "\nAcceptance for NumJets>=2: " << A2 << "\nAcceptance for NumJets>=3: " << A3 << "\nAcceptance for NumJets>=4: " << A4 << "\nAcceptance for NumJets>=5: " << A5 << "\nAcceptance for NumJets>=6: " << A6 << endl;
	
	double s0 = sigma118*A0;
	double s1 = sigma118*A1;
	double s2 = sigma118*A2;
	double s3 = sigma118*A3;
	double s4 = sigma118*A4;
	double s5 = sigma118*A5;
	double s6 = sigma118*A6;
	
	cout << "\nCross section for NumJets>=0: " << s0 << " pb\nCross section for NumJets>=1: " << s1 << " pb\nCross section for NumJets>=2: " << s2 << " pb\nCross section for NumJets>=3: " << s3 << " pb\nCross section for NumJets>=4: " << s4 << " pb\nCross section for NumJets>=5: " << s5 << " pb\nCross section for NumJets>=6: " << s6 << " pb\n" << endl;
	
	double err_A0 = (sqrt(Nsel[0]))/Ngen;
	double err_A1 = (sqrt(Nsel[1]))/Ngen;
	double err_A2 = (sqrt(Nsel[2]))/Ngen;
	double err_A3 = (sqrt(Nsel[3]))/Ngen;
	double err_A4 = (sqrt(Nsel[4]))/Ngen;
	double err_A5 = (sqrt(Nsel[5]))/Ngen;
	double err_A6 = (sqrt(Nsel[6]))/Ngen;
	
	double err_s0 = sqrt((sigma118*sigma118*err_A0*err_A0)+(A0*A0*err_sigma118*err_sigma118));
	double err_s1 = sqrt((sigma118*sigma118*err_A1*err_A1)+(A1*A1*err_sigma118*err_sigma118));
	double err_s2 = sqrt((sigma118*sigma118*err_A2*err_A2)+(A2*A2*err_sigma118*err_sigma118));
	double err_s3 = sqrt((sigma118*sigma118*err_A3*err_A3)+(A3*A3*err_sigma118*err_sigma118));
	double err_s4 = sqrt((sigma118*sigma118*err_A4*err_A4)+(A4*A4*err_sigma118*err_sigma118));
	double err_s5 = sqrt((sigma118*sigma118*err_A5*err_A5)+(A5*A5*err_sigma118*err_sigma118));
	double err_s6 = sqrt((sigma118*sigma118*err_A6*err_A6)+(A6*A6*err_sigma118*err_sigma118));
	
	ofstream outfile ("SigmaNumJets10_0118.txt");
	
	outfile << 0 << "\t" << s0 << "\t" << err_s0 << "\n" << 1 << "\t" << s1 << "\t" << err_s1 << "\n" << 2 << "\t" << s2 << "\t" << err_s2 << "\n" << 3 << "\t" << s3 << "\t" << err_s3 << "\n" << 4 << "\t" << s4 << "\t" << err_s4 << "\n" << 5 << "\t" << s5 << "\t" << err_s5 << "\n" << 6 << "\t" << s6 << "\t" << err_s6 << endl;
	
}

void SigmaNumJetsPlot() {
	
	TGraphErrors *g3 = new TGraphErrors("SigmaNumJets10_0110.txt","%lg %lg %lg");
	TGraphErrors *g4 = new TGraphErrors("SigmaNumJets10_0130.txt","%lg %lg %lg");
	TGraphErrors *g5 = new TGraphErrors("SigmaNumJets10_0118.txt","%lg %lg %lg");
	
	TCanvas *c = new TCanvas("c","SigmaNumJets");
	
	g3->SetName("g3");
	g4->SetName("g4");
	g5->SetName("g5");
	
	g3->SetLineColor(kBlue+2);
	g3->SetMarkerColor(kBlue+2);
	g3->SetMarkerStyle(20);
	g3->SetMarkerSize(1.2);
	g3->SetFillColor(kBlue+2);
	g3->SetFillStyle(3003);
	
	g4->SetLineColor(kRed);
	g4->SetMarkerColor(kRed);
	g4->SetMarkerStyle(20);
	g4->SetMarkerSize(1.2);
	g4->SetFillColor(kRed);
	g4->SetFillStyle(3003);
	
	g5->SetLineColor(kGreen+2);
	g5->SetMarkerColor(kGreen+2);
	g5->SetMarkerStyle(20);
	g5->SetMarkerSize(1.2);
	g5->SetFillColor(kGreen+2);
	g5->SetFillStyle(3003);
	
	c->cd();
	
	TMultiGraph *mg = new TMultiGraph();
	
	mg->SetTitle("Cross section weighted by acceptance - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;NumJetsPerEvent;#sigma_{nj} [pb]");
	mg->Add(g3,"P");
	mg->Add(g4,"P");
	mg->Add(g5,"P");
	
	mg->GetXaxis()->SetNdivisions(10);
	
	mg->Draw("A");
	
	TLegend leg(.5,.7,.9,.9); // (x1,y1,x2,y2)
	
	leg.AddEntry("g3","CT10nnlo #alpha_{s}=0.1100","EP");
	leg.AddEntry("g5","CT10nnlo #alpha_{s}=0.1180","EP");
	leg.AddEntry("g4","CT10nnlo #alpha_{s}=0.1300","EP");
	leg.DrawClone("SAME");
	
}
