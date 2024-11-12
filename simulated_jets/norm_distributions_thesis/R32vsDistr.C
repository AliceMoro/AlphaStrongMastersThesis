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

using namespace std;

void R32vsQ() {
	
	TFile *f = new TFile("jets_ct10_10gev_0130.root","READ");
	
	TTree *t = (TTree*)f->Get("jets_ct10_10gev_0130");
	
	TCanvas *c1 = new TCanvas("c1","R32vsQ");
	
	int N = t->GetEntries();
	
	cout << "\nNumber of entries: " << N << endl;
	cout << endl;
	
	// variables to store branch content
	double jet2pt;
	double jet3pt;
	double Q;
	
	t->SetBranchAddress("jet2pt",&jet2pt);
	t->SetBranchAddress("jet3pt",&jet3pt);
	t->SetBranchAddress("Q",&Q);
	
	double Qmin = 20.0; // Q with 8 bins wide 2.5 GeV in [20.0,40.0] GeV
	double Qmax = 40.0; // Q with 6 bins wide 10 GeV in (40.0,100.0] GeV
	
	cout << "\nQmin is: " << Qmin << endl;
	cout << "Qmax is: " << Qmax << endl;
	
	int nbins = 8;
	
	double deltabin = (Qmax-Qmin)/nbins;
	
	cout << "\nBin width is: " << deltabin << endl;
	cout << endl;
	
	// arrays of counters: I want a counter for each of the 8+6=14 bins
	double count2[14] = {0.0}; // ATTENTION! MANUALLY INSERT nbins!
	double count3[14] = {0.0};
	double Qcentr[14] = {0.0};
	
	for (int i=0; i<N; i++) {
		
		t->GetEntry(i);
		
		// I want to consider only events with at least two jets
		if (jet2pt == 0.0) continue;
		
		for (int k=0; k<nbins ; k++) { // k runs from 0 to nbins, ok
			
			if (Q>=(Qmin + k*deltabin) && Q<(Qmin + (k+1)*deltabin)) {
				
				if (jet2pt>20) {
					
					count2[k]++;
					
				}
				
				if (jet3pt>20) {
					
					count3[k]++;
					
				}
				
			}
			
			Qcentr[k] = Qmin + (k+0.5)*deltabin; // I compute the centroid of each bin
			
		}
		
		if (Q>40 && Q<=50) {
			if (jet2pt>20) {
				count2[8]++;
			}
			if (jet3pt>20) {	
				count3[8]++;
			}
		}
		
		if (Q>50 && Q<=60) {
			if (jet2pt>20) {
				count2[9]++;
			}
			if (jet3pt>20) {	
				count3[9]++;
			}
		}
		
		if (Q>60 && Q<=70) {
			if (jet2pt>20) {
				count2[10]++;
			}
			if (jet3pt>20) {	
				count3[10]++;
			}
		}
		
		if (Q>70 && Q<=80) {
			if (jet2pt>20) {
				count2[11]++;
			}
			if (jet3pt>20) {	
				count3[11]++;
			}
		}
		
		if (Q>80 && Q<=90) {
			if (jet2pt>20) {
				count2[12]++;
			}
			if (jet3pt>20) {	
				count3[12]++;
			}
		}
		
		if (Q>90 && Q<=100) {
			if (jet2pt>20) {
				count2[13]++;
			}
			if (jet3pt>20) {	
				count3[13]++;
			}
		}
		
	}
	
	Qcentr[8] = (40 + 50)/2;
	Qcentr[9] = (50 + 60)/2;
	Qcentr[10] = (60 + 70)/2;
	Qcentr[11] = (70 + 80)/2;
	Qcentr[12] = (80 + 90)/2;
	Qcentr[13] = (90 + 100)/2;
	
	for (int n2=0; n2<nbins+6; n2++) {
		
		cout << "Number of events with 2 jets in Q bin number " << n2+1 << ": " << count2[n2] << endl;
		
	}
	
	cout << endl;
	
	for (int n3=0; n3<nbins+6; n3++) {
		
		cout << "Number of events with 3 jets in Q bin number " << n3+1 << ": " << count3[n3] << endl;
		
	}
	
	cout << endl;
	
	int sum2 = count2[0] + count2[1] + count2[2] + count2[3] + count2[4] + count2[5] + count2[6] + count2[7] + count2[8] + count2[9] + count2[10] + count2[11] + count2[12] + count2[13];
	
	cout << "The sum of events with 2 jets is: " << sum2 << endl;
	
	int sum3 = count3[0] + count3[1] + count3[2] + count3[3] + count3[4] + count3[5] + count3[6] + count3[7] + count3[8] + count3[9] + count3[10] + count3[11] + count3[12] + count3[13];
	
	cout << "The sum of events with 3 jets is: " << sum3 << endl;
	cout << endl;
	
	for (int c=0; c<nbins+6; c++) {
		
		cout << "Centroid of bin number " << c+1 << " is:\t" << Qcentr[c] << endl;
		
	}
	
	cout << endl;
	
	// array for R_32 in each bin
	double R_32[14] = {0.0}; // ATTENTION! MANUALLY INSERT nbins!
	double R_32err[14] = {0.0};
	
	for (int r=0; r<nbins+6; r++) {
		
		if (count2[r]==0) {
			R_32[r] = 0;
		}
		if (count2[r]!=0) {
			R_32[r] = count3[r]/count2[r];
			R_32err[r] = (1/pow(count2[r],2)) * sqrt(count2[r]*pow(count3[r],2) + count3[r]*pow(count2[r],2));
		}
		
		cout << "R_32 in bin number " << r+1 << " is:\t" << R_32[r] << endl; 
		
	}
	
	ofstream outfile ("R32vsQ_10gev_0130.txt");
	
	if (outfile.is_open()) {
		
		for (int f=0; f<nbins+6; f++) {
			
			outfile << Qcentr[f] << "\t" << R_32[f] << "\t" << R_32err[f] << endl;
			
		}
		
		outfile.close();
		
	}
	
	else cout << "Unable to open the output file";
	
	TGraphErrors *R32vsQ = new TGraphErrors("R32vsQ_10gev_0130.txt","%lg %lg %lg");
	
	c1->cd();
	
	R32vsQ->SetMarkerColor(2);
	R32vsQ->SetMarkerStyle(20);
	R32vsQ->SetMarkerSize(1);
	
	R32vsQ->SetTitle("R_{32} vs Q - #alpha_{s} = 0.110, p_{T, min}(j) = 10 GeV, p_{T}(jet) = 20 GeV;Q [GeV];R_{32}"); // ("title;xlabel;ylabel")
	R32vsQ->GetXaxis()->SetTitleSize(0.05);
	R32vsQ->GetYaxis()->SetTitleSize(0.05);
	R32vsQ->GetYaxis()->SetTickLength(0.02);
	
//	R32vsQ->GetXaxis()->SetRangeUser(5,47); // xrange
//	R32vsQ->GetYaxis()->SetRangeUser(0,1.2); // yrange
	
	R32vsQ->Draw("AP");
	
	c1->SetGridx();
	c1->SetGridy();
	
	
	
	f->Close();
	
}

void R32vsm12() {
	
	TFile *f = new TFile("jets_ct10_10gev_0130.root","READ");
	
	TTree *t = (TTree*)f->Get("jets_ct10_10gev_0130");
	
	TCanvas *c2 = new TCanvas("c2","R32vsm12");
	
	int N = t->GetEntries();
	
	cout << "\nNumber of entries: " << N << endl;
	
	// variables to store branches contents
	double jet1mass;
	double jet2mass;
	double jet1px;
	double jet1py;
	double jet1pz;
	double jet2px;
	double jet2py;
	double jet2pz;
	double jet2pt;
	double jet3pt;
		
	t->SetBranchAddress("jet1mass",&jet1mass);
	t->SetBranchAddress("jet2mass",&jet2mass);
	t->SetBranchAddress("jet1px",&jet1px);
	t->SetBranchAddress("jet1py",&jet1py);
	t->SetBranchAddress("jet1pz",&jet1pz);
	t->SetBranchAddress("jet2px",&jet2px);
	t->SetBranchAddress("jet2py",&jet2py);
	t->SetBranchAddress("jet2pz",&jet2pz);
	t->SetBranchAddress("jet2pt",&jet2pt);
	t->SetBranchAddress("jet3pt",&jet3pt);
	
	// array to store mass12 values
	double mass12[622422] = {0.0}; // 0.110: 619972, 0.118:621619, 0.130: 622422
	
	for (int i=0; i<N; i++) {
		
		t->GetEntry(i);
		mass12[i] = pow(((jet1mass*jet1mass) + (jet2mass*jet2mass) + (2*pow((((jet1mass*jet1mass) + (jet1px*jet1px) + (jet1py*jet1py) + (jet1pz*jet1pz))*((jet2mass*jet2mass) + (jet2px*jet2px) + (jet2py*jet2py) + (jet2pz*jet2pz))),0.5)) - (2*((jet1px*jet2px) + (jet1py*jet2py) + (jet1pz*jet2pz)))),0.5);
	}
	
	double mass12min = 0.0;   // m12 with 12 bins wide 10 GeV in [0.0,120.0] GeV
	double mass12max = 120.0; // m12 with 4 bins wide 20 GeV in (120.0,200.0] GeV
	
	cout << "\nmass12_min is: " << mass12min << endl;
	cout << "mass12_max is: " << mass12max << endl;
	
	int nbins = 12;
	
	double deltabin = (mass12max-mass12min)/nbins;
	
	cout << "\nBin width is: " << deltabin << endl;
	cout << endl;
	
	// arrays of counters: I want a counter for each of the 12+4=16 bins
	double count2[16] = {0.0}; // ATTENTION! MANUALLY INSERT nbins!
	double count3[16] = {0.0};
	double mass12centr[16] = {0.0};
	
	for (int i=0; i<N; i++) {
		
		t->GetEntry(i);
		
		// I want to consider only events with at least two jets
		if (jet2pt == 0.0) continue;
		
		for (int k=0; k<nbins ; k++) { // k runs from 0 to nbins, ok
			
			if (mass12[i]>=(mass12min + k*deltabin) && mass12[i]<(mass12min + (k+1)*deltabin)) {
				
				if (jet2pt>20) {
					
					count2[k]++;
					
				}
				
				if (jet3pt>20) {
					
					count3[k]++;
					
				}
				
			}
			
			mass12centr[k] = mass12min + (k+0.5)*deltabin; // here I compute the centroid of each bin
			
		}
		
		if (mass12[i]>120 && mass12[i]<=140) {
			if (jet2pt>20) {
				count2[12]++;
			}
			if (jet3pt>20) {	
				count3[12]++;
			}
		}
		
		if (mass12[i]>140 && mass12[i]<=160) {
			if (jet2pt>20) {
				count2[13]++;
			}
			if (jet3pt>20) {	
				count3[13]++;
			}
		}
		
		if (mass12[i]>160 && mass12[i]<=180) {
			if (jet2pt>20) {
				count2[14]++;
			}
			if (jet3pt>20) {	
				count3[14]++;
			}
		}
		
		if (mass12[i]>180 && mass12[i]<=200) {
			if (jet2pt>20) {
				count2[15]++;
			}
			if (jet3pt>20) {	
				count3[15]++;
			}
		}
		
	}
	
	mass12centr[12] = (120 + 140)/2;
	mass12centr[13] = (140 + 160)/2;
	mass12centr[14] = (160 + 180)/2;
	mass12centr[15] = (180 + 200)/2;
	
	for (int n2=0; n2<nbins+4; n2++) {
		
		cout << "Number of events with 2 jets in bin number " << n2+1 << ": " << count2[n2] << endl;
		
	}
	
	cout << endl;
	
	for (int n3=0; n3<nbins+4; n3++) {
		
		cout << "Number of events with 3 jets in bin number " << n3+1 << ": " << count3[n3] << endl;
		
	}
	
	cout << endl;
	
	int sum2 = count2[0] + count2[1] + count2[2] + count2[3] + count2[4] + count2[5] + count2[6] + count2[7] + count2[8] + count2[9] + count2[10] + count2[11] + count2[12] + count2[13] + count2[14] + count2[15];
	
	cout << "The sum of events with 2 jets is: " << sum2 << endl;
	
	int sum3 = count3[0] + count3[1] + count3[2] + count3[3] + count3[4] + count3[5] + count3[6] + count3[7] + count3[8] + count3[9] + count3[10] + count3[11] + count3[12] + count3[13] + count3[14] + count3[15];
	
	cout << "The sum of events with 3 jets is: " << sum3 << endl;
	cout << endl;
	
	for (int c=0; c<nbins+4; c++) {
		
		cout << "Centroid of bin number " << c+1 << " is:\t" << mass12centr[c] << endl;
		
	}
	
	cout << endl;
	
	// array for R_32 in each bin
	double R_32[16] = {0.0}; // ATTENTION! MANUALLY INSERT nbins!
	double R_32err[16] = {0.0};
	
	for (int r=0; r<nbins+4; r++) {
		
		if (count2[r]==0) {
			R_32[r] = 0;
		}
		if (count2[r]!=0) {
			R_32[r] = count3[r]/count2[r];
			R_32err[r] = (1/pow(count2[r],2)) * sqrt(count2[r]*pow(count3[r],2) + count3[r]*pow(count2[r],2));
		}
		
		cout << "R_32 in bin number " << r+1 << " is:\t" << R_32[r] << endl; 
		
	}
	
	ofstream outfile ("R32vsm12_10gev_0130.txt");
	
	if (outfile.is_open()) {
		
		for (int f=0; f<nbins+4; f++) {
			
			outfile << mass12centr[f] << "\t" << R_32[f] << "\t" << R_32err[f] << endl;
			
		}
		
		outfile.close();
		
	}
	
	else cout << "Unable to open the output file";
	
	TGraphErrors *R32vsmass12 = new TGraphErrors("R32vsm12_10gev_0130.txt","%lg %lg %lg");
	
	c2->cd();
	
	R32vsmass12->SetMarkerColor(2);
	R32vsmass12->SetMarkerStyle(20);
	R32vsmass12->SetMarkerSize(1);
	
	R32vsmass12->SetTitle("R_{32} vs m_{12} - #alpha_{s} = 0.110, p_{T, min}(j) = 10 GeV, p_{T}(jet) = 20 GeV;m_{12} [GeV];R_{32}"); // ("title;xlabel;ylabel")
	R32vsmass12->GetXaxis()->SetTitleSize(0.05);
	R32vsmass12->GetYaxis()->SetTitleSize(0.05);
	R32vsmass12->GetYaxis()->SetTickLength(0.02);
	
//	R32vsmass12->GetXaxis()->SetRangeUser(5,47); // xrange
//	R32vsmass12->GetYaxis()->SetRangeUser(0,1.2); // yrange
	
	R32vsmass12->Draw("AP");
	
	c2->SetGridx();
	c2->SetGridy();
	
	
	
	f->Close();
	
}
