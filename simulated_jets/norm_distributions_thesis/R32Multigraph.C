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

void multigraph() {
	
	TGraphErrors *g110 = new TGraphErrors("R32vsQ_10gev_0110.txt","%lg %lg %lg");
	TGraphErrors *g118 = new TGraphErrors("R32vsQ_10gev_0118.txt","%lg %lg %lg");
	TGraphErrors *g130 = new TGraphErrors("R32vsQ_10gev_0130.txt","%lg %lg %lg");
	
	TCanvas *c = new TCanvas("c","R32vsQ");
	
	c->cd();
	
	g110->SetName("g110");
	g118->SetName("g118");
	g130->SetName("g130");
	
	g110->SetLineColor(kBlue+2);
	g110->SetMarkerColor(kBlue+2);
	g110->SetMarkerStyle(20);
	g110->SetMarkerSize(0.8);
	g110->SetFillColor(kBlue+2);
	g110->SetFillStyle(3003);
	
	g118->SetLineColor(kGreen+1);
	g118->SetMarkerColor(kGreen+2);
	g118->SetMarkerStyle(20);
	g118->SetMarkerSize(0.8);
	g118->SetFillColor(kGreen+2);
	g118->SetFillStyle(3003);
	
	g130->SetLineColor(kRed);
	g130->SetMarkerColor(kRed);
	g130->SetMarkerStyle(20);
	g130->SetMarkerSize(0.8);
	g130->SetFillColor(kRed);
	g130->SetFillStyle(3003);
	
	TMultiGraph *mg = new TMultiGraph();
	
	mg->SetTitle("R_{32} vs Q - p_{T, min}(j) = 10 GeV, p_{T}(jet) > 20 GeV;Q [GeV];R_{32}");
	mg->Add(g110,"4P");
	mg->Add(g130,"4P");
	mg->Add(g118,"4P");
	
	mg->Draw("A");
	
	TLegend leg(.1,.7,.5,.9); // (x1,y1,x2,y2)
	
	leg.AddEntry("g110","CT10nnlo #alpha_{s}=0.1100","FP");
	leg.AddEntry("g118","CT10nnlo #alpha_{s}=0.1180","FP");
	leg.AddEntry("g130","CT10nnlo #alpha_{s}=0.1300","FP");
	leg.DrawClone("SAME");
	
}
