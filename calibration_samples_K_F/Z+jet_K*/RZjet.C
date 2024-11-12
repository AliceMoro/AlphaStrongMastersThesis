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

void RZjet() { // R = pt(jet)/pt(Z)
	
	double jetpt_min = 2e4; // MeV
	double jetpt_max = 10e4; // MeV
	
	TFile *fdata = new TFile("Zmumu_data.root");
	TFile *fmc = new TFile("Zmumu_mc.root");
	
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
	
	TH1D *hRdata = new TH1D("hRdata","Calibration sample Z+jet - R_{Z} = p_{T}(jet) / p_{T}(Z^{0}) - 20 GeV < p_{T}(jet) < 100 GeV;R_{Z};A.U.",100,0,10);
	TH1D *hRmc = new TH1D("hRmc","Calibration sample Z+jet - R_{Z} = p_{T}(jet) / p_{T}(Z^{0}) - 20 GeV < p_{T}(jet) < 100 GeV;R_{Z};A.U.",100,0,10);
	
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
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
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
			
			Rmc[m] = jetpt_mc / Zpt_mc;
			
			hRmc->Fill(Rmc[m]);
			
		}
		
	}
	
	double Nsel_data = hRdata->Integral();
	double Nsel_mc = hRmc->Integral();
		
	cout << "\nNsel data = " << Nsel_data << endl;
	cout << "\nNsel MC = " << Nsel_mc << "\n" << endl;
	
	hRdata->Scale(1/(hRdata->Integral()));
	hRmc->Scale(1/(hRmc->Integral()));
	
	TCanvas *c1 = new TCanvas("c1","RZjet");
	
	hRdata->GetXaxis()->SetRangeUser(0,5);
	hRdata->SetLineColor(kRed);
	hRdata->SetMarkerStyle(20);
	hRdata->SetMarkerSize(1.2);
	hRdata->SetMarkerColor(kRed);
	hRmc->GetXaxis()->SetRangeUser(0,5);
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
	
}








void phiZjet() {
	
	TFile *fdata = new TFile("Zmumu_data.root");
	TFile *fmc = new TFile("Zmumu_mc.root");
	
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
	
	double phidata[190849] = {0.0};
	double phimc[374981] = {0.0};
	
	TH1D *hphidata = new TH1D("hphidata","Calibration sample Z^{0}+jet - #Delta#phi_{Z^{0}+jet} - data;#Delta#phi_{Z^{0}+jet};Events",35,0,3.5);
	TH1D *hphimc = new TH1D("hphimc","Calibration sample Z^{0}+jet - #Delta#phi_{Z^{0}+jet} - MC;#Delta#phi_{Z^{0}+jet};Events",35,0,3.5);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double mupeta_data = -log(tan(acos(muppz_data/sqrt(muppx_data*muppx_data+muppy_data*muppy_data+muppz_data*muppz_data))/2));
		
		double mumeta_data = -log(tan(acos(mumpz_data/sqrt(mumpx_data*mumpx_data+mumpy_data*mumpy_data+mumpz_data*mumpz_data))/2));
		
		if (totCand_data != 1) continue; // select events with only one jet
		
		if (jetpt_data < 2e4 || jetpt_data > 1e5) continue; // cut in pT of jet
		
		if (jeteta_data < 2.2 || jeteta_data > 4.2) continue; // cut in eta of jet
		
		if (muppt_data < 2e4 || mumpt_data < 2e4) continue; // cut in pT of muons
		
		if (mupeta_data < 2 || mupeta_data > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_data < 2 || mumeta_data > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_data == 1 && mupHLT1MuonHighPt_data == 1 && mupHLT2_data == 1) || (mumL0MuonEW_data == 1 && mumHLT1MuonHighPt_data == 1 && mumHLT2_data == 1)) { // cuts in trigger lines on muons
		
			phidata[d] = acos(((Zpx_data*jetpx_data) + (Zpy_data*jetpy_data))/(pow((((Zpx_data*Zpx_data) + (Zpy_data*Zpy_data))*((jetpx_data*jetpx_data) + (jetpy_data*jetpy_data))),0.5)));
			
			hphidata->Fill(phidata[d]);
		
		}
		
	}
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
		
		double mupeta_mc = -log(tan(acos(muppz_mc/sqrt(muppx_mc*muppx_mc+muppy_mc*muppy_mc+muppz_mc*muppz_mc))/2));
			
		double mumeta_mc = -log(tan(acos(mumpz_mc/sqrt(mumpx_mc*mumpx_mc+mumpy_mc*mumpy_mc+mumpz_mc*mumpz_mc))/2));
		
		if (totCand_mc != 1) continue; // select events with only one jet
		
		if (jetpt_mc < 2e4 || jetpt_mc > 1e5) continue; // cut in pT of jet
		
		if (jeteta_mc < 2.2 || jeteta_mc > 4.2) continue; // cut in eta of jet
		
		if (muppt_mc < 2e4 || mumpt_mc < 2e4) continue; // cut in pT of muons
		
		if (mupeta_mc < 2 || mupeta_mc > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_mc < 2 || mumeta_mc > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_mc == 1 && mupHLT1MuonHighPt_mc == 1 && mupHLT2_mc == 1) || (mumL0MuonEW_mc == 1 && mumHLT1MuonHighPt_mc == 1 && mumHLT2_mc == 1)) { // cuts in trigger lines on muons
			
			phimc[m] = acos(((Zpx_mc*jetpx_mc) + (Zpy_mc*jetpy_mc))/(pow((((Zpx_mc*Zpx_mc) + (Zpy_mc*Zpy_mc))*((jetpx_mc*jetpx_mc) + (jetpy_mc*jetpy_mc))),0.5)));
			
			hphimc->Fill(phimc[m]);
		
		}
		
	}
	
	double Nsel_data = hphidata->Integral();
	double Nsel_mc = hphimc->Integral();
	
	cout << "\nNsel data = " << Nsel_data << endl;
	cout << "\nNsel MC = " << Nsel_mc << "\n" << endl;
	
	hphidata->Scale(1/(hphidata->Integral()));
	hphimc->Scale(1/(hphimc->Integral()));
	
	TCanvas *c2 = new TCanvas("c2","phiZjet");
	
	hphidata->GetXaxis()->SetRangeUser(0,3.4);
	hphimc->GetXaxis()->SetRangeUser(0,3.4);
	hphidata->GetYaxis()->SetRangeUser(0,0.23);
	hphimc->GetYaxis()->SetRangeUser(0,0.23);
	hphidata->SetTitle("Calibration sample Z+jet - #Delta#phi_{Z+jet};#Delta#phi_{Z+jet};A.U.");
	hphimc->SetTitle("Calibration sample Z+jet - #Delta#phi_{Z+jet};#Delta#phi_{Z+jet};A.U.");
	hphidata->SetLineColor(kRed);
	hphidata->SetMarkerStyle(20);
	hphidata->SetMarkerSize(1);
	hphidata->SetMarkerColor(kRed);
	hphimc->SetLineColor(kBlue+2);
	hphimc->SetMarkerStyle(20);
	hphimc->SetMarkerSize(1);
	hphimc->SetMarkerColor(kBlue+2);
	
	TLine *xline = new TLine(2.8,0,2.8,0.23);
	
	xline->SetLineColor(kBlack);
	xline->SetLineWidth(3);
	
	gStyle->SetOptStat("");
	
	c2->cd();
	
	hphidata->Draw();
	hphimc->Draw("SAME");
	xline->Draw("SAME");
	
	TLegend leg2(.1,.65,.25,.9); // (x1,y1,x2,y2)
	
	leg2.AddEntry(hphidata,"Data","PL");
	leg2.AddEntry(hphimc,"MC","PL");
	leg2.AddEntry(xline,"Cut","L");
	leg2.DrawClone("SAME");
	
}








void GetJetpT() {
	
	TFile *fdata = new TFile("Zmumu_data.root");
	TFile *fmc = new TFile("Zmumu_mc.root");
	
	fdata->cd("DecayTreeTuple");
	TTree *tdata = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	fmc->cd("DecayTreeTuple");
	TTree *tmc = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	int Ndata = tdata->GetEntries();
	int Nmc = tmc->GetEntries();
	
	cout << "\nNumber of entries of data tuple = " << Ndata << "\nNumber of entries of MC tuple   = " << Nmc << "\n" << endl;
	
	double jetpt_data = 0.0;
	double jeteta_data = 0.0;
	
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
	
	TH1D *hpTdata = new TH1D("hpTdata","Calibration sample Z^{0}+jet - p_{T}(jet) - data;p_{T}(jet) [GeV];A.U.",60,0,300);
	TH1D *hpTmc = new TH1D("hpTmc","Calibration sample Z^{0}+jet - p_{T}(jet) - MC;p_{T}(jet) [GeV];A.U.",60,0,300);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double phi_data = acos(((Zpx_data*jetpx_data) + (Zpy_data*jetpy_data))/(pow((((Zpx_data*Zpx_data) + (Zpy_data*Zpy_data))*((jetpx_data*jetpx_data) + (jetpy_data*jetpy_data))),0.5)));
		
		double mupeta_data = -log(tan(acos(muppz_data/sqrt(muppx_data*muppx_data+muppy_data*muppy_data+muppz_data*muppz_data))/2));
		
		double mumeta_data = -log(tan(acos(mumpz_data/sqrt(mumpx_data*mumpx_data+mumpy_data*mumpy_data+mumpz_data*mumpz_data))/2));
		
//		if (jetpt_data < 1e4 || jetpt_data > 1e5) continue; // cut in pT of jet
		
		if (totCand_data != 1) continue; // select events with only one jet
		
		if (jeteta_data < 2.2 || jeteta_data > 4.2) continue; // cut in eta of jet
		
		if (phi_data < 2.8) continue; // cut in DeltaPhi
		
		if (muppt_data < 2e4 || mumpt_data < 2e4) continue; // cut in pT of muons
		
		if (mupeta_data < 2 || mupeta_data > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_data < 2 || mumeta_data > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_data == 1 && mupHLT1MuonHighPt_data == 1 && mupHLT2_data == 1) || (mumL0MuonEW_data == 1 && mumHLT1MuonHighPt_data == 1 && mumHLT2_data == 1)) { // cuts in trigger lines on muons
		
			hpTdata->Fill(jetpt_data/1000);
			
		}
		
	}
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
			
		double phi_mc = acos(((Zpx_mc*jetpx_mc) + (Zpy_mc*jetpy_mc))/(pow((((Zpx_mc*Zpx_mc) + (Zpy_mc*Zpy_mc))*((jetpx_mc*jetpx_mc) + (jetpy_mc*jetpy_mc))),0.5)));
		
		double mupeta_mc = -log(tan(acos(muppz_mc/sqrt(muppx_mc*muppx_mc+muppy_mc*muppy_mc+muppz_mc*muppz_mc))/2));
		
		double mumeta_mc = -log(tan(acos(mumpz_mc/sqrt(mumpx_mc*mumpx_mc+mumpy_mc*mumpy_mc+mumpz_mc*mumpz_mc))/2));
		
//		if (jetpt_mc < 1e4 || jetpt_mc > 1e5) continue; // cut in pT of jet

		if (totCand_mc != 1) continue; // select events with only one jet
		
		if (jeteta_mc < 2.2 || jeteta_mc > 4.2) continue; // cut in eta of jet
		
		if (phi_mc < 2.8) continue; // cut in DeltaPhi
			
		if (muppt_mc < 2e4 || mumpt_mc < 2e4) continue; // cut in pT of muons
		
		if (mupeta_mc < 2 || mupeta_mc > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_mc < 2 || mumeta_mc > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_mc == 1 && mupHLT1MuonHighPt_mc == 1 && mupHLT2_mc == 1) || (mumL0MuonEW_mc == 1 && mumHLT1MuonHighPt_mc == 1 && mumHLT2_mc == 1)) { // cuts in trigger lines on muons
			
			hpTmc->Fill(jetpt_mc/1000);
			
		}
		
	}
	
	double Nsel_data = hpTdata->Integral();
	double Nsel_mc = hpTmc->Integral();
		
	cout << "\nNsel data = " << Nsel_data << endl;
	cout << "\nNsel MC = " << Nsel_mc << "\n" << endl;
	
	hpTdata->Scale(1/(hpTdata->Integral()));
	hpTmc->Scale(1/(hpTmc->Integral()));
	
	TCanvas *c3 = new TCanvas("c3","pTjet");
	
	hpTdata->SetLineColor(kRed);
	hpTmc->SetTitle("Calibration sample Z^{0}+jet - p_{T}(jet);p_{T}(jet) [GeV];A.U.");
	hpTdata->SetTitle("Calibration sample Z^{0}+jet - p_{T}(jet);p_{T}(jet) [GeV];A.U.");
	gStyle->SetOptStat("");
	
	c3->cd();
	
	hpTmc->Draw();
	hpTdata->Draw("SAME");
	
	TLegend leg3(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg3.AddEntry("hpTdata","Data");
	leg3.AddEntry("hpTmc","MC");
	leg3.DrawClone("SAME");
	
}








void GetJeteta() {
	
	TFile *fdata = new TFile("Zmumu_data.root");
	TFile *fmc = new TFile("Zmumu_mc.root");
	
	fdata->cd("DecayTreeTuple");
	TTree *tdata = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	fmc->cd("DecayTreeTuple");
	TTree *tmc = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	int Ndata = tdata->GetEntries();
	int Nmc = tmc->GetEntries();
	
	cout << "\nNumber of entries of data tuple = " << Ndata << "\nNumber of entries of MC tuple   = " << Nmc << "\n" << endl;
	
	double jetpt_data = 0.0;
	double jeteta_data = 0.0;
	
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
	
	TH1D *hetadata = new TH1D("hetadata","Calibration sample Z^{0}+jet - #eta(jet) - data;#eta(jet);A.U.",35,1.5,5);
	TH1D *hetamc = new TH1D("hetamc","Calibration sample Z^{0}+jet - #eta(jet) - MC;#eta(jet);A.U.",35,1.5,5);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double phi_data = acos(((Zpx_data*jetpx_data) + (Zpy_data*jetpy_data))/(pow((((Zpx_data*Zpx_data) + (Zpy_data*Zpy_data))*((jetpx_data*jetpx_data) + (jetpy_data*jetpy_data))),0.5)));
		
		double mupeta_data = -log(tan(acos(muppz_data/sqrt(muppx_data*muppx_data+muppy_data*muppy_data+muppz_data*muppz_data))/2));
		
		double mumeta_data = -log(tan(acos(mumpz_data/sqrt(mumpx_data*mumpx_data+mumpy_data*mumpy_data+mumpz_data*mumpz_data))/2));
		
		if (totCand_data != 1) continue; // select events with only one jet
		
		if (jetpt_data < 1e4 || jetpt_data > 1e5) continue; // cut in pT of jet
		
//		if (jeteta_data < 2.2 || jeteta_data > 4.2) continue; // cut in eta of jet
		
		if (phi_data < 2.8) continue; // cut in DeltaPhi
		
		if (muppt_data < 2e4 || mumpt_data < 2e4) continue; // cut in pT of muons
		
		if (mupeta_data < 2 || mupeta_data > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_data < 2 || mumeta_data > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_data == 1 && mupHLT1MuonHighPt_data == 1 && mupHLT2_data == 1) || (mumL0MuonEW_data == 1 && mumHLT1MuonHighPt_data == 1 && mumHLT2_data == 1)) { // cuts in trigger lines on muons
		
			hetadata->Fill(jeteta_data);
			
		}
		
	}
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
			
		double phi_mc = acos(((Zpx_mc*jetpx_mc) + (Zpy_mc*jetpy_mc))/(pow((((Zpx_mc*Zpx_mc) + (Zpy_mc*Zpy_mc))*((jetpx_mc*jetpx_mc) + (jetpy_mc*jetpy_mc))),0.5)));
		
		double mupeta_mc = -log(tan(acos(muppz_mc/sqrt(muppx_mc*muppx_mc+muppy_mc*muppy_mc+muppz_mc*muppz_mc))/2));
		
		double mumeta_mc = -log(tan(acos(mumpz_mc/sqrt(mumpx_mc*mumpx_mc+mumpy_mc*mumpy_mc+mumpz_mc*mumpz_mc))/2));
		
		if (totCand_mc != 1) continue; // select events with only one jet
		
		if (jetpt_mc < 1e4 || jetpt_mc > 1e5) continue; // cut in pT of jet

//		if (jeteta_mc < 2.2 || jeteta_mc > 4.2) continue; // cut in eta of jet
		
		if (phi_mc < 2.8) continue; // cut in DeltaPhi
			
		if (muppt_mc < 2e4 || mumpt_mc < 2e4) continue; // cut in pT of muons
		
		if (mupeta_mc < 2 || mupeta_mc > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_mc < 2 || mumeta_mc > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_mc == 1 && mupHLT1MuonHighPt_mc == 1 && mupHLT2_mc == 1) || (mumL0MuonEW_mc == 1 && mumHLT1MuonHighPt_mc == 1 && mumHLT2_mc == 1)) { // cuts in trigger lines on muons
			
			hetamc->Fill(jeteta_mc);
			
		}
		
	}
	
	double Nsel_data = hetadata->Integral();
	double Nsel_mc = hetamc->Integral();
		
	cout << "\nNsel data = " << Nsel_data << endl;
	cout << "\nNsel MC = " << Nsel_mc << "\n" << endl;
	
	hetadata->Scale(1/(hetadata->Integral()));
	hetamc->Scale(1/(hetamc->Integral()));
	
	TCanvas *c4 = new TCanvas("c4","etajet");
	
	hetadata->SetLineColor(kRed);
	hetamc->SetTitle("Calibration sample Z^{0}+jet - #eta(jet);#eta(jet);A.U.");
	hetadata->SetTitle("Calibration sample Z^{0}+jet - #eta(jet);#eta(jet);A.U.");
	gStyle->SetOptStat("");
	
	c4->cd();
	
	hetamc->Draw();
	hetadata->Draw("SAME");
	
	TLegend leg4(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg4.AddEntry("hetadata","Data");
	leg4.AddEntry("hetamc","MC");
	leg4.DrawClone("SAME");
	
}








void GetJetphi() {
	
	TFile *fdata = new TFile("Zmumu_data.root");
	TFile *fmc = new TFile("Zmumu_mc.root");
	
	fdata->cd("DecayTreeTuple");
	TTree *tdata = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	fmc->cd("DecayTreeTuple");
	TTree *tmc = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	int Ndata = tdata->GetEntries();
	int Nmc = tmc->GetEntries();
	
	cout << "\nNumber of entries of data tuple = " << Ndata << "\nNumber of entries of MC tuple   = " << Nmc << "\n" << endl;
	
	double jetpt_data = 0.0;
	double jeteta_data = 0.0;
	
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
	
	TH1D *hphidata = new TH1D("hphidata","Calibration sample Z^{0}+jet - #phi(jet) - data;#phi(jet);A.U.",32,0,3.2);
	TH1D *hphimc = new TH1D("hphimc","Calibration sample Z^{0}+jet - #phi(jet) - MC;#phi(jet);A.U.",32,0,3.2);
	
	for (int d=0; d<Ndata+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double phi_data = acos(((Zpx_data*jetpx_data) + (Zpy_data*jetpy_data))/(pow((((Zpx_data*Zpx_data) + (Zpy_data*Zpy_data))*((jetpx_data*jetpx_data) + (jetpy_data*jetpy_data))),0.5)));
		
		double mupeta_data = -log(tan(acos(muppz_data/sqrt(muppx_data*muppx_data+muppy_data*muppy_data+muppz_data*muppz_data))/2));
		
		double mumeta_data = -log(tan(acos(mumpz_data/sqrt(mumpx_data*mumpx_data+mumpy_data*mumpy_data+mumpz_data*mumpz_data))/2));
		
		if (totCand_data != 1) continue; // select events with only one jet
		
		if (jetpt_data < 1e4 || jetpt_data > 1e5) continue; // cut in pT of jet
		
		if (jeteta_data < 2.2 || jeteta_data > 4.2) continue; // cut in eta of jet
		
		if (phi_data < 2.8) continue; // cut in DeltaPhi
		
		if (muppt_data < 2e4 || mumpt_data < 2e4) continue; // cut in pT of muons
		
		if (mupeta_data < 2 || mupeta_data > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_data < 2 || mumeta_data > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_data == 1 && mupHLT1MuonHighPt_data == 1 && mupHLT2_data == 1) || (mumL0MuonEW_data == 1 && mumHLT1MuonHighPt_data == 1 && mumHLT2_data == 1)) { // cuts in trigger lines on muons
		
			hphidata->Fill(acos(jetpx_data/(pow((jetpx_data*jetpx_data+jetpy_data*jetpy_data),0.5))));
			
		}
		
	}
	
	for (int m=0; m<Nmc+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
			
		double phi_mc = acos(((Zpx_mc*jetpx_mc) + (Zpy_mc*jetpy_mc))/(pow((((Zpx_mc*Zpx_mc) + (Zpy_mc*Zpy_mc))*((jetpx_mc*jetpx_mc) + (jetpy_mc*jetpy_mc))),0.5)));
		
		double mupeta_mc = -log(tan(acos(muppz_mc/sqrt(muppx_mc*muppx_mc+muppy_mc*muppy_mc+muppz_mc*muppz_mc))/2));
		
		double mumeta_mc = -log(tan(acos(mumpz_mc/sqrt(mumpx_mc*mumpx_mc+mumpy_mc*mumpy_mc+mumpz_mc*mumpz_mc))/2));
		
		if (totCand_mc != 1) continue; // select events with only one jet
		
		if (jetpt_mc < 1e4 || jetpt_mc > 1e5) continue; // cut in pT of jet

		if (jeteta_mc < 2.2 || jeteta_mc > 4.2) continue; // cut in eta of jet
		
		if (phi_mc < 2.8) continue; // cut in DeltaPhi
			
		if (muppt_mc < 2e4 || mumpt_mc < 2e4) continue; // cut in pT of muons
		
		if (mupeta_mc < 2 || mupeta_mc > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_mc < 2 || mumeta_mc > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_mc == 1 && mupHLT1MuonHighPt_mc == 1 && mupHLT2_mc == 1) || (mumL0MuonEW_mc == 1 && mumHLT1MuonHighPt_mc == 1 && mumHLT2_mc == 1)) { // cuts in trigger lines on muons
			
			hphimc->Fill(acos(jetpx_mc/(pow((jetpx_mc*jetpx_mc+jetpy_mc*jetpy_mc),0.5))));
			
		}
		
	}
	
	double Nsel_data = hphidata->Integral();
	double Nsel_mc = hphimc->Integral();
		
	cout << "\nNsel data = " << Nsel_data << endl;
	cout << "\nNsel MC = " << Nsel_mc << "\n" << endl;
	
	hphidata->Scale(1/(hphidata->Integral()));
	hphimc->Scale(1/(hphimc->Integral()));
	
	TCanvas *c4 = new TCanvas("c4","etajet");
	
	hphidata->SetLineColor(kRed);
	hphimc->SetTitle("Calibration sample Z^{0}+jet - #phi(jet);#phi(jet);A.U.");
	hphidata->SetTitle("Calibration sample Z^{0}+jet - #phi(jet);#phi(jet);A.U.");
	gStyle->SetOptStat("");
	
	c4->cd();
	
	hphidata->Draw();
	hphimc->Draw("SAME");
	
	TLegend leg4(.75,.75,.9,.9); // (x1,y1,x2,y2)
	
	leg4.AddEntry("hphidata","Data");
	leg4.AddEntry("hphimc","MC");
	leg4.DrawClone("SAME");
	
}








void Zmass() {
	
	TFile *fdata = new TFile("Zmumu_data.root");
	TFile *fmc = new TFile("Zmumu_mc.root");
	
	fdata->cd("DecayTreeTuple");
	TTree *tdata = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	fmc->cd("DecayTreeTuple");
	TTree *tmc = dynamic_cast<TTree*>(gDirectory->Get("tuple"));
	
	// jet transverse momentum and pseudorapidity
	double jetpt_data = 0.0;
	double jeteta_data = 0.0;
	double jetpt_mc = 0.0;
	double jeteta_mc = 0.0;
	
	tdata->SetBranchAddress("Jet_PT",&jetpt_data);
	tdata->SetBranchAddress("Jet_Eta",&jeteta_data);
	tmc->SetBranchAddress("Jet_PT",&jetpt_mc);
	tmc->SetBranchAddress("Jet_Eta",&jeteta_mc);
	
	// Z mass
	double ZM_data = 0.0;
	double ZM_mc = 0.0;
	
	tdata->SetBranchAddress("Z0_MM",&ZM_data);
	tmc->SetBranchAddress("Z0_MM",&ZM_mc);
	
	// muon variables
	double muppt_data = 0.0;
	double muppx_data = 0.0;
	double muppy_data = 0.0;
	double muppz_data = 0.0;
	
	double mumpt_data = 0.0;
	double mumpx_data = 0.0;
	double mumpy_data = 0.0;
	double mumpz_data = 0.0;
	
	double muppt_mc = 0.0;
	double muppx_mc = 0.0;
	double muppy_mc = 0.0;
	double muppz_mc = 0.0;
	
	double mumpt_mc = 0.0;
	double mumpx_mc = 0.0;
	double mumpy_mc = 0.0;
	double mumpz_mc = 0.0;
	
	tdata->SetBranchAddress("mup_PT",&muppt_data);
	tdata->SetBranchAddress("mup_PX",&muppx_data);
	tdata->SetBranchAddress("mup_PY",&muppy_data);
	tdata->SetBranchAddress("mup_PZ",&muppz_data);
	
	tdata->SetBranchAddress("mum_PT",&mumpt_data);
	tdata->SetBranchAddress("mum_PX",&mumpx_data);
	tdata->SetBranchAddress("mum_PY",&mumpy_data);
	tdata->SetBranchAddress("mum_PZ",&mumpz_data);
	
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
	
	// to select only events with one jet
	ULong64_t totCand_data = 0;
	ULong64_t totCand_mc = 0;
	
	tdata->SetBranchAddress("totCandidates",&totCand_data);
	tmc->SetBranchAddress("totCandidates",&totCand_mc);
	
	TH1D *hZM_data = new TH1D("hZM_data","Calibration sample Z+jet - invariant mass;m_{#mu^{+}#mu^{-}} [GeV];A.U.",650,0,1300);
	TH1D *hZM_mc = new TH1D("hZM_mc","Calibration sample Z+jet - invariant mass;m_{#mu^{+}#mu^{-}} [GeV];A.U.",650,0,1300);
	
	for (int d=0; d<tdata->GetEntries()+1; d++) {
		
		if (d % 1000 == 0) cout << "Data event " << d << endl;
		
		tdata->GetEntry(d);
		
		double mupeta_data = -log(tan(acos(muppz_data/sqrt(muppx_data*muppx_data+muppy_data*muppy_data+muppz_data*muppz_data))/2));
		
		double mumeta_data = -log(tan(acos(mumpz_data/sqrt(mumpx_data*mumpx_data+mumpy_data*mumpy_data+mumpz_data*mumpz_data))/2));
		
		if (totCand_data != 1) continue; // select events with only one jet
		
		if (jetpt_data < 2e4) continue; // cut in pT of jet
		
		if (jeteta_data < 2.2 || jeteta_data > 4.2) continue; // cut in eta of jet
		
		if (ZM_data < 4000) continue; // neglect peak at 0
		
		if (muppt_data < 2e4 || mumpt_data < 2e4) continue; // cut in pT of muons
		
		if (mupeta_data < 2 || mupeta_data > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_data < 2 || mumeta_data > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_data == 1 && mupHLT1MuonHighPt_data == 1 && mupHLT2_data == 1) || (mumL0MuonEW_data == 1 && mumHLT1MuonHighPt_data == 1 && mumHLT2_data == 1)) { // cuts in trigger lines on muons
			
			hZM_data->Fill(ZM_data/1000);
			
		}
		
	}
	
	for (int m=0; m<tmc->GetEntries()+1; m++) {
		
		if (m % 1000 == 0) cout << "MC event " << m << endl;
		
		tmc->GetEntry(m);
			
		double mupeta_mc = -log(tan(acos(muppz_mc/sqrt(muppx_mc*muppx_mc+muppy_mc*muppy_mc+muppz_mc*muppz_mc))/2));
		
		double mumeta_mc = -log(tan(acos(mumpz_mc/sqrt(mumpx_mc*mumpx_mc+mumpy_mc*mumpy_mc+mumpz_mc*mumpz_mc))/2));
		
		if (totCand_mc != 1) continue; // select events with only one jet
		
		if (jetpt_mc < 2e4) continue; // cut in pT of jet
		
		if (jeteta_mc < 2.2 || jeteta_mc > 4.2) continue; // cut in eta of jet
		
		if (ZM_mc < 4000) continue; // neglect peak at 0
		
		if (muppt_mc < 2e4 || mumpt_mc < 2e4) continue; // cut in pT of muons
		
		if (mupeta_mc < 2 || mupeta_mc > 4.5) continue; // cut in eta of mu+
		
		if (mumeta_mc < 2 || mumeta_mc > 4.5) continue; // cut in eta of mu-
		
		if ((mupL0MuonEW_mc == 1 && mupHLT1MuonHighPt_mc == 1 && mupHLT2_mc == 1) || (mumL0MuonEW_mc == 1 && mumHLT1MuonHighPt_mc == 1 && mumHLT2_mc == 1)) { // cuts in trigger lines on muons
			
			hZM_mc->Fill(ZM_mc/1000);
			
		}
		
	}
	
	double Nsel_data = hZM_data->Integral();
	double Nsel_mc = hZM_mc->Integral();
		
	cout << "\nNsel data = " << Nsel_data << endl;
	cout << "\nNsel MC = " << Nsel_mc << "\n" << endl;
	
	hZM_data->Scale(1/(hZM_data->Integral()));
	hZM_mc->Scale(1/(hZM_mc->Integral()));
	
	TCanvas *cm = new TCanvas("cm","Z invariant mass");
	
	hZM_data->GetYaxis()->SetRangeUser(0,0.35);
	hZM_data->SetLineColor(kRed);
	hZM_data->SetMarkerStyle(20);
	hZM_data->SetMarkerSize(1);
	hZM_data->SetMarkerColor(kRed);
	hZM_data->GetXaxis()->SetNdivisions(504,kFALSE);
	hZM_mc->GetYaxis()->SetRangeUser(0,0.35);
	hZM_mc->SetLineColor(kBlue+2);
	hZM_mc->SetMarkerStyle(20);
	hZM_mc->SetMarkerSize(1);
	hZM_mc->SetMarkerColor(kBlue+2);
	
	TLine *xline = new TLine(91.1880,0,91.1880,0.35);
	
	xline->SetLineColor(kBlack);
	xline->SetLineWidth(3);
	
	gStyle->SetOptStat("");
	
	cm->cd();
	
	hZM_mc->Draw();
	hZM_data->Draw("SAME");
	xline->Draw("SAME");
	
	TLegend leg(.6,.7,.9,.9); // (x1,y1,x2,y2)
	
	leg.AddEntry(hZM_data,"Data","PL");
	leg.AddEntry(hZM_mc,"MC","PL");
	leg.AddEntry(xline,"Nominal Z^{0} mass","L");
	leg.DrawClone("SAME");
	
}
