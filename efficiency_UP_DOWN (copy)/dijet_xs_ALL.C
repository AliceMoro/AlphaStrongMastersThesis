#include <math.h>

using namespace std;

void dijet_xs() {
	
	TFile* f_data = new TFile("/home/alice/Desktop/TESI/calibration_samples/dijet/more_F/dijet_data_2016.root");
	TTree* t_data = (TTree*) f_data->Get("Tuple");
	
	TFile* f_mc = new TFile("/home/alice/Desktop/TESI/calibration_samples/dijet/more_F/qq.root");
	TTree* t_mc = (TTree*) f_mc->Get("Tuple");
	
	double jet1pt_mc = 0.0;
	double jet2pt_mc = 0.0;
	double weight = 0.0;
	double dijet_dphi_mc = 0.0;
	double dijet_dphi_data = 0.0;
	double jet1pt_data = 0.0;
	double jet2pt_data = 0.0;
	
	t_mc->SetBranchAddress("jet0_pt",&jet1pt_mc);
	t_mc->SetBranchAddress("jet1_pt",&jet2pt_mc);
	t_mc->SetBranchAddress("weight",&weight);
	t_mc->SetBranchAddress("dijet_dphi",&dijet_dphi_mc);
	t_data->SetBranchAddress("dijet_dphi",&dijet_dphi_data);
	t_data->SetBranchAddress("jet0_pt",&jet1pt_data);
	t_data->SetBranchAddress("jet1_pt",&jet2pt_data);
	
	// number of events
	int NdataOrig = t_data->GetEntries();
	int NmcOrig = t_mc->GetEntries();
	
	cout << "\nNdata = " << NdataOrig << "\nNmc   = " << NmcOrig << endl;
	
	// number of events after Δϕ cut
	int Ndata = t_data->GetEntries("dijet_dphi>1");
	int Nmc = t_mc->GetEntries("dijet_dphi>1");
	
	cout << "\nNdata(Δϕ>1) = " << Ndata << "\nNmc(Δϕ>1)   = " << Nmc << "\n" << endl;
	
	// efficiency
	
	double Nsel = 0.0;
	double Ntot = 0.0;
	
	// new ROOT file with tree for corrected jet1pt and jet2pt
	TFile* newfile = new TFile("correctedALL.root","RECREATE");
	TTree* newtree = t_mc->CloneTree(-1); // -1 means I am copying all the entries
	
	double jet1pt_mc_corr = 0.0;
	double jet2pt_mc_corr = 0.0;
	
	newtree->Branch("jet1pt_mc_corr",&jet1pt_mc_corr,"jet1pt_mc_corr/D");
	newtree->Branch("jet2pt_mc_corr",&jet2pt_mc_corr,"jet2pt_mc_corr/D");
	
	// array to store jetpt bins minima
	double jetpt_min[16] = {0.0};
	jetpt_min[0] = 2e4;
	
	for (int min=0; min<16; min++) jetpt_min[min+1] = jetpt_min[min] + 5e3;
	
	// array to store jetpt bins maxima
	double jetpt_max[16] = {0.0};
	jetpt_max[0] = 2.5e4;
	
	for (int max=0; max<16; max++) {
		
		jetpt_max[max+1] = jetpt_max[max] + 5e3;
		
		cout << "Bin number " << max+1 << ": [" << jetpt_min[max]/1000 << "," << jetpt_max[max]/1000 << "] GeV" << endl;
		
	}
	
	// arrays to store JES and JER values
	double JES[16] = {0.0};
	double JER[16] = {0.0};
	double coljes1 = 0.0;
	double coljes2 = 0.0;
	double coljes3 = 0.0;
	double coljer1 = 0.0;
	double coljer2 = 0.0;
	double coljer3 = 0.0;
	
	ifstream infile1("/home/alice/Desktop/TESI/calibration_samples/Z+jet/Zmumu/JES_more_bins_fit_wo_errors/JESvsJetpT.txt");
	ifstream infile2("/home/alice/Desktop/TESI/calibration_samples/dijet/more_F/JERvsJetpT.txt");
	
	int jes = 0;
	int jer = 0;
	
	for (int j=0; j<2; j++) { // skip first 2 JES lines (dijet sample starts at 20 GeV)
		
		string line;
		getline(infile1,line);
		
	}
	
	while (infile1 >> coljes1 >> coljes2 >> coljes3) JES[jes++] = coljes2;
	
	while (infile2 >> coljer1 >> coljer2 >> coljer3) JER[jer++] = coljer2;
	
	// random generator for smearing
	TRandom3 *RanGen = new TRandom3(0);
	
	for (int n=0; n<16; n++) {
		
		cout << "\nJES of bin " << n+1 << " = " << JES[n] << endl;
		cout << "JER of bin " << n+1 << " = " << JER[n] << endl;
		
		for (int i=0; i<t_mc->GetEntries(); i++) {
		
			t_mc->GetEntry(i);
			
			if (jet1pt_mc < jetpt_min[n] || jet2pt_mc < jetpt_min[n] || jet1pt_mc > jetpt_max[n] || jet2pt_mc > jetpt_max[n]) {
				
				jet1pt_mc_corr = jet1pt_mc;
				
				jet2pt_mc_corr = jet2pt_mc;
				
			}
			
			
			jet1pt_mc_corr = RanGen->Gaus((jet1pt_mc*JES[n]),JER[n]*(jet1pt_mc*JES[n]));
			newtree->GetBranch("jet1pt_mc_corr")->Fill();
			
			jet2pt_mc_corr = RanGen->Gaus((jet2pt_mc*JES[n]),JER[n]*(jet2pt_mc*JES[n]));
			newtree->GetBranch("jet2pt_mc_corr")->Fill();
			
		}
		
	}
	
	newfile->Write();
	
	for (int a=0; a<newtree->GetEntries(); a++){
		
		newtree->GetEntry(a);
		
		if (dijet_dphi_mc < 1) continue;
		
		if (jet1pt_mc_corr < 2e4 || jet2pt_mc_corr < 2e4 || jet1pt_mc_corr > 10e4 || jet2pt_mc_corr > 10e4) continue;
		
		Nsel += weight;
		
	}
	
	Ntot = 2e6*1+2e6*0.22+2e6*0.12+2e6*0.0035;
	
	double eff_ALL = Nsel/Ntot;
	
	cout << "\nNtot = " << Ntot << "\nNsel = " << Nsel << endl;
	
	
	
	
	
	
	
	
	
	
	
		
	// plot MC dijet sample uncorrected and corrected
	
	TH1D* pt1 = new TH1D("pt1","p_{T}(jet_{1}) corr. with JES+JER - MC dijet sample;p_{T}(jet_{1}) [GeV];A.U.",32,0,160);
	TH1D* pt2 = new TH1D("pt2","p_{T}(jet_{2}) corr. with JES+JER - MC dijet sample;p_{T}(jet_{2}) [GeV];A.U.",32,0,160);
	TH1D* pt1corr = new TH1D("pt1corr","p_{T}(jet_{1}) corr. with JES+JER - MC dijet sample;p_{T}(jet_{1}) [GeV];A.U.",32,0,160);
	TH1D* pt2corr = new TH1D("pt2corr","p_{T}(jet_{2}) corr. with JES+JER - MC dijet sample;p_{T}(jet_{2}) [GeV];A.U.",32,0,160);
	
	for (int r=0; r<t_mc->GetEntries(); r++) {
		t_mc->GetEntry(r);
		if (jet1pt_mc >= 2e4) pt1->Fill(jet1pt_mc/1000,weight);
		if (jet2pt_mc >= 2e4) pt2->Fill(jet2pt_mc/1000,weight);
	}
	
	for (int s=0; s<newtree->GetEntries(); s++) {
		newtree->GetEntry(s);
		if (jet1pt_mc_corr >= 2e4) pt1corr->Fill(jet1pt_mc_corr/1000,weight);
		if (jet2pt_mc_corr >= 2e4) pt2corr->Fill(jet2pt_mc_corr/1000,weight);
	}
	
	pt1->Scale(1/(pt1->Integral()));
	pt2->Scale(1/(pt2->Integral()));
	pt1corr->Scale(1/(pt1corr->Integral()));
	pt2corr->Scale(1/(pt2corr->Integral()));
	
	TCanvas *c = new TCanvas("c","pt w/o and w/ corr. JES");
	c->Divide(2,1);
	
	pt1->SetTitle("p_{T}(jet_{1}) corr. with JES+JER - MC dijet sample;p_{T}(jet_{1}) [GeV];A.U.");
	pt2->SetTitle("p_{T}(jet_{2}) corr. with JES+JER - MC dijet sample;p_{T}(jet_{2}) [GeV];A.U.");
	
	pt1->SetLineColor(kBlue+2);
	pt1->SetLineWidth(1);
	pt1->SetMarkerStyle(20);
	pt1->SetMarkerSize(1);
	pt1->SetMarkerColor(kBlue+2);
	
	pt2->SetLineColor(kBlue+2);
	pt2->SetLineWidth(1);
	pt2->SetMarkerStyle(20);
	pt2->SetMarkerSize(1);
	pt2->SetMarkerColor(kBlue+2);
	
	pt1corr->SetLineColor(kViolet+6);
	pt1corr->SetLineWidth(1);
	pt1corr->SetMarkerStyle(20);
	pt1corr->SetMarkerSize(1);
	pt1corr->SetMarkerColor(kViolet+6);
	
	pt2corr->SetLineColor(kViolet+6);
	pt2corr->SetLineWidth(1);
	pt2corr->SetMarkerStyle(20);
	pt2corr->SetMarkerSize(1);
	pt2corr->SetMarkerColor(kViolet+6);
	
	pt1->GetYaxis()->SetRangeUser(0,0.43);
	pt2->GetYaxis()->SetRangeUser(0,0.43);
	
	gStyle->SetOptStat("");
	
	c->cd(1);
	
	pt1->Draw("");
	pt1corr->Draw("SAME");
	
	TLegend *leg1 = new TLegend(0.5, 0.65, 0.9, 0.9);
	leg1->SetFillColor(0);
	leg1->SetBorderSize(1);
	leg1->AddEntry((TObject*)0,"JES from Z+jet sample","");
	leg1->AddEntry((TObject*)0,"JER from dijet sample","");
	leg1->AddEntry(pt1,"Uncorrected","PL");
	leg1->AddEntry(pt1corr,"Corrected","PL");
	leg1->Draw();
	
	c->cd(2);
	
	pt2->Draw("");
	pt2corr->Draw("SAME");
	
	TLegend *leg2 = new TLegend(0.5, 0.65, 0.9, 0.9);
	leg2->SetFillColor(0);
	leg2->SetBorderSize(1);
	leg2->AddEntry((TObject*)0,"JES from Z+jet sample","");
	leg2->AddEntry((TObject*)0,"JER from dijet sample","");
	leg2->AddEntry(pt2,"Uncorrected","PL");
	leg2->AddEntry(pt2corr,"Corrected","PL");
	leg2->Draw();
	
	
	
	
	
	
	
	
	
	
	
	
	// plot data and corrected MC from dijet sample
	
	TH1D* hpt1 = new TH1D("hpt1","p_{T}(jet_{1}) - data dijet sample;p_{T}(jet_{1}) [GeV];A.U.",32,0,160);
	TH1D* hpt2 = new TH1D("hpt2","p_{T}(jet_{2}) - data dijet sample;p_{T}(jet_{2}) [GeV];A.U.",32,0,160);
	TH1D* hpt1corr = new TH1D("hpt1corr","p_{T}(jet_{1}) corr. with JES+JER - MC dijet sample;p_{T}(jet_{1}) [GeV];A.U.",32,0,160);
	TH1D* hpt2corr = new TH1D("hpt2corr","p_{T}(jet_{2}) corr. with JES+JER - MC dijet sample;p_{T}(jet_{2}) [GeV];A.U.",32,0,160);
	
	for (int r=0; r<t_data->GetEntries(); r++) {
		t_data->GetEntry(r);
		if (jet1pt_data >= 2e4) hpt1->Fill(jet1pt_data/1000);
		if (jet2pt_data >= 2e4) hpt2->Fill(jet2pt_data/1000);
	}
	
	for (int s=0; s<newtree->GetEntries(); s++) {
		newtree->GetEntry(s);
		if (jet1pt_mc_corr >= 2e4) hpt1corr->Fill(jet1pt_mc_corr/1000,weight);
		if (jet2pt_mc_corr >= 2e4) hpt2corr->Fill(jet2pt_mc_corr/1000,weight);
	}
	
	hpt1->Scale(1/(hpt1->Integral()));
	hpt2->Scale(1/(hpt2->Integral()));
	hpt1corr->Scale(1/(hpt1corr->Integral()));
	hpt2corr->Scale(1/(hpt2corr->Integral()));
	
	TCanvas *c1 = new TCanvas("c1","pt data and corrected MC");
	c1->Divide(2,1);
	
	hpt1->SetTitle("p_{T}(jet_{1}) - data and corr. MC dijet sample;p_{T}(jet_{1}) [GeV];A.U.");
	hpt2->SetTitle("p_{T}(jet_{2}) - data and corr. MC dijet sample;p_{T}(jet_{2}) [GeV];A.U.");
	
	hpt1->SetLineColor(kRed);
	hpt1->SetLineWidth(1);
	hpt1->SetMarkerStyle(20);
	hpt1->SetMarkerSize(1);
	hpt1->SetMarkerColor(kRed);
	
	hpt2->SetLineColor(kRed);
	hpt2->SetLineWidth(1);
	hpt2->SetMarkerStyle(20);
	hpt2->SetMarkerSize(1);
	hpt2->SetMarkerColor(kRed);
	
	hpt1corr->SetLineColor(kViolet+6);
	hpt1corr->SetLineWidth(1);
	hpt1corr->SetMarkerStyle(20);
	hpt1corr->SetMarkerSize(1);
	hpt1corr->SetMarkerColor(kViolet+6);
	
	hpt2corr->SetLineColor(kViolet+6);
	hpt2corr->SetLineWidth(1);
	hpt2corr->SetMarkerStyle(20);
	hpt2corr->SetMarkerSize(1);
	hpt2corr->SetMarkerColor(kViolet+6);
	
	hpt1->GetYaxis()->SetRangeUser(0,0.43);
	hpt2->GetYaxis()->SetRangeUser(0,0.43);
	
	gStyle->SetOptStat("");
	
	c1->cd(1);
	
	hpt1->Draw("");
	hpt1corr->Draw("SAME");
	
	TLegend *leg3 = new TLegend(0.5, 0.65, 0.9, 0.9);
	leg3->SetFillColor(0);
	leg3->SetBorderSize(1);
	leg3->AddEntry((TObject*)0,"JES from Z+jet sample","");
	leg3->AddEntry((TObject*)0,"JER from dijet sample","");
	leg3->AddEntry(hpt1,"Data","PL");
	leg3->AddEntry(hpt1corr,"Corrected MC","PL");
	leg3->Draw();
	
	c1->cd(2);
	
	hpt2->Draw("");
	hpt2corr->Draw("SAME");
	
	TLegend *leg4 = new TLegend(0.5, 0.65, 0.9, 0.9);
	leg4->SetFillColor(0);
	leg4->SetBorderSize(1);
	leg4->AddEntry((TObject*)0,"JES from Z+jet sample","");
	leg4->AddEntry((TObject*)0,"JER from dijet sample","");
	leg4->AddEntry(hpt2,"Data","PL");
	leg4->AddEntry(hpt2corr,"Corrected MC","PL");
	leg4->Draw();
	
	
	
	
	
	
	
	
	
	
	
	
	// phase space correction
	
	// sigma(MadGraph phase space)/sigma(LHCb phase space) computed with MadGraph
	double sigma_lhcb = 8.6e8;       // pb, pt_hat>10 GeV and theta<400 mrad
	double sigma_mad  = 4.6643252e9; // pb, as=0.118, m>15 GeV and pt>10 GeV
	
	// acceptance from MadGraph for NumJets>=2 with as=0.118, jetpt>20 GeV, 2.2<jeteta<4.2
	double A = 0.002933; // 2933/1000000 events with 2 or more jets
	
	double ps = A*sigma_mad/sigma_lhcb; // *((2*M_PI)/(M_PI-2.8));
	
	// luminosity (2% uncertainty)
	double lumi = 1600; // pb^-1
	
	// GEC efficiency
	double egec = 0.6;
	
	// pre-scale (HLT*stripping)
	double pre = 0.001*0.013;
	
	// errors of efficiencies
	double Delta_eff_JESup = 4.7869e-05;
	double Delta_eff_JESdown = 0.000183574;
	double Delta_eff_JERup = 0.0001128;
	double Delta_eff_JERdown = 0.000151388;
	
	double Delta_eff_up = sqrt((Delta_eff_JESup*Delta_eff_JESup)+(Delta_eff_JERup*Delta_eff_JERup));
	double Delta_eff_down = sqrt((Delta_eff_JESdown*Delta_eff_JESdown)+(Delta_eff_JERdown*Delta_eff_JERdown));
	
	double Delta_eff = (Delta_eff_up+Delta_eff_down) / 2.0;
	
	// relative errors
	double ER_stat = 1/sqrt(Ndata);
	
	double ER_eff = Delta_eff/eff_ALL;
	
	double ER_L = 0.02;
	
	// cross section
	double xs_ALL = Ndata/(eff_ALL*egec*lumi*pre)*ps;
	
	// errors of cross section
	double err_stat = ER_stat*xs_ALL;
	
	double err_eff = ER_eff*xs_ALL;
	
	double err_L = ER_L*xs_ALL;
	
	double err_ALL = sqrt((err_stat*err_stat)+(err_eff*err_eff)+(err_L*err_L));
	
	cout << "\n  εALL = " << eff_ALL << "\n ΔεALL = " << Delta_eff << "\n  ΔεUP = " << Delta_eff_up << "\nΔεDOWN = " << Delta_eff_down << "\n\nPhase space correction: " << ps << "\n\nσjj with JES+JER correction = " << xs_ALL << " ± " << err_ALL << " pb\n" << "\nerr_stat = " << err_stat << "\nerr_eff = " << err_eff << "\nerr_L = " << err_L << endl;
	
	return;
	
}
