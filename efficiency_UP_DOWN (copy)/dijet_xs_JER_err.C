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
	
	t_mc->SetBranchAddress("jet0_pt",&jet1pt_mc);
	t_mc->SetBranchAddress("jet1_pt",&jet2pt_mc);
	t_mc->SetBranchAddress("weight",&weight);
	t_mc->SetBranchAddress("dijet_dphi",&dijet_dphi_mc);
	t_data->SetBranchAddress("dijet_dphi",&dijet_dphi_data);
	
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
	double NselUP = 0.0;
	double NselDOWN = 0.0;
	double Ntot = 0.0;
	
	// new ROOT file with tree for corrected jet1pt and jet2pt
	TFile* newfile = new TFile("correctedJER.root","RECREATE");
	TTree* newtree = t_mc->CloneTree(-1); // -1 means I am copying all the entries
	
	double jet1pt_mc_corr = 0.0;
	double jet2pt_mc_corr = 0.0;
	double jet1ptUP_mc_corr = 0.0;
	double jet2ptUP_mc_corr = 0.0;
	double jet1ptDOWN_mc_corr = 0.0;
	double jet2ptDOWN_mc_corr = 0.0;
	
	newtree->Branch("jet1pt_mc_corr",&jet1pt_mc_corr,"jet1pt_mc_corr/D");
	newtree->Branch("jet2pt_mc_corr",&jet2pt_mc_corr,"jet2pt_mc_corr/D");
	newtree->Branch("jet1ptUP_mc_corr",&jet1ptUP_mc_corr,"jet1ptUP_mc_corr/D");
	newtree->Branch("jet2ptUP_mc_corr",&jet2ptUP_mc_corr,"jet2ptUP_mc_corr/D");
	newtree->Branch("jet1ptDOWN_mc_corr",&jet1ptDOWN_mc_corr,"jet1ptDOWN_mc_corr/D");
	newtree->Branch("jet2ptDOWN_mc_corr",&jet2ptDOWN_mc_corr,"jet2ptDOWN_mc_corr/D");
	
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
	
	// arrays to store JES, JER, errJER, JERup, JERdown values
	double JES[16] = {0.0};
	double JER[16] = {0.0};
	double errJER[16] = {0.0};
	double JERup[16] = {0.0};
	double JERdown[16] = {0.0};
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
	
	while (infile2 >> coljer1 >> coljer2 >> coljer3) {
		
		JER[jer] = coljer2;
		errJER[jer] = coljer3;
		jer++;
		
	}
	
	// random generator for smearing
	TRandom3 *RanGen = new TRandom3(0);
	
	for (int n=0; n<16; n++) {
		
		JERup[n] = JER[n] + errJER[n];
		JERdown[n] = JER[n] - errJER[n];
		
		cout << "\n    JES of bin " << n+1 << " = " << JES[n] << endl;
		cout << "    JER of bin " << n+1 << " = " << JER[n] << " ± " << errJER[n] << "\n  JERup of bin " << n+1 << " = " << JERup[n] << "\nJERdown of bin " << n+1 << " = " << JERdown[n] << endl;
		
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
			
			jet1ptUP_mc_corr = RanGen->Gaus((jet1pt_mc*JES[n]),JERup[n]*(jet1pt_mc*JES[n]));
			newtree->GetBranch("jet1ptUP_mc_corr")->Fill();
			
			jet2ptUP_mc_corr = RanGen->Gaus((jet2pt_mc*JES[n]),JERup[n]*(jet2pt_mc*JES[n]));
			newtree->GetBranch("jet2ptUP_mc_corr")->Fill();
			
			jet1ptDOWN_mc_corr = RanGen->Gaus((jet1pt_mc*JES[n]),JERdown[n]*(jet1pt_mc*JES[n]));
			newtree->GetBranch("jet1ptDOWN_mc_corr")->Fill();
			
			jet2ptDOWN_mc_corr = RanGen->Gaus((jet2pt_mc*JES[n]),JERdown[n]*(jet2pt_mc*JES[n]));
			newtree->GetBranch("jet2ptDOWN_mc_corr")->Fill();
			
		}
		
	}
	
	newfile->Write();
	
	for (int a=0; a<newtree->GetEntries(); a++){
		
		newtree->GetEntry(a);
		
		if (dijet_dphi_mc < 1) continue;
		
		if (jet1pt_mc_corr < 2e4 || jet2pt_mc_corr < 2e4 || jet1pt_mc_corr > 10e4 || jet2pt_mc_corr > 10e4) continue;
		
		Nsel += weight;
		
		if (jet1ptUP_mc_corr < 2e4 || jet2ptUP_mc_corr < 2e4 || jet1ptUP_mc_corr > 10e4 || jet2ptUP_mc_corr > 10e4) continue;
		
		NselUP += weight;
		
		if (jet1ptDOWN_mc_corr < 2e4 || jet2ptDOWN_mc_corr < 2e4 || jet1ptDOWN_mc_corr > 10e4 || jet2ptDOWN_mc_corr > 10e4) continue;
		
		NselDOWN += weight;
		
	}
	
	Ntot = 2e6*1+2e6*0.22+2e6*0.12+2e6*0.0035;
	
	double eff_ALL = Nsel/Ntot;
	double eff_JERup = NselUP/Ntot;
	double eff_JERdown = NselDOWN/Ntot;
	
	double Delta_eff_JERup = abs(eff_ALL-eff_JERup);
	double Delta_eff_JERdown = abs(eff_ALL-eff_JERdown);
	
	cout << "\n     Ntot = " << Ntot << "\n     Nsel = " << Nsel << "\n   NselUP = " << NselUP << "\n NselDOWN = " << NselDOWN << endl;
	cout << "\n     εALL = " << eff_ALL << "\n   εJERup = " << eff_JERup << "\n εJERdown = " << eff_JERdown << endl;
	cout << "\n  ΔεJERup = " << Delta_eff_JERup << "\nΔεJERdown = " << Delta_eff_JERdown << endl;
	
	// phase space correction
	
	// sigma(MadGraph phase space)/sigma(LHCb phase space) computed with MadGraph
	double sigma_lhcb = 8.6e8;       // pb, pt_hat>10 GeV and theta<400 mrad
	double sigma_mad  = 4.6643252e9; // pb, as=0.118, m>15 GeV and pt>10 GeV
	
	// acceptance from MadGraph for NumJets>=2 with as=0.118, jetpt>20 GeV, 2.2<jeteta<4.2
	double A = 0.002933; // 2933/1000000 events with 2 or more jets
	
	double ps = A*sigma_mad/sigma_lhcb;
	
	// luminosity (2% uncertainty)
	double lumi = 1600; // pb^-1
	
	// GEC efficiency
	double egec = 0.6;
	
	// pre-scale (HLT*stripping)
	double pre = 0.001*0.0013;
	
	// cross section
	double xs_ALL = Ndata/(eff_ALL*egec*lumi*pre*ps);
	
	cout << "\nsigma_lhcb = " << sigma_lhcb << " pb" << "\n\nA*sigma_mad (to be compared with sigma_lhcb) = " << A*sigma_mad << " pb" << "\n\nεALL = " << eff_ALL << "\n\nPhase space correction: " << ps << "\n\nσjj with JES+JER correction = " << xs_ALL << "\n" << endl;
	
//	" ± " << ER_xs_ALL << " pb\n" << endl;
	
	return;
	
}
