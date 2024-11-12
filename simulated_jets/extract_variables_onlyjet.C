#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

using namespace std;

void extract_variables() {
	
	// create and initialize arrays to store variables
	double px[100] = {0.0};
	double py[100] = {0.0};
	double pz[100] = {0.0};
	double E[100] = {0.0};
	
	// create addresses to branches
	TBranch *b_px;
	TBranch *b_py;
	TBranch *b_pz;
	TBranch *b_E;
	
	// open the input ROOT file
	TFile *fin = new TFile("pp_jj_ct10_05gev_0130.root");
	
	// get the tree from the input ROOT file
	TTree *tin = (TTree*)fin->Get("hepmc3_tree");
	
	int N = tin->GetEntries();
	
	cout << "\nNumber of entries: " << N << endl;
	
	// set branch addresses from input tree
	tin->SetBranchAddress("particles.momentum.m_v1",px,&b_px);
	tin->SetBranchAddress("particles.momentum.m_v2",py,&b_py);
	tin->SetBranchAddress("particles.momentum.m_v3",pz,&b_pz);
	tin->SetBranchAddress("particles.momentum.m_v4",E,&b_E);
	
	// create new output file
	TFile *fout = new TFile("jets_ct10_05gev_0130.root","RECREATE");
	
	// create new output tree
	TTree *tout = new TTree("jets_ct10_05gev_0130","jets_ct10_05gev_0130");
	
	// create output branches
	double jet1pt = 0.0;
	double jet1eta = 0.0;
	double jet1phi = 0.0;
	double jet1px = 0.0;
	double jet1py = 0.0;
	double jet1pz = 0.0;
	double jet1E = 0.0;
	double jet1mass = 0.0;
	
	double jet2pt = 0.0;
	double jet2eta = 0.0;
	double jet2phi = 0.0;
	double jet2px = 0.0;
	double jet2py = 0.0;
	double jet2pz = 0.0;
	double jet2E = 0.0;
	double jet2mass = 0.0;
	
	double jet3pt = 0.0;
	double jet3eta = 0.0;
	double jet3phi = 0.0;
	double jet3px = 0.0;
	double jet3py = 0.0;
	double jet3pz = 0.0;
	double jet3E = 0.0;
	double jet3mass = 0.0;
	
	double Q = 0.0;
	
	tout->Branch("jet1pt",&jet1pt,"jet1pt/D");
	tout->Branch("jet1eta",&jet1eta,"jet1eta/D");
	tout->Branch("jet1phi",&jet1phi,"jet1phi/D");
	tout->Branch("jet1px",&jet1px,"jet1px/D");
	tout->Branch("jet1py",&jet1py,"jet1py/D");
	tout->Branch("jet1pz",&jet1pz,"jet1pz/D");
	tout->Branch("jet1E",&jet1E,"jet1E/D");
	tout->Branch("jet1mass",&jet1mass,"jet1mass/D");
	
	tout->Branch("jet2pt",&jet2pt,"jet2pt/D");
	tout->Branch("jet2eta",&jet2eta,"jet2eta/D");
	tout->Branch("jet2phi",&jet2phi,"jet2phi/D");
	tout->Branch("jet2px",&jet2px,"jet2px/D");
	tout->Branch("jet2py",&jet2py,"jet2py/D");
	tout->Branch("jet2pz",&jet2pz,"jet2pz/D");
	tout->Branch("jet2E",&jet2E,"jet2E/D");
	tout->Branch("jet2mass",&jet2mass,"jet2mass/D");
	
	tout->Branch("jet3pt",&jet3pt,"jet3pt/D");
	tout->Branch("jet3eta",&jet3eta,"jet3eta/D");
	tout->Branch("jet3phi",&jet3phi,"jet3phi/D");
	tout->Branch("jet3px",&jet3px,"jet3px/D");
	tout->Branch("jet3py",&jet3py,"jet3py/D");
	tout->Branch("jet3pz",&jet3pz,"jet3pz/D");
	tout->Branch("jet3E",&jet3E,"jet3E/D");
	tout->Branch("jet3mass",&jet3mass,"jet3mass/D");
	
	tout->Branch("Q",&Q,"Q/D");
	
	// loop over events
	for (int ev=0; ev<N; ev++) {
		
		// if event_number/100 gives no remainder -> display multiples of 100
		if (ev%1==0) {
			cout << "\nEVENT NUMBER: " << ev << endl;
		}
		
		// reset arrays for each event to avoid data carryover
		memset(px,0.0,sizeof(px));
		memset(py,0.0,sizeof(py));
		memset(pz,0.0,sizeof(pz));
		memset(E,0.0,sizeof(E));
		
		// access input tree
		tin->GetEntry(ev);
		
		// flag to determine if the event has valid jets
		bool ValidJets = false;
		
		// counter to track the number of found jets
		int jet_count = 0;
		
		// loop over particles per event
		for (int part=0; part<100; part++) {
			
			double eta = -log(tan(acos(pz[part]/sqrt(px[part]*px[part]+py[part]*py[part]+pz[part]*pz[part]))/2));
			
			if (isnan(eta)) continue;
			
			cout << "Eta: " << eta << endl;
			
			// reject particles outside acceptance 2.2 < eta < 4.2
			if (eta < 2.2 || eta > 4.2) continue;
			
			// compute jet variables
			jet_count++;
				
			// first 3 jets are already sorted in pT for each event
			cout << "Jet number: " << jet_count << endl;
			
			if (jet_count == 1) {
				cout << "First jet is: " << part << endl;
				
				// first jet
				jet1pt = sqrt(px[part]*px[part]+py[part]*py[part]);
				cout << "jet1pt = " << jet1pt << endl;
				jet1eta = -log(tan(acos(pz[part]/sqrt(px[part]*px[part]+py[part]*py[part]+pz[part]*pz[part]))/2));
				jet1phi = acos(px[part]/(sqrt(px[part]*px[part]+py[part]*py[part])));
				jet1px = px[part];
				jet1py = py[part];
				jet1pz = pz[part];
				jet1E = E[part];
				jet1mass = sqrt(E[part]*E[part]-px[part]*px[part]-py[part]*py[part]-pz[part]*pz[part]);
			}
			
			else if (jet_count == 2) {
				cout << "Second jet is: " << part << endl;
				
				// second jet
				jet2pt = sqrt(px[part]*px[part]+py[part]*py[part]);
				cout << "jet2pt = " << jet2pt << endl;
				jet2eta = -log(tan(acos(pz[part]/sqrt(px[part]*px[part]+py[part]*py[part]+pz[part]*pz[part]))/2));
				jet2phi = acos(px[part]/(sqrt(px[part]*px[part]+py[part]*py[part])));
				jet2px = px[part];
				jet2py = py[part];
				jet2pz = pz[part];
				jet2E = E[part];
				jet2mass = sqrt(E[part]*E[part]-px[part]*px[part]-py[part]*py[part]-pz[part]*pz[part]);
				Q = (jet1pt+jet2pt)/2;
			}
			
			else if (jet_count == 3) {
				cout << "Third jet is: " << part << endl;
				
				// third jet
				jet3pt = sqrt(px[part]*px[part]+py[part]*py[part]);
				cout << "jet3pt = " << jet3pt << endl;
				jet3eta = -log(tan(acos(pz[part]/sqrt(px[part]*px[part]+py[part]*py[part]+pz[part]*pz[part]))/2));
				jet3phi = acos(px[part]/(sqrt(px[part]*px[part]+py[part]*py[part])));
				jet3px = px[part];
				jet3py = py[part];
				jet3pz = pz[part];
				jet3E = E[part];
				jet3mass = sqrt(E[part]*E[part]-px[part]*px[part]-py[part]*py[part]-pz[part]*pz[part]);
				
				break; // found three jets, exit the loop
			}
			
			// set ValidJets flag when at least one jet is found
			ValidJets = true;
		}
		
		// check if the event has valid jets before filling the tree
		if (ValidJets) {
			// fill output tree
			tout->Fill();
		}
		
		else {
			// reset jet variables only if the event has no valid jets
			// this ensures no overwriting of previous event variables
			jet1pt = 0.0;
			jet1eta = 0.0;
			jet1phi = 0.0;
			jet1px = 0.0;
			jet1py = 0.0;
			jet1pz = 0.0;
			jet1E = 0.0;
			jet1mass = 0.0;
			
			jet2pt = 0.0;
			jet2eta = 0.0;
			jet2phi = 0.0;
			jet2px = 0.0;
			jet2py = 0.0;
			jet2pz = 0.0;
			jet2E = 0.0;
			jet2mass = 0.0;
			Q = 0.0;
			
			jet3pt = 0.0;
			jet3eta = 0.0;
			jet3phi = 0.0;
			jet3px = 0.0;
			jet3py = 0.0;
			jet3pz = 0.0;
			jet3E = 0.0;
			jet3mass = 0.0;
		}
	}
		
	fout->Write();
	fout->Close();
	fin->Close();
	
}
