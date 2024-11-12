// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file rootIOTree_example_write.cc
 *  @brief Basic example of use of root I/O: writing events to file
 *
 *  @author Witold Pokorski/Andrii Verbytskyi
 *  @date   29/10/15
 */
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/WriterRootTree.h"
#include "HepMC3/Print.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"


#include <iostream>

using namespace HepMC3;

/** Main */
int main(int argc, char **argv) {

    if( argc<3 ) {
        std::cout << "Usage: " << argv[0] << " <input_hepmc2_file> <output_root_file>" << std::endl;
        exit(-1);
    }

    ReaderAsciiHepMC2 text_input(argv[1]);
    WriterRootTree  root_output(argv[2]);

    std::vector <fastjet::PseudoJet> fjInputs, inclusiveJets, sortedJets;
    fastjet::JetDefinition *jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.5); //, fastjet::E_scheme, fastjet::Best);

    int events_parsed = 0;
    FourVector *jetmom;
    GenParticle *jet;

    while( !text_input.failed() ) {

        GenEvent evt(Units::GEV,Units::MM);

        text_input.read_event(evt);
        if( text_input.failed() ) break;

        if( events_parsed == 0 ) {
            std::cout << "First event: " << std::endl;
        }

	fjInputs.resize(0);
	inclusiveJets.resize(0);
	sortedJets.resize(0);

	const std::vector<GenParticlePtr> particles = evt.particles();
	for(unsigned long ip = 0; ip < particles.size(); ip++){
		unsigned isfinal = 0, ishadron = 0, isMuon=0, isPhoton=0, isElectron=0;
		GenParticle p = *particles[ip];
		std::vector<GenParticlePtr> children = p.children();
		if(children.size() == 0) isfinal++;
		if(abs(p.pid()) > 60) ishadron++;
        if(abs(p.pid()) == 13) isMuon++;
        if(abs(p.pid()) == 11) isElectron++;
        if(abs(p.pid()) == 22) isPhoton++;
		if( isfinal && (ishadron || isMuon || isElectron || isPhoton) ){
			fjInputs.push_back(fastjet::PseudoJet(p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.momentum().e()));
		}
	}

	// Run Fastjet algorithm
	// Extract inclusive jets sorted by pT (min pT of 5.0 GeV)
	fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
	inclusiveJets = clustSeq.inclusive_jets(5.0);
	sortedJets = sorted_by_pt(inclusiveJets);

	// for(unsigned i = 0; i < sortedJets.size();i++){
	// 	jetmom = new FourVector(sortedJets[i].px(), sortedJets[i].py(), sortedJets[i].pz(), sortedJets[i].e());
	// 	jet = new GenParticle(*jetmom, 999, 1);
	// 	evt.add_particle(jet);
	// }


	// // Write
    //     root_output.write_event(evt);
    //     ++events_parsed;

    // Create a new event for storing only jets
    GenEvent jetEvent(Units::GEV, Units::MM);

    for (unsigned i = 0; i < sortedJets.size(); ++i) {
        // Only add jets to the jetEvent
        FourVector jetmom(sortedJets[i].px(), sortedJets[i].py(), sortedJets[i].pz(), sortedJets[i].e());
        GenParticlePtr jet = std::make_shared<GenParticle>(jetmom, 999, 1);
        jetEvent.add_particle(jet);
    }

    // Write only the jetEvent to the ROOT file
    root_output.write_event(jetEvent);
    ++events_parsed;

        if((events_parsed % 1000) == 0) std::cout << "Event: " << events_parsed << std::endl;
	
    }






    

    text_input.close();
    root_output.close();

    std::cout << "Events parsed and written: " << events_parsed << std::endl;

    delete jetDef;

    return 0;
}
