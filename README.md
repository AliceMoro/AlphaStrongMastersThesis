# Alice Moro's Masters Thesis work

All useful ROOT files can be found here: https://drive.google.com/drive/folders/1TAwn5ll-_hyRoDiqyBTvuqv1rN44elUB?usp=sharing

Here is a description of the content of the folders. To better understand the analysis flow, they should be employed in the order proposed in the following.

1. simulated_jets: extraction of jet variables (e.g. from pp_jj_ct10_10gev_0118.root - that is the output of ana.cc, which converts HepMC files to ROOT files - to jets_ct10_10gev_0118.root), computating and plotting distributions (pT, Q, mass, eta, phi, R_32, etc.), acceptances and minimum number of jets per event, used to get the simulated weighted cross sections to be fitted and compared to the experimental one in order to measure alpha strong;
2. calibration_samples_K*_F*: studying R_Z and R_12 ratios with dijet and Z+jet calibration samples, extracting K* (JES) and F*(JER) calibration factors;
3. cross_sec_alpha_s_extraction: correcting jets pT with K* and F* factors, determining efficiencies relative to K* and F* factors, measuring data dijet cross section, obtaining alpha strong estimate.
