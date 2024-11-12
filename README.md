## "Measurement of the strong coupling constant with the LHCb detector"

*Alice Moro's masters thesis @ unipd*

All useful ROOT files can be found here: https://drive.google.com/drive/folders/1TAwn5ll-_hyRoDiqyBTvuqv1rN44elUB?usp=sharing

Here a description of the content of the folders is presented. To better understand the analysis flow, they should be employed in the order proposed in the following.

1. **simulated_jets**: extraction of jet variables (e.g. from pp_jj_ct10_10gev_0118.root - that is the output of ana.cc, which converts HepMC files to ROOT files - to jets_ct10_10gev_0118.root), computating and plotting distributions (p<sub>T</sub>, Q, mass, eta, phi, R<sub>32</sub>, etc.), acceptances and minimum number of jets per event, used to get the simulated weighted cross sections to be fitted and compared to the experimental one in order to measure α<sub>s</sub>;

2. **calibration_samples_K_F**: studying R<sub>Z</sub> and R<sub>12</sub> ratios with Z+jet and dijet calibration samples respectively, extracting K* (JES) and F* (JER) calibration factors;

3. **cross_sec_alpha_s_extraction**: correcting jets p<sub>T</sub> with K* and F* factors, determining efficiencies relative to K* and F* factors, measuring data dijet cross section, obtaining α<sub>s</sub> estimate.
