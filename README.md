# A dataset and a matlab toolbox for morphological analysis of Fetal Heart Rate signal

All information on the wiki https://github.com/utsb-fmm/FHRMA/wiki

To start, type **FHRMA-menu** on Matlab or launch the app

![Presentation image](https://ercf.univ-catholille.fr/fhrviewer-demo/Presentation.png)
=======
FIG. 1. Illustration of the main *fhrmorpho* interface. (A) main window, (B) discordance evaluation between WMFB method and expert consensus (C) WMFB method analysis (D) Expert consensus analysis.

**This toolbox is related to several papers. Please cite those papers if you use any of the data or source code of this repository.**  [4] must be cited if you use the toolbox. [1]  must be cited if you use the morphological analysis (baseline, Acceleration, deceleration) [3] must be cited if you use the morphological analysis dataset. [5] must be cited if you use the WMFB method (current best) for morphological analysis. [6] must be cited if you use the false signal detection, method, interface and/or dataset.

[1] Houzé de l’Aulnoit, A., Boudet, S., Demailly, R., Delgranche, A., Peyrodie, L., Beuscart, R., Houzé de l’Aulnoit,D. - Automated fetal heart rate analysis for baseline determination and acceleration/deceleration detection: A comparison of 11 methods versus expert consensus. Biomedical Signal Processing and Control 49:113 -123,2019, DOI:10.1016/j.bspc.2018.10.002 [Download on Researchgate](https://www.researchgate.net/publication/329718625_Automated_fetal_heart_rate_analysis_for_baseline_determination_and_accelerationdeceleration_detection_A_comparison_of_11_methods_versus_expert_consensus)

[2] Houzé de l'Aulnoit, Agathe, Boudet, Samuel, Demailly, Romain, Peyrodie, Laurent, Beuscart, Regis, Houzé de l'Aulnoit, Denis - Baseline fetal heart rate analysis: eleven automatic methods versus expert consensus. Engineering in Medicine and Biology Society (EMBC), 2016 IEEE 38th Annual International Conference of the pp. 3576--3581,2016, DOI:10.1109/EMBC.2016.7591501 [Download on researchgate](https://www.researchgate.net/publication/309349819_Baseline_fetal_heart_rate_analysis_Eleven_automatic_methods_versus_expert_consensus)

[3] Boudet, S., Houzé de l’Aulnoit, A., Demailly, R., Delgranche, A., Peyrodie, L., Beuscart, R., Houzé de l’Aulnoit,D. - Fetal heart rate signal dataset for training morphological analysis methods and evaluating them against an expert consensus. Preprints pp. Submitted to data in brief,2019, DOI:10.20944/preprints201907.0039.v1 [Download on researchgate](https://www.researchgate.net/publication/334164380_Fetal_Heart_Rate_Signal_Dataset_for_Training_Morphological_Analysis_Methods_and_Evaluating_them_Against_an_Expert_Consensus)

[4] Boudet, S., Houzé de l’Aulnoit, A., Demailly, R., Delgranche, A., Peyrodie, L., Beuscart, R., Houzé de l’Aulnoit,D. - A fetal heart rate morphological analysis toolbox for MATLAB. SoftwareX. 2020 Jan 1;11:100428. DOI:10.1016/j.softx.2020.100428 [Download on researchgate](https://www.researchgate.net/publication/339535549_A_fetal_heart_rate_morphological_analysis_toolbox_for_MATLAB)

[5] Boudet, S., Houzé de l’Aulnoit, A., Demailly, R., Peyrodie, L., Beuscart, R., Houzé de l’Aulnoit,D. - Fetal heart rate baseline computation with a weighted median filter. Computers in biology and medicine. 2019 Nov 1;114:103468. DOI:10.1016/j.compbiomed.2019.103468 [Download on researchgate](https://www.researchgate.net/publication/336035977_Fetal_heart_rate_baseline_computation_with_a_weighted_median_filter)

[6] Boudet, S., Houzé de l’Aulnoit, A., Demailly, R., Peyrodie, L., Houzé de l’Aulnoit,D. - Use of deep learning to detect the maternal heart rate and false signals on fetal heart rate recordings. Biosensors 2022; 12(9):691. DOI:10.3390/bios12090691 [Download on researchgate](https://www.researchgate.net/publication/363027076_Use_of_Deep_Learning_to_Detect_the_Maternal_Heart_Rate_and_False_Signals_on_Fetal_Heart_Rate_Recordings)

Patch note V2.1 - July 21th 2022
- Adding python source code for training FS methods
- Updating documentation

Patch note V2.0 - July 5th 2022
- Introduced the new sub-project of False Signal (FS) detection
- fhropen and fhrsave functions have change parameters to read other data available from CTG monitors (Sensor types and signal quality). Parameters are not in the same order than before.
- Train test and analyses data for baseline and A/D called FHRMA dataset, are now in FHRMAdataset subfolder
- interpolFHR now works also when Sig loss are coded NaN (or 0)
- Correct an incompatibility with other system than Windows. multisigfilter c function was compiled only for Windows and we replace with the Matlab built-in (but slower) function filtfilt in case this cases of other systems.
- Make a new start function FHRMA-menu and an APP to give an overview of this toolbox
