NEURON codes from the paper:

Non-invasive Suppression of Essential Tremor via Phase-Locked Disruption of its Temporal Coherence, by Schreglmann et al. (2020) Nature Communications.

This code was contributed by Xu Zhang. If you need platform specific or more help to run the model than provided here please consult this web page:

https://senselab.med.yale.edu/ModelDB/NEURON_DwnldGuide.cshtml

Usage:

For simulation of essential tremor without tACS, please refer to http://modeldb.yale.edu/257028.

A. CCTC model (Scaled-up CCTC model in the original Zhang & Santaniello 2019 model; see the link above) with phase-locked tACS: inside "CCTC_model/" folder

1) Compile all the .mod files in the folder "modfiles". For Windows users, use mknrndll.exe, then move the generated "nrnmech.dll" out of the modfiles directory (to where all the .hoc files are located). For Unix users, use the nrnivmodl command, then move the generated "x86_64" folder out of the modfiles directory (to where all the .hoc files are located).

2) Ensure the directories "recordings_full/" and "temporary/" exist (create them if not), so data can be stored after simulation ends.

3) The MATLAB script run_phase_locked_tACS.m runs one simulation of phase-locked stimulation with a particular pair of parameters (tACS amplitude & phase), and calls the operating system to run NEURON scripts (INIT_ET_tACS.hoc and ET_tACS.hoc). In the process, the MATLAB program iterates through each time step, estimates the tremor phase from the recordings of Vim TC neurons and adjusts the tACS phase accordingly in real time. Note that each simulation may take 2~3 hours on a desktop computer to finish.

4) The MATLAB script run_non_phase_locked_tACS.m sets up the parameter (tACS amplitude) and runs one simulation of non-phase-locked stimulation.

B. A further 5-fold scaled-up model with phase-locked tACS: inside "Scaleup_5x_all/" folder

The usage is exactly the same as above.

C. Example data for generating plots g(i-ii)~j(i-ii) of the corresponding supplementary figure. Execute plot_Example.m to reproduce the plots. Note that the MATLAB Wavelet toolbox is required for plotting the spectrograms in g(iV)~j(iv). Both plots can be previewed from Effective_phase.png and Ineffective_phase.png, respectively.