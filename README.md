# dosecombinations

Dose Finding Studies for Combination Therapies: Guide to R Code Files


Designs:

KEY (Keyboard)
BOIN
PIPE
SF (Surface-Free)
BLRM

Files:

For each design, the prior calibrations are executed in the file ‘<METHOD>_priorcal.R’ and the simulations are executed in the file ‘<METHOD>_sims.R’.

(For the BLRM there are the additional ‘BLRM_priorcal_ewoc.R’ and ‘BLRM_priorcal_out.R’ files. Due to the computational intensity of the prior calibrations, separate files are required to calibrate the prior distribution parameters, overdose control and to determine the best values.)

Graphs are drawn in the file ‘graphs.R’

The case study is explored in the files ‘case_study_gandhi_stop.R’ and ‘case_study_gandhi_no_stop.R’, corresponding to using a stopping rule and not respectively. 'patient_genration_case_study.R' contains functions fo patient generation for the case study. All methods are implemented in each file.
