# TTA-UC_simulation
MATLAB scripts for computational modelling in support of the paper titled: "Identifying efficiency-loss pathways in triplet-triplet annihilation upconversion systems".
The two MATLAB scripts that are provided in this repository are designed to simulate the steady-state and time-dependent emission characteristics of a TTA-UC system

Steady-state analysis

1. Open 'FullFssShow.m' in MATLAB

2. Define all relevant variables in the first section. Input the value of the variable directly in the appropriate line. For example a value for the annihilator intrinsic triplet quenching rate constant of '1000 /s' can be written as '1e3' in line 12. If you want to disregard a particular process, please set the relevant rate constant to 0. For example, if you do not want to consider FRET, then set kfret to 0.
 
3. Run lines 24 to 57 to obtain steady-state upconverted fluorescence rates, steady-state sensitizer phosphorescence rates, emission ratios, and quantum yields as a function of irradiance. The irradiance range you want to simulate is specified in line 28. You can increase or decrease the range, as well as change the number of datapoints here. The default is 0.001 mW/cm2 to 1000000 mW/cm2. 

4. You can then plot any of the values I mention in step 3 as a function of irradiance. 

5. For simulating parasitic sensitizer reabsorption effects, run lines 59 to 92 in the next section.

Time-dep analysis

1. Define variables just as you did previously for the steady-state analysis. Note that the timestep dt is a new variable. My suggestion is to use a value of 1 ns (1e-9) as a start. 

2. Run section 2. Here, you can change the value of irradiance in the first 2 for-loop blocks as needed. Additionally for-loop blocks can also be added to simulate the continued decay in UCPL or sensitizer phosphorescence. The length of the for-loop can also be changed to increase the duration of the simulation. Note that the run-time of this section can vary based on the complexity of the final model (e.g.: whether or not all UEL mechanisms are active). It may take several minutes for the MATLAB to complete running this section.
 
3. Run section 3 to plot your data and also calculate effective annihilator triplet lifetimes. 

**Additional notes: Choosing a smaller timestep dt improves the accuracy of the simulation and allows you to resolve faster kinetic processes such as prompt fluorescence or rapid TET. To compensate for the short duration of the simulation with small dt, the length the for-loop blocks may be increased, however doing so will also greatly increase time it takes for MATLAB to complete running the simulation. 
