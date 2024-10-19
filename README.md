The accompanying code is used to analyze output from the ACCLIP convective transport diagnostic and produce figures for an upcoming publication in JGR-Atmospheres.

Code is currently preliminary but a final version of the code will be made available prior to final publication.

The process_airplane_traj3d.py script is designed to read the output from the convective transport diagnostic (available from https://doi.org/10.26023/DP1P-C32K-YJ05) and package it for use by the plotting scripts.

The OBSmatch.py script is provided because it contains useful subroutines.  

Figures are produced using the Python scripts which are named for the figure they produce.  The only exceptions are:
- Figures 2a and 2c are produced using the plot_airplane_traj3d.py script
- Figure 6 is produced using plot_tracertracer.py
