# VFN-Lab #

Code used for running the VFN bench and analyzing the results.

## \AnalysisCode\ ##
    * Contains a library of functions (AnalysisLib) for 
    reading the VFN scan data and analyzing it. 
    * Also contains sample scripts showing how an analysis 
    is done.
    * One of the sample scripts is used for doing a LSQ Fit.
    * Contains another folder (ZernikeAnalysis) for analyzing
    wavefronts and decomposing them into a Zernike basis. 
    This is useful for processing Zygo or other WF data.

## \ControlCode\ ##
    * Contains control code for the various devices in the 
    the VFN bench.
    * Also has control scripts used to do the VFN scans.
    * Also has control scripts for miscellaneous tests.
    
## \IDSandThorlabsCameraCode\ ##
    * Code for using the IDS or Thorlabs camera. This includes
    the simple script (command-line) and GUI.
