This folder contains the Fortran code needed to implement an ABC-SMC algorithm
to estimate transmission parameters for a high mortality disease using
mortality data.

The methods are methods described in Guinat et al. (2018) Transboundary and 
Emerging Diseases 65(2), e264-e271 (doi: 10.1111/tbed.12748).

The individual files and how to run the compiled code is described in
"ABC-SMC Code Description.pdf"

The data for African swine fever are in the "ASFData" branch:
DayOfConfirmation_Russia.txt - day of mortality data when ASFV confirmed in herd
InitialHerdSizes_Russia.txt - initial size for each herd
MortalityData_Russia_Herd[n].txt - mortality data for herd [n]
SMCInitFile.txt - example of file to initialise SMC algorithm
