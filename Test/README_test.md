The test case for the N2 gas temperature fitting script 'Kris_Gas_Fit.m' consists of three data files and three corresponding reference data files. The goal is to run a test case of the script that should produce the two included images ('N2_Fit_test_iterations.jpeg' and 'N2_Fit_Test_Result.jpeg'). 
The data was taken using a 13.56 MHz atmospheric pressure plasma jet operating at 8W plasma power using helium gas. The source of N2 is the ambient atmosphere. The spectrometer used is a Princeton Instruments ACTON2500 with a PIXIS 256 ICCD camera attached.

The general operation of the script is fairly straightforward. These scripts are written for Matlab and therefore a copy of Matlab is a prerequisite. Be sure that all data files and Matlab scripts are in the same folder. 
First, run 'testdata_test.m' which will import the test datasets and compile them into the 'testdata.mat' file. This file will act as the input for 'Kris_Gas_Fit_test.m'.
Then run 'Kris_Gas_Fit_test.m', it should take upwards of a minute to complete. While running some text will be displayed on the console reading:
working1working2working3working4working5working6working7working8working9working10
This is an internal update letting the user know that the temperature range selected is being converged on. Text 'failingup' or 'failingdown' will appear if the range of gas temperatures being simulated is changing in search of a best fit. This should let the user know if the correct range has been selected. 
In the test the correct range has been selected, after running the test case the user can change the variable 'T_rotation' on line 37 to see how the temperature range will automatically shift to find the best fit. 

Some variables that can be changed if needed: 

  num_bootstrap :: changing this will change the number of iterations to be simulated according to the monte-carlo bootstrap bethod. For details on this method please see Kristopher Fords thesis.
  
  Jmax :: changing Jmax will alter the number of wavenumbers that are simulated for the rotational bands. It is not suggested to change this.
  
  T_rotation :: the temperature range to begin the simulation, including the granularity of the simulation. 
  
  fwhm :: the FWHM of the spectrometer used to collect the data in nm.

LIMITS: 
If the resulting rotational temperature is >1000K or <200K the error 'temperature_limit_exceeded' will appear. This is arbitrary and can be user selected. Line: 273
When the temperature range changes, the default step change for the range is set to 150K. This can be changes on lines (261,265).

PLOTTING: 
The bottom portion of the script has a number of pre-built plots that the user can select from. Feel free to change what is needed. 
