Hello Dr. Wakin!

Welcome to the code from our project. 
There are many files in here, but to reproduce our results you need only run test_signals. 
This file tests how each of our algorithms reconstructs a signal, giving the number of terms and the error.

Knobs to turn in test_signals:
t: 
	specifically the middle value defining the timestep (Fixed_Coefficients Algorithm takes a LONG time, comment out the sections that run/plot it if needed)

signals: 
	comment or uncomment different signals we tested with, or try making your own

metric_def_list: 
	performance for all algorithms is highly dependent on what basis functions are used. Maybe you can find some that work better than what we chose.

percent_error_threshold: 
	change the error threshold each of the algorithms aims for in their reconstruction. Default is 10%. Fixed_Coefficients has issues with low error thresholds and will take care of itself if this occurs, it simply won't meet the error specified in some cases.

Note that changing the code might yield unexpected results, this hasn't been tested thoroughly for useability, apologies.