Read-me

This repository serves as an addendum to the following peer-reviewed article pulished in the AAPS journal (2021):

---------------------------------------------------------------------------------------------------------------------------
A novel approach for quantifying the pharmacological activity of T-cell engagers utilizing in vitro time-course experiments 
and streamlined data analysis

Arthur J Van De Vyver, Miro Eigenmann, Meric Ovacik, Christian Pohl, Sylvia Herter, Tina Weinzierl, Tanja Fauti, 
Christian Klein, Thorsten Lehr, Marina Bacac, Antje-Christine Walz

---------------------------------------------------------------------------------------------------------------------------

This repository contains two Python files that serve as tools to explore the pharmacology of CD3-bispecific antibodies

1) time-independent PKPD analysis - AAPS.py

	This file allows the reader to upload their own experimental dataset and perform a time-independent analysis.
	The dataset needs to be prepared in Excel (.xlsx) format and contain the following information in separate columns:
		a) measured time points 	(column name 'Time')
		b) tested drug concentrations 	(column name 'Conc')
		c) tested compounds 		(column name 'Drug')
		d) tested cell lines 		(column name 'CellLine')
		e) readouts of interest 	(column name 'Readout')
		f) observed data 		(column name 'Obs') 
		
	The path to the data file can be specified on line 231
	It is encouraged to explicitly name the experimental conditions and compounds in the dataset as to allow easier 
	navigation through and in the output files 
	Please provide an empty folder ('Figures') to store the output figures 

	The outputs of the time-independent analysis include:
		a) plots of raw data per experimental conditions 
			(denoted by the readout and experimental condition)
		b) EC50 estimations at different time points 
			(denoted ‘Potency Change Over Time’) 
		c) EC50 estimation across different experimental conditions including different cells lines and different drugs 
			(denoted ‘Cumulative Potency’) 
		d) sigmoidal model fits for each experimental condition 
			(denoted ‘Model Fit’)
		e) estimated sigmoidal model parameters for the model fits including RSE% values 
			('Sigmoidal Model Parameters.txt')
		f) simulations for sigmoidal fits for each experimental condition 
			('Sigmoidal Model Simulations.xls') 
		g) calculated AUCE values (AUCE.xls) as well as Hockey-Stick model fits 
			(denoted with 'HOCKEYSTICK') 
		h) estimated model parameters including RSE% values
			('Hockey Stick Parameters.txt') 
		i) simulations for Hockey stick fits for each experimental condition 
			('Hockey Stick Model Simulations.xls')

2) QE_synapse_calculation.py

	Based on the equations derived by Schropp and colleagues (1) to calculate the number of immune synapses 
	(i.e., trimeric complexes) between CD3-bispecific antibodies, 
	tumor cells, and T-cells under quasi-equilibrium (QE) assumptions. 
	The equations are summarised in 'QE equations summary.pdf' 

	The code enables the calculation of the concentration of immune synapses formed under the experimental conditions 
	specified by the user.

	If desired, a plot can be generated showing immune synapse formation over a range of concentrations and 
	visualise the bell-shape. Please provide an empty folder ('Figures') to store the output figure.
	An additional plot can be generated to compare the simulations under QE assumptions with those of the full ODE model
	
	      

References:

(1) 	Schropp J, Khot A, Shah DK, Koch G. Target-Mediated Drug Disposition Model for Bispecific Antibodies: Properties, 
	Approximation, and Optimal Dosing Strategy. CPT: Pharmacometrics & Systems Pharmacology. 2019;8(3):177-87. 
	doi: 10.1002/psp4.12369.
