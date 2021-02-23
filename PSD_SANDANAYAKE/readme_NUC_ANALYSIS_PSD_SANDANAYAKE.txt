NUC_ANALYSIS.C by P.S.D. Sandanayake
************************************

This program is a ROOT macro which takes user inputs on how to handle the level scheme data files
(ex. TE128POS) and spectrum data files (ex. te128pos_m1) and calculates the level density plots,
B(XL) plots and combine them to produce f_xl plots for Te128 and Sc44 nuclei.

For the program to function correctly, 
	- The data files must be of no extension or .txt
	- The data files must be located in the same directory as the program
To run it,
	- enter the command "root -l" into the terminal (all commands must be entered without quotes (")).
	- when root has started, type ".x NUC_ANALYSIS.C"

******* Start of phase I: LEVEL DENSITY plots *******

	- Now the user will see a question "Input the STARTING STRING of the files you need to analyze
	  (IN UPPERCASE. For this project, please enter either TE128 or SC44):". The input MUST BE
	  in uppercase in this project because both files either start with "TE128" or "SC44" 
	  concerning level densities.
	- Then it requires the user to provide the bin width (which is a sensible value of course). In
	  this project, it is kept at 0.2MeV as a standard.
	- Then it asks the user "Do you wish to make plots specific in J for level density (1 if YES, 
	  0 if NO)?". Thus input 1 or 0 must be input as per user choice.
	- Once this is done, a ROOT canvas will pop up with the total level density plot by default.
	  Keep hitting ENTER on the active canvas until it reveals the Jmax values for either parities
	  for the selected nucleus. This serves as a reference of limits of J for which the user will
	  make the next choice.
	- If the user choice was "1", he will be asked how many J_i values does he wish to analyze,
	  and for the given no. of choices, the program will ask J_i and parity one after the other.
	  For parity, he must enter "1" for + parity and "-1" for - parity.
	- Following the previous print, as he keeps hitting ENTER, the level density plots will be
	  revealed. All relevant data files for all histograms will be automatically saved in the 
	  directory along with the plots in PDF format.

******* Start of phase II: B(M1) plots *******

	- Now, depending on the nucleus, user will be asked differently here. For Te128, he will be
	  asked to choose a specific parity to continue, whereas for Sc44 it will not be the case.
	- For Te128, user has to input "1" for + parity, "-1" for - parity.
	- Next question will be regarding making cuts in excitation energy. user must input "1" if YES
	  or "0" if NO. If YES, he will be asked for the number of cuts he needs to make, and will be
	  cycled through to input the upper and lower limits of energy for the cuts in order.
	- Then he will be asked if he wishes to make selections on J_i for which again he must answer
	  "1" if YES, "0" if NO. If YES, he will be asked for the no. of J_i and also the list of J_i.
	- At this point a plot will be displayed on the canvas again. Keep pressing ENTER to reveal 
	  all the plots as per the user's choices. Again, all relevant data files will be stored in 		  
	  the directory along with a PDF file with the plot.

******* Start of phase III: nuclear strength plots *******

	- When you enter here, the canvas will directly display the de-excitation function plots. For
	  Te128, de-excitation strength functions for + parity, - parity and combined will be revealed
	  as you hit ENTER on the active canvas. For Sc44 only the total de-excitation strength 
	  function will be plotted because it is an E1 case. Again, all relevant data wil be printed
	  to files and plots to PDF files accordingly in the directory.

*****IMPORTANT*****

	- In case you observe weird plots (or any unexpected errors), this is due to memory leaks.
	  Therefore it is expected that the program is executed a second time only after exiting
	  ROOT environment and repeating the above procedure.
	  
