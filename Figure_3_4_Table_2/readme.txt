Process raw NOAA OI data to regional subsets that contain SSTA and MHW categories:

	parallel_process_data_to_regional_blocks.sh
	Contacternates files into continuous time series in blocks of 30o x 20o (lon x lat)
	Automatically spans multuple parallel jobs, 1 for each blocl


	parallel_process_MHW_regional_blocks.sh
	Parallel processing of regional blocks to calculate MHW statistics using:
	Uses: marineHeatWaves90.py

	parallel_process_MHW_regional_blocks_98pc.sh
	As above but uses the 98% criteria used for Table 2
	Uses: marineHeatWaves98.py


Generation of dat for Table 2:
	MHWControlGUI_updated.m
	MHWControlGUI_updated.fig
	Above code uses the processed region MHW data and opens it in a GUI that calculates the size, maximum intensity etc of manually selected MHW. Data derived from the output of this function are used to populate table 2


Figure 3:
summary_fig_detrended_corrected.m
Collates information for significant increases/decreases in MHW day occurrence associated with different climate indices and generates fig. 4
summaryFig_modesVSmhw_detrended_corrected_aug2018.m 
Generates data for the above script. Calculates increase/decrease in marine heatwave days at each grid point and tests if the change is significant (based on Monte Carlo test)
load_modes.m
Used in the above. Loads climate indices (based on the following data files)
NAO.txt
AMO.txt
nino34.txt
PDO.txt
TPI_IPO.txt
ANino.txt
SAM.txt
MODOKI.txt
DMI.txt
NPGO.txt

parallel_process_MHW_regional_blocks_2degree.sh
Calculate regional MHW statistics from raw NOA OI SSST data for use in summary_fig_detrended_corrected.m
Uses: regional_MHW_pc90_reducedFileSize_2degree.py


Figure 4:
plot_regional_anomaly_drivers_corrected.m
Processes data for Figure 4: MHW days associated with each mode/region
regional_anomaly_drivers_cummulatice_stats_detrended_corrected.m
Generates data for above script.Calculates increase/decrease in marine heatwave days for each region and tests if the change is significant (based on Monte Carlo test)

