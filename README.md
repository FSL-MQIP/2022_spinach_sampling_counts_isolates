# 2022_spinach_sampling_counts_isolates

Code for the paper titled: A longitudinal study on the bacterial quality of baby spinach cultivated in Arizona and California

1.	The data and code for this study are divided into a “counts” folder (containing bacterial enumeration data and code) and an “isolates” folder (containing bacterial isolates’ data and code).
2.	All data, in the “counts” and “isolates” folders, are provided in the “raw” form (i.e., collected data) and the “wrangled” form (i.e., after structuring the data for analysis).
3.	The R files are numbered in sequence of analysis. For example, in the “isolates” folder, the first file “01_parsing_data.R” includes code for wrangling the raw data; the wrangled data is then used for analysis in the subsequent R files (i.e., “02_statistical_test.R” and “03_plotting_data.R”). Thus, the R files need to be run in order, starting from the file numbered 1, when repeating these analyses.
4.	The raw weather data and corresponding code (“02_weather_wrangling.R”) are not uploaded to Github as they contain information that can identify industry collaborators
5.	Code for key analyses in the paper can be found in:
5a.	“counts/03_regression.R”, which contains code for the H, APC and PC models.
5b.	“counts/06_analysis_of_growth_parameters.R”, which contains code for analyzing the growth parameters fit to the packaged spinach.
5c.	“isolates/02_statistical_test.R”, which contains code for all analyses conducted on the isolates’ data.
5d.	Other analyses can be found in the remaining R files in this repository.
