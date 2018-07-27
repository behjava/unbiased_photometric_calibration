# unbiased_photometric_calibration
Unbiased photometric calibration

Uses maximum likelihood method to calibrate photometry against a set of standard magnitudes. The sample code is for luminance filter.
The calibration starts with a first estimation of a,b, and c and rejects outliers from a given standard deviation (2 or 3). The fitting and the outlier rejection continues until there is no outlier. The final fitting is performed with higher precision (i.e. finer parameter grids).
Python modules: numpy, scipy

# It was used for calibration of data for this research: http://adsabs.harvard.edu/abs/2016A%26A...588A..89J

You are required to acknowledge the use of this code and also cite the above publication if you use the contents of this repository. 
Please contact me in case you had questions.

