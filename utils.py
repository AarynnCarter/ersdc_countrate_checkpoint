import glob, os

import numpy as np
import astropy.io.fits as fits
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt

def load_data(filepath, sci_ext='SCI', err_ext='ERR', suffix='calints.fits', x_ax=2, y_ax=1, t_ax=0):
	'''
	Function to load JWST pipeline Stage 2 equivalent data (e.g. calints.fits files)
	'''

	# Check if filepath is a directory, or a 
	if os.path.isdir(filepath):
		all_files = sorted(glob.glob(filepath+'*'+suffix))
	elif os.path.isfile(filepath):
		all_files = [filepath]
	else:
		raise ValueError('File path not recognised as a file or directory.')

	# Loop over provided files
	for i, file in enumerate(all_files):
		# Open FITS file
		with fits.open(file) as hdul:
			# Read in science and error data using default / provided extensions. 
			sci = hdul[sci_ext].data
			err = hdul[err_ext].data

			if i == 0:
				# If this is the first file, let's initialise some arrays to start
				# saving data into 
				all_sci = sci
				all_err = err
			else:
				# Append data to the existing arrays (probably a smarter way to do this)
				all_sci = np.append(all_sci, sci, axis=t_ax)
				all_err = np.append(all_err, err, axis=t_ax)


	# Create a Data class that we can use for analysis.
	data = Data(all_sci, all_err, x_ax=2, y_ax=1, t_ax=0)

	return data

class Data:
	'''
	A simple class to hold data
	'''
	def __init__(self, sci, err, x_ax=2, y_ax=1, t_ax=0):
		self.sci = sci
		self.err = err
		self.x_ax = x_ax
		self.y_ax = y_ax
		self.t_ax = t_ax

	def residual_analysis():
		# TODO
		return

	def bad_pixels(self, nints=5, sigma=10):

		# Get number of pixels in image
		npixels = self.sci[0].shape[0] * self.sci[0].shape[1]

		# Create an empty array to assign the number of outliers
		if nints == 'all':
			nints = self.sci.shape[self.t_ax]

		all_n_outliers = np.empty(nints)

		# Function to calculate median and standard deviation of pixel values.  
		def identify_std(values):
			std = np.nanstd(values)
			return std
		def identify_median(values):
			median = np.nanmedian(values)
			return median

		# Describe a footprint to apply function (i.e. neighbouring pixels)
		footprint =  np.array([[1,1,1],[1,0,1],[1,1,1]])

		# Run filter functions
		for i, image in enumerate(self.sci[0:nints]):
			# Get standard deviation and median images
			stds = ndimage.generic_filter(image, identify_std, footprint=footprint, mode='constant', cval=np.nan)
			meds = ndimage.generic_filter(image, identify_median, footprint=footprint, mode='constant', cval=np.nan)

			# Identify outliers
			outlier_image = np.divide(np.subtract(image,meds), stds)

			# Count outliers
			n_outliers = np.count_nonzero(outlier_image > 10)

			# Save to array
			all_n_outliers[i] = n_outliers

		plt.scatter(range(1, len(all_n_outliers)+1), all_n_outliers/npixels*100)
		plt.xlabel('Integration #', fontsize=16)
		plt.ylabel('Bad Pixel % (10$\sigma$)', fontsize=16)

		return

	def background_trend():
		# TODO

		return

	def summary():
		# TODO

		return





