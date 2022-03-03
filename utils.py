import glob, os

import numpy as np
import astropy.io.fits as fits
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats

class Data:
	'''
	A simple class to hold data

	Using **kwargs to keep things future proof, see parameters in load_data 
	'''
	def __init__(self, filepath, **kwargs):
		# Load all relevant arguments that have been provided, set to defaults if not. 
		self.sci_ext = kwargs.get('sci_ext', 'SCI')
		self.err_ext = kwargs.get('err_ext', 'ERR')
		self.suffix = kwargs.get('suffix', 'calints.fits')
		self.x_ax = kwargs.get('x_ax', 2)
		self.y_ax = kwargs.get('y_ax', 1)
		self.t_ax = kwargs.get('t_ax', 0)

		# Load the data 
		self.load_data(filepath, self.sci_ext, self.err_ext, self.suffix, self.x_ax, self.y_ax, self.t_ax)

		# Load a seed image if possible 
		# if (self.instrument == 'NIRCAM') or (self.instrument == 'MIRI'):
		# 	self.load_seed()
		# else:
		# 	self.seed = None

		# Create a dictionary that we can save outputs to
		self.save_dict = {}

	def load_data(self, filepath, sci_ext='SCI', err_ext='ERR', suffix='calints.fits', x_ax=2, y_ax=1, t_ax=0):
		'''
		Function to load JWST pipeline Stage 2 equivalent data (e.g. calints.fits files)
		'''

		# Check if filepath is a directory, or a file
		if os.path.isdir(filepath):
			all_files = sorted(glob.glob(filepath+'*'+suffix))
		elif os.path.isfile(filepath):
			all_files = [filepath]
		else:
			raise ValueError('File path not recognised as a file or directory.')

		# Raise error if no files found.
		if not all_files:
			raise ValueError('No files found at provided filepath!')

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

					#Also, let's take this moment to get other information from the headers
					phead = hdul[0].header

					self.instrument = phead['INSTRUME'].upper()
				else:
					# Append data to the existing arrays (probably a smarter way to do this)
					all_sci = np.append(all_sci, sci, axis=t_ax)
					all_err = np.append(all_err, err, axis=t_ax)

		# If necessary, transpose array to a general format that the rest of the code will understand.
		if (x_ax != 2) or (y_ax != 1) or (t_ax != 0):
			all_sci = np.transpose(all_sci, axes=(t_ax, y_ax, x_ax))
			all_err = np.transpose(all_err, axes=(t_ax, y_ax, x_ax))

		# We can't use all the integrations (it would take too long), 
		# so pick out specific ones depending on the dataset. 
		if self.instrument == 'NIRCAM':
			self.use_ints = [1,20,32,47,-1]
			plot_aspect = 10
			vl, vh = 0.01, 0.6
		elif self.instrument == 'NIRISS':
			self.use_ints = [1,2,3,4,-1]
			plot_aspect=2
			vl, vh = 0.001, 0.05
		elif self.instrument == 'NIRSPEC':
			self.use_ints = [1,10,15,20,-1]
			plot_aspect = 5
			vl, vh = 0.01, 0.6
		elif self.instrument == 'MIRI':
			self.use_ints = [1,2,3,4,-1]
			plot_aspect=2
			vl, vh = 0.01, 0.6
		else:
			raise ValueError('Instrument:{} not recognised!'.format(self.instrument))

		# Save those specific integrations to the class
		self.sci = all_sci[self.use_ints,:,:]
		self.err = all_err[self.use_ints,:,:]

		# Explicitly free the memory being used to hold all of the data.
		# This might not be necessary, but who knows.
		del all_sci
		del all_err

		# Make a quick plot of the first integration image 
		plt.figure(figsize=(24,6))
		ax = plt.gca()

		ax.imshow(self.sci[0], aspect=plot_aspect, norm=LogNorm(vmin=vl*np.nanmax(self.sci[0]), vmax=vh*np.nanmax(self.sci[0])))
		ax.set_title('Wow, what a lovely {} spectrum!'.format(self.instrument), fontsize=16)
		ax.tick_params(axis='both', labelsize=16)

		return

	def basic_properties(self):
		# Function to calculate a range of basic properties for the images. 

		# Plot some profiles of the dispersion/cross-dispersion axes
		self.profiles(self.sci)

		# Get a histogram of all pixel values
		self.histogram(self.sci)

		# Get some simple quantitative measures from the images
		self.quantitatives(self.sci)

		return

	def seed_comparison(self):
		# Function to calculate comparisons to the seed images. 

		return

	def load_seed(self):
		# Function to load the seed image in as part of the class. 

		return


	def profiles(self, images):
		# Create 1D profiles from a given image
		fig, axs = plt.subplots(1, 2, figsize=(24,6), gridspec_kw={'width_ratios': [2, 1]})
		axs[0].set_title('Summed Dispersion Profile', fontsize=24)
		axs[1].set_title('Summed X-Dispersion Profile', fontsize=24)

		#Loop over the images we are interested in
		for image in images:
			disp = np.nansum(image, axis=0)
			xdisp = np.nansum(image, axis=1)

			axs[0].plot(range(len(disp)), disp)
			axs[1].plot(range(len(xdisp)), xdisp)

		plt.subplots_adjust(wspace=0.1)
		for ax in axs:
			ax.yaxis.get_offset_text().set_fontsize(18)
			ax.tick_params(axis='both', labelsize=18)

		axs[0].set_ylabel('Counts', fontsize=20)
		plt.show()

		return 

	def histogram(self, images):
		fig = plt.figure()
		ax = plt.gca()

		#Loop over the images we are interested in
		for image in images:
			image_1d = image.flatten()

			bins = np.logspace(np.log10(1e-3), np.log10(np.nanmax(image_1d)), 100)

			plt.hist(image_1d, bins=bins, histtype='step')

		# Use the last image to set the x limits
		ax.set_title('Histogram of Pixel Values', fontsize=16)
		ax.set_xscale('log')
		ax.tick_params(axis='both', labelsize=14)
		ax.set_xlabel('Pixel Flux / Counts', fontsize=14)
		ax.set_ylabel('Number of Pixels', fontsize=14)
		plt.show()

	def quantitatives(self):
		# Set up arrays to save things into 
		means = np.empty(len(self.use_ints))
		medians = np.empty(len(self.use_ints))
		vmaxs = np.empty(len(self.use_ints))
		vmins = np.empty(len(self.use_ints))
		stds = np.empty(len(self.use_ints))

		# Loop over the images we are interested in
		for i, image in enumerate(images):

			means[i] = np.nanmean(image)
			medians[i] = np.nanmedian(image)
			vmaxs[i] = np.nanmax(image)
			vmins[i] = np.nanmin(image)
			stds[i] = np.nanstd(image)

		# Print things out for people to see
		print('Integration #\'s: ', self.use_ints)
		print('Mean Values: ',means)
		print('Median Values: ',medians)
		print('Min Values: ',vmins)
		print('Max Values: ',vmaxs)
		print('Standard Dev Values: ',stds)

		# Now save to the global dictionary
		self.save_dict['means'] = means
		self.save_dict['medians'] = medians
		self.save_dict['vmins'] = mins
		self.save_dict['vmaxs'] = maxs
		self.save_dict['stds'] = stds

		return

	def residual_analysis():
		# TODO
		return

	def bad_pixels(self, sigma=10):

		# Get number of pixels in a single image
		npixels = self.sci[0].shape[0] * self.sci[0].shape[1]

		# Create an empty array to assign the number of outliers
		nints = self.sci.shape[0]
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
		for i, image in enumerate(self.sci):
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

		#Save to global dict
		self.save_dict['n_outliers'] = all_n_outliers

		return

	def background_trend():
		# TODO

		return

	def summary():
		# TODO

		return





