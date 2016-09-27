"""K-Index Calculator
   ==================

    A module for calculation of the K-Index from horizontal geomagnetic components.

	This module uses the FMI method, details of which can be found here:
	http://swans.meteo.be/sites/default/files/documentation/TN-RMI-2010-01_K-LOGIC.pdf

	Functions Included
	-----------------------------------------------------------------

	Main Functions
	-------------------

		KIndexSuperCalc
		- Calculates the K-index using the FMI method.

		KIndexPlotter
		- Plots K-Index values in a nice Bar Plot

		BxByBzPlotter
		- Plots Bx, By and Bz nicely

		OtherPlotter
		- Plots Declination, Horizontal and derivative of Horizontal nicely

	Secondary Functions
	-------------------

		Time2Float
		- Converts datetime object or array of datetime objects to floats.

		Float2Time
		- Converts float or array of floats to datetime objects.

		MinuteBin
		- Bin second data into minutes.

		KIndex
		- Calculate K-Index according to the FMI Method.

		FMISmooth
		- Calculates the Solar Regular (Sr) curve according to the FMI Method.

		SrSmooth
		- Smooths the solar regular curves.

		Subtracted
		- Subtracts the smoothed solar regular curves from the minute binned data.

		KIndexBarColor
		-Colours K-index bar plot so it looks nice.

		H
		- Gets horizontal component from bx and by

		Declination
		- gets declination in radians and degrees

		Slope
		- gets derivative of magnetic components


	Usage
	-------------------
	Assuming you have 4 days of geomagnetic data in the following format:
	
	Times (array of datetime objects in seconds)
	Bx, By, Bz (array of geomagnetic data in seconds)

	First, convert the Times array into second floats:

	>>> Times_float = Time2Float(Times)

	Now convert all of the data to minute bins:

	>>> minute_time, minute_bx, minute_by, minute_bz = MinuteBin(Times_float, Bx, By, Bz)

	To get the K-Index for these 4 days:

	>>> k_index, k_time, order = KIndexSuperCalc(minute_time, minute_bx, minute_by, n)

	where n is the maximum threshold for a K9 event (dependent on latitude).

	NOTE: The first few values in the array k_index are likely to be inaccurate.
	This is because of the way the K-Index is calculated. It is safest to dismiss the first 8
	calculated K-Index values.

	To plot the K-Index:

	>>> KIndexPlotter(k_index, k_time, m)

	where m is the figure number.


	Author
	-------------------

	Written by Sean Blake in Trinity College Dublin, 2014-2016.

	Email: blakese@tcd.ie

	GITHUB: https://github.com/TerminusEst

	Uses the MIT license.

"""

import datetime
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
##########################################################################
##########################################################################
##########################################################################

def H(minute_bx, minute_by):
	"""Returns the horizontal magnetic component

		Parameters
		-----------
		minute_bx, minute_by = arrays of NS and EW magnetic components

		Returns
		-----------
		minute_bh = array of horizontal magnetic component
		
		-----------------------------------------------------------------
	"""
	return np.sqrt(array(minute_bx)**2 + array(minute_by)**2)

def Slope(minute_time, *args):
	"""Returns dBx/dt, dBy/dt of geomagnetic components (in nT/min)

		Parameters
		-----------
		minute_time = array of FMI cleaned time
		*args = array(s) of magnetic components

		Returns
		-----------
		output = list of derivative magnetic components

		-----------------------------------------------------------------
	"""
	output = []
	new_time = arange(minute_time[0], minute_time[-1], 60)

	for value_array in args:
		a = numpy.interp(new_time, minute_time, value_array)
		sloped = np.gradient(a)
		sloped = np.concatenate((sloped, array([sloped[-1]])))

		output.append(sloped)

	return output

def Declination(minute_bx, minute_by):
	"""Returns geomagnetic declination in radians and degrees

		Parameters
		-----------
		minute_bx, minute_by = arrays of horizontal magnetic components

		Returns
		-----------
		D_rad = array of declination in radians
		D_deg = array of declination in degrees

		-----------------------------------------------------------------
	"""
	minute_bh = H(minute_bx, minute_by)
	D_rad = [math.asin(y/h) for y, h in zip(minute_by, minute_bh)]
	D_deg = [(180*x/math.pi) for x in D_rad]

	return D_rad, D_deg

def OtherPlotter(minute_time, minute_bx, minute_by, n, titlez = "", stamp = ""):
	"""Plots Declination, Horizontal and Derivative of Horizontal in a nice plot

		NOTE: Plots 3 days of data.

		Parameters
		-----------
		minute_time = array of float time values
		minute_bx, minute_by = arrays of magnetic components
		n = number of figure
		titlez = title of figure
		stamp = string of stamp on bottom left of plot (e.g., "TCD/DIAS")

		Returns
		-----------
		fig = Plot of data
		-----------------------------------------------------------------
	"""

	# First get declination, bh and dbh:
	D_rad, D_deg = Declination(minute_bx, minute_by)
	minute_bh = H(minute_bx, minute_by)
	dbh = Slope(minute_time, minute_bh)[0]


	# Start with the figure
	fig = plt.figure(n)
	plt.clf()
	plt.subplots_adjust(bottom = 0.1, top = 0.93, left = 0.1, right = 0.94, hspace=0.1)

	# This is for the timestamp
	nowz = datetime.datetime.now()
	timestamp_string = "Plot updated {}/{}/{} {} UT".format(nowz.day, nowz.month, nowz.year, str(nowz)[10:19]) 

	# This is the first day of the data (shows last three days)
	start_day = minute_time[0] - minute_time[0]%(24*60*60)

	start = Float2Time(start_day) + datetime.timedelta(1)
	middle1 = start + datetime.timedelta(1)
	middle2 = start + datetime.timedelta(1)
	end = start + datetime.timedelta(3)

	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
	# First subplot: Declination

	a1 = plt.subplot(311)
	plt.ticklabel_format(useOffset=False)
	plt.plot(Float2Time(minute_time), D_deg, 'r')

	plt.locator_params(axis='y', nbins=4)
	plt.ylabel("D (degrees)", fontsize = 14)
	plt.title(titlez, fontsize = 16)

	plt.setp(a1.get_xticklabels(), visible=False)
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:00'))
	plt.gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,6)))
	plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
	plt.tick_params('x', length=5, width=1, which='both')

	plt.grid(True)

	for i in range(4):
		plt.axvline(Float2Time(start_day + (i*24*60*60)), lw = 2, color = "black")

	xlim([start, end])

	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
	# Second subplot: Horizontal component

	a2 = plt.subplot(312)
	plt.ticklabel_format(useOffset=False)
	plt.plot(Float2Time(minute_time), minute_bh, 'b')

	plt.locator_params(axis='y', nbins=5)
	plt.ylabel("H (nT)", fontsize = 14)

	plt.setp(a2.get_xticklabels(), visible=False)
	H_up_lim = (int(max(minute_bh))/10)*10+10
	H_low_lim = (int(min(minute_bh))/10)*10-10
	plt.ylim([H_low_lim, H_up_lim])

	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:00'))
	plt.gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,6)))
	plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
	plt.tick_params('x', length=5, width=1, which='both')

	for i in range(4):
		plt.axvline(Float2Time(start_day + (i*24*60*60)), lw = 2, color = "black")

	xlim([start, end])

	plt.grid(True)

	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
	# Third subplot: Derivative of Horizontal component

	a3 = plt.subplot(313)
	plt.ticklabel_format(useOffset=False)
	plt.plot(Float2Time(minute_time), dbh, 'g')

	plt.locator_params(axis='y', nbins=4)
	plt.xlabel('{}-{}-{}                         {}-{}-{}                         {}-{}-{}'.format(start.day, start.month, start.year, middle1.day, middle1.month, middle1.year, middle2.day, middle2.month, middle2.year), fontsize = 14)
	plt.ylabel("dH/dt (nT/min)", fontsize = 14)

	plt.xlim([start, end])
	plt.ylim([-5, 5])

	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:00'))
	plt.gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,6)))
	plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
	plt.tick_params('x', length=5, width=1, which='both')

	for i in range(4):
		plt.axvline(Float2Time(start_day + (i*24*60*60)), lw = 2, color = "black")

	plt.grid(True)

	plt.figtext(0.65, 0.02, timestamp_string, fontsize = 11, style='italic')
	plt.figtext(0.02, 0.02, stamp, fontsize = 11, style = 'italic')

	return fig

def BxByBzPlotter(minute_time, minute_bx, minute_by, minute_bz, n, titlez = "", stamp = ""):
	"""Plots Bx, By and Bz in a nice plot

		NOTE: Plots 3 days of data.

		Parameters
		-----------
		minute_time = array of float time values
		minute_bx, minute_by, minute_bz = arrays of magnetic components
		n = number of figure
		titlez = title of figure
		stamp = string of stamp on bottom left of plot (e.g., "TCD/DIAS")

		Returns
		-----------
		fig = Plot of data
		-----------------------------------------------------------------
	"""

	# This is for the timestamp
	nowz = datetime.datetime.now()
	timestamp_string = "Plot updated {}/{}/{} {} UT".format(nowz.day, nowz.month, nowz.year, str(nowz)[10:19]) 

	# This is the first day of the data (shows last three days)
	start_day = minute_time[0] - minute_time[0]%(24*60*60)

	start = Float2Time(start_day) + datetime.timedelta(1)
	middle1 = start + datetime.timedelta(1)
	middle2 = start + datetime.timedelta(1)
	end = start + datetime.timedelta(3)

	# Start with the figure
	fig = plt.figure(n)
	plt.clf()
	plt.subplots_adjust(bottom = 0.1, top = 0.93, left = 0.1, right = 0.94, hspace=0.1)

	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
	# First subplot: Bx

	a1 = plt.subplot(311)
	plt.ticklabel_format(useOffset=False)
	plt.plot(Float2Time(minute_time), minute_bx, 'r')

	plt.locator_params(axis='y', nbins=6)
	plt.ylabel("Bx (nT)", fontsize = 14)
	plt.title(titlez, fontsize = 16)

	plt.setp(a1.get_xticklabels(), visible=False)
	plt.xlim([start, end])
	Bx_up_lim = (int(max(minute_bx))/10)*10+10
	Bx_low_lim = (int(min(minute_bx))/10)*10-10
	plt.ylim([Bx_low_lim, Bx_up_lim])

	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:00'))
	plt.gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,6)))
	plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
	plt.tick_params('x', length=5, width=1, which='both')
	for i in range(4):
		plt.axvline(Float2Time(start_day + (i*24*60*60)), lw = 2, color = "black")

	plt.grid(True)

	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
	# Second subplot: By

	a2 = plt.subplot(312)
	plt.ticklabel_format(useOffset=False)
	plt.plot(Float2Time(minute_time), minute_by, 'b')

	plt.locator_params(axis='y', nbins=5)
	plt.ylabel("By (nT)", fontsize = 14, labelpad = 0)

	plt.setp(a2.get_xticklabels(), visible=False)
	plt.xlim([start, end])
	By_up_lim = (int(max(minute_by))/10)*10+10
	By_low_lim = (int(min(minute_by))/10)*10-10
	plt.ylim([By_low_lim, By_up_lim])

	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:00'))
	plt.gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,6)))
	plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
	plt.tick_params('x', length=5, width=1, which='both')
	for i in range(4):
		plt.axvline(Float2Time(start_day + (i*24*60*60)), lw = 2, color = "black")

	plt.grid(True)

	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
	# Third subplot: Bz

	a3 = plt.subplot(313)
	plt.ticklabel_format(useOffset=False)
	plt.plot(Float2Time(minute_time), minute_bz, 'g')

	plt.locator_params(axis='y', nbins=5)
	plt.xlabel('{}-{}-{}                         {}-{}-{}                         {}-{}-{}'.format(start.day, start.strftime("%B")[0:3], start.year, middle1.day, middle1.strftime("%B")[0:3], middle1.year, middle2.day, middle2.strftime("%B")[0:3], middle2.year), fontsize = 14)
	plt.ylabel("Bz (nT)", fontsize = 14)

	plt.xlim([start, end])
	Bz_up_lim = (int(max(minute_bz))/10)*10+10
	Bz_low_lim = (int(min(minute_bz))/10)*10-10
	plt.ylim([Bz_low_lim, Bz_up_lim])

	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:00'))
	plt.gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,6)))
	plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
	plt.tick_params('x', length=5, width=1, which='both')
	for i in range(4):
		plt.axvline(Float2Time(start_day + (i*24*60*60)), lw = 2, color = "black")

	plt.grid(True)

	plt.figtext(0.65, 0.02, timestamp_string, fontsize = 11, style='italic')
	plt.figtext(0.02, 0.02, stamp, fontsize = 11, style = 'italic')

	return fig

