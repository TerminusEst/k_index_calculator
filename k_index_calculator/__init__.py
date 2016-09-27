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
import math
import datetime
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator

################################################################################
################################################################################

def Time2Float(x):
	"""Converts datetime object or array of datetime objects to floats.

		Parameters
		-----------
		x = single or array of datetime objects

		Returns
		-----------
		y = single or array of floats

		-----------------------------------------------------------------
	"""
	if (type(x) == np.ndarray) or (type(x) == list):
		emptyarray = []
		for i in x:
			z = (i - datetime.datetime(1970, 1, 1, 0)).total_seconds()
			emptyarray.append(z)
		emptyarray = np.array([emptyarray])
		return emptyarray[0]
	else:
		return (x - datetime.datetime(1970, 1, 1, 0)).total_seconds()

#-------------------------------------------------------------------------------

def Float2Time(x):
	"""Converts float or array of floats to datetime objects.

		Parameters
		-----------
		x = single or array of floats

		Returns
		-----------
		y = single or array of datetime objects

		-----------------------------------------------------------------
	"""
	if (type(x) == np.ndarray) or (type(x) == list):
		emptyarray = []
		for i in x:
			z = datetime.datetime.utcfromtimestamp(i)
			emptyarray.append(z)
		emptyarray = np.array([emptyarray])
		return emptyarray[0]
	else:
		return datetime.datetime.utcfromtimestamp(x)

#-------------------------------------------------------------------------------

def MinuteBin(timedate_float, bx, by, bz):
	"""Bin second data into minutes

		Parameters
		-----------
		timedate_float = array of second data in float
		bx, by, bz = arrays of magnetic data

		Returns
		-----------
		minute_time = array of minute data in float
		minute_bx, minute_by, minute_bz = binned magnetic data

		-----------------------------------------------------------------
	"""
	# n = number of days in data
	n = int((timedate_float[-1] - timedate_float[0])/(24*60*60)) + 1

	# Gets the start of the day in seconds
	day_seconds = int(timedate_float[0])-int(timedate_float[0])%(24*3600)

	# Creates array of minutes
	minutes = np.arange(0, n * 1440)
	minutes = (minutes * 60) + day_seconds

	# master is numpy array with columns for bx, by, bz, count and times
	master = np.zeros((n*1440, 5))
	master[:,-1] = minutes

	# loop over times
	for i, v in enumerate(timedate_float):
		# check which master row it belongs to
		index = int((v - day_seconds)/60) #- 1
		# add to each column
		try:
			master[index][3] += 1
			master[index][0] += bx[i]
			master[index][1] += by[i]
			master[index][2] += bz[i]
		except:
			continue

	# get non-zero indices
	indi = master[:,3] > 0

	minute_bx = master[:,0][indi] / master[:,3][indi]
	minute_by = master[:,1][indi] / master[:,3][indi]
	minute_bz = master[:,2][indi] / master[:,3][indi]
	minute_time = master[:,4][indi]
	
	return minute_time, minute_bx, minute_by, minute_bz

#-------------------------------------------------------------------------------

def KIndex(minute_time, minute_bx, minute_by, k9):
	"""Calculate K-Index according to the FMI Method
		
		Detailed instructions can be found here:
		http://swans.meteo.be/sites/default/files/documentation/TN-RMI-2010-01_K-LOGIC.pdf

		Parameters
		-----------
		minute_time = array of minute data in float
		minute_bx, minute_by = arrays of horizontal magnetic data
		k9 = lower limit for K9 event

		Returns
		-----------
		k_index = array of k_index values (0's are represented as 0.25)
		timestamp = array of time floats
		order = array of number order of 3-hour blocks

		-----------------------------------------------------------------
	"""
	# lists to be populated
	timestamp, variation, k_index, order = [], [], [], []

	# start of the day in seconds
	day_seconds = int(minute_time[0])-int(minute_time[0])%(24*3600)

	#loop over minute_array and sort them according to 3-hour block
	start = 0
	hour_block1 = int((minute_time[0] - day_seconds)/(3*60*60))
	for index, value in enumerate(minute_time):
		hour_block2 = int((value - day_seconds)/(3*60*60))

		# if hr1 and hr2 and not equal, we have entered a new 3-hr block
		if hour_block2 != hour_block1:
			try:
				varx = max(minute_bx[start:index-1]) - min(minute_bx[start:index-1])
				vary = max(minute_by[start:index-1]) - min(minute_by[start:index-1])
	
				# append max variation for that block
				variation.append(max(varx, vary))	
				timestamp.append(day_seconds + (hour_block1*3*60*60))
				order.append(hour_block1+1)
	
				hour_block1 = hour_block2
				start = index
			except:
				continue

	# add last entry
	varx = max(minute_bx[start:-1]) - min(minute_bx[start:-1])
	vary = max(minute_by[start:-1]) - min(minute_by[start:-1])
	variation.append(max(varx, vary))
	timestamp.append(day_seconds + (hour_block1*3*60*60))
	order.append(hour_block1+1)

	# now to use these variations to calculate the k-index value
	niemegk = np.array([500, 330, 200, 120, 70, 40, 20, 10, 5, 0])	# reference
	thresh = niemegk * k9/500.0	

	k_index = []		# k_index list to be populated
	for i in variation:
		for index, j in enumerate(thresh):
			if i >= j:
				z = 9-index
				if z == 0:
					z = 0.25
				k_index.append(z)
				break

	return np.array(k_index), np.array(timestamp), np.array(order)

#-------------------------------------------------------------------------------

def FMISmooth(minute_time, minute_bx, minute_by, k_index, k_time):
	"""Calculates the Solar Regular (Sr) curve according to the FMI Method
		
		Detailed instructions can be found here:
		http://swans.meteo.be/sites/default/files/documentation/TN-RMI-2010-01_K-LOGIC.pdf

		Parameters
		-----------
		minute_time = array of minute data in float
		minute_bx, minute_by = arrays of horizontal magnetic data
		k_index = list of k_index values
		k_time = list of float times corresponding to k_index input

		Returns
		-----------
		clean_time = array of time floats
		Srx, Sry = arrays of NS and EW solar regular curves

		-----------------------------------------------------------------
	"""
	extra_time = [120, 120, 120, 60, 60, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
							  60, 60, 60, 120, 120, 120]

	start_time = minute_time[0]
	start_day = start_time - start_time%(24*60*60)
	end_time = minute_time[-1]

	# number of hours to look at:
	hour_number = int((end_time - start_day)/3600) + 1

	hour_indices =[[] for i in range(hour_number)]

	k_index_thing, k_time_thing = [], []
	for i, j in zip(k_index, k_time):
		if i == 0.25:
			k_index_thing.extend((0, 0, 0))
		else:
			k_index_thing.extend((i, i, i))

		l = j*3 
		k_time_thing.extend((l-3, l-2, l-1))

	master = np.zeros((hour_number, 7))
	blank =[]
	for i in range(hour_number):
		try:
			k = k_index_thing[k_time_thing.index(i)]
			n = k**3.3
		except:
			n = 0
		m = extra_time[i%24]

		master[i][0] = i	# hour number
	

		avg_time =start_day + (i*60*60) + (30*60) # half an hour time
		blank.append(avg_time)

		master[i][1] = avg_time

		start_time = avg_time - (30*60) - ((n + m)*60)
		master[i][2] = start_time

		end_time = avg_time + (30*60) + ((n + m)*60)
		master[i][3] = end_time


	for index, value in enumerate(minute_time):

		for i, v in enumerate(master):
			start = master[i][2]
			end = master[i][3]

			if start <= value <= end:
				master[i][4] += minute_bx[index]
				master[i][5] += minute_by[index]
				master[i][6] += 1
			

	nz = (master == 0).sum(1)
	q = master[nz == 0, :]

	smoothed_time = master[:,1]
	smoothed_bx = master[:,4]/master[:,6]
	smoothed_by = master[:,5]/master[:,6]

	clean_time, clean_bx, clean_by = [], [], []
	for t, x, y in zip(smoothed_time,smoothed_bx, smoothed_by):
		if (np.isnan(x) == False) and (np.isnan(y) == False):
			clean_time.append(t)
			clean_bx.append(x)
			clean_by.append(y)


	smt, Srx, Sry = [], [], []
	for i in clean_time:
         x = np.interp(i, clean_time, clean_bx)
         y = np.interp(i, clean_time, clean_by)
         
         Srx.append(x)
         Sry.append(y)
         

	return clean_time, Srx, Sry

#-------------------------------------------------------------------------------

def SrSmooth(clean_time, Srx, Sry, z):
	"""Smooths the solar regular curves

		Parameters
		-----------
		clean_time = array of FMI cleaned time
		Srx, Sry = arrays of Solar regular curves
		z = smoothing order (<=5)

		Returns
		-----------
		smooth_time = array of smooth time floats
		smooth_Srx, smooth_Sry = smoothed solar regular curves

		-----------------------------------------------------------------
	"""

	x = clean_time
	smooth_time = np.linspace(x[0], x[-1], 10*len(x))

	y1 = Srx
	ius1 = InterpolatedUnivariateSpline(x, y1, k=z)
	smooth_Srx = ius1(smooth_time)

	y2 = Sry
	ius2 = InterpolatedUnivariateSpline(x, y2, k=z)
	smooth_Sry = ius2(smooth_time)

	return smooth_time, smooth_Srx, smooth_Sry


#-------------------------------------------------------------------------------

def Subtracted(minute_time, smooth_time, minute_bx, smooth_Srx, minute_by, smooth_Sry):
	"""Subtracts the smoothed solar regular curves from the minute binned data.

		Parameters
		-----------
		minute_time = array of minute binned time floats
		smooth_time = array of solar regular time floats
		minute_bx, minute_by = arrays of minute binned horizontal magnetic floats
		smooth_Srx, smooth_Sry = arrays of smoothed solar regular curves

		Returns
		-----------
		subtracted_bx, subtracted_by = arrays of subtracted magnetic floats

		-----------------------------------------------------------------
	"""

	subtracted_bx, subtracted_by = [], []
	for index, value in enumerate(minute_time):
		x = np.interp(value, smooth_time, smooth_Srx)
		y = np.interp(value, smooth_time, smooth_Sry)
		subtracted_bx.append(minute_bx[index]-x)
		subtracted_by.append(minute_by[index]-y)

	return subtracted_bx, subtracted_by

#-------------------------------------------------------------------------------

def KIndexSuperCalc(minute_time, minute_bx, minute_by, n):
	"""Calculates the K-index using the FMI method.

		Details of the method used can be found here:
		http://swans.meteo.be/sites/default/files/documentation/TN-RMI-2010-01_K-LOGIC.pdf

		NOTE: This algorithm can take from 1 to 4 days of data. However many days of data are used,
		it is safer to ignore the values of the first days output.

		This is because the first day in the output are calculated only to inform the subsequent
		values.

		For example, if you input 4 days of data, only the final 3 days of the output are deemed
		"correct".

		Parameters
		-----------
		minute_time = array of minute binned time floats
		minute_bx, minute_by = arrays of minute binned horizontal magnetic components
		n = upper threshold for K9 event in nT

		Returns
		-----------
		k_index3 = array of K-index values
		k_time3 = array of K-index times (in float)
		k_order3 = order of K-index values (i.e., which 3-hour block they belong to)

		-----------------------------------------------------------------
	"""

	# First K-Index
	k_index1, k_time1, order1 = KIndex(minute_time, minute_bx, minute_by, n)

	# Second K-Index
	smoothed_time, smoothed_bx, smoothed_by = FMISmooth(minute_time, minute_bx,
	minute_by, k_index1, k_time1)

	smooth_time1, smooth_bx1, smooth_by1 = SrSmooth(smoothed_time, smoothed_bx, smoothed_by, 3)

	subtracted_bx1, subtracted_by1 = Subtracted(minute_time, smooth_time1,
	minute_bx, smooth_bx1, minute_by, smooth_by1)

	k_index2, k_time2, order2 = KIndex(minute_time, subtracted_bx1, subtracted_by1, n)

	# Third K-Index
	smoothed_time2, smoothed_bx2, smoothed_by2 = FMISmooth(minute_time, minute_bx, 
											minute_by, k_index2, k_time2)

	smooth_time2, smooth_bx2, smooth_by2 = SrSmooth(smoothed_time2, smoothed_bx2, smoothed_by2, 3)

	subtracted_bx2, subtracted_by2 = Subtracted(minute_time, smooth_time2,
	minute_bx, smooth_bx2, minute_by, smooth_by2)

	k_index3, k_time3, order3 = KIndex(minute_time, subtracted_bx2, subtracted_by2, n)

	return k_index3, k_time3, order3

#-------------------------------------------------------------------------------

def KIndexBarColor(k_index, barlist):
	"""Colours yer k-index plots so it looks nice.

	EXAMPLE:
	barlist = ax3.bar(k_timestamps, k_index, width = 0.124)
	colored(k_index)
	"""

	for i in np.arange(0, len(k_index), 1):
		if k_index[i] >= 8:
			barlist[i].set_color('deeppink')
			barlist[i].set_edgecolor('k')
			continue
		if k_index[i] >= 6:
			barlist[i].set_color('r')
			barlist[i].set_edgecolor('k')			
			continue
		if k_index[i] >= 5:
			barlist[i].set_color('orange')
			barlist[i].set_edgecolor('k')			
			continue
		if k_index[i] >= 4:
			barlist[i].set_color('g')
			barlist[i].set_edgecolor('k')			
			continue
		if k_index[i] >= 2:
			barlist[i].set_color('c')
			barlist[i].set_edgecolor('k')			
			continue
		if k_index[i] >= 0:
			barlist[i].set_color('b')
			barlist[i].set_edgecolor('k')			
			continue

#-------------------------------------------------------------------------------

def KIndexPlotter(k_index, k_time, n):
	"""Plots K-Index values in a nice Bar Plot

		NOTE: Plots 3 days of data.

		Parameters
		-----------
		k_index = list of k_index values
		k_time = list of float times corresponding to k_index input

		Returns
		-----------
		figure(n) = Plot of data
		-----------------------------------------------------------------
	"""

	fig = plt.figure(n)

	start_day = k_time[0] - k_time[0]%(24*60*60)

	start = Float2Time(start_day)
	middle1 = start + datetime.timedelta(1)
	middle2 = start + datetime.timedelta(2)
	end = start + datetime.timedelta(3)

	plt.clf()

	barlist = plt.bar(Float2Time(k_time), k_index, width = 0.124, edgecolor = "black")
	KIndexBarColor(k_index, barlist)

	for i in range(4):
		plt.axvline(Float2Time(start_day + (i*24*60*60)), lw = 2, color = "black")

	plt.ylim([0, 9])
	plt.ylabel("K-Index", fontsize = 24)

	plt.grid(True)

	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:00'))
	plt.gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,6)))
	
	plt.xlabel('{}-{}-{}                         {}-{}-{}                         {}-{}-{}'.format(start.day, start.month, start.year, middle1.day, middle1.month, middle1.year, middle2.day, middle2.month, middle2.year), fontsize = 14)

	plt.show()


	return fig

#-------------------------------------------------------------------------------

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
	return np.sqrt(np.array(minute_bx)**2 + np.array(minute_by)**2)

#-------------------------------------------------------------------------------

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
	new_time = np.arange(minute_time[0], minute_time[-1], 60)

	for value_array in args:
		a = np.interp(new_time, minute_time, value_array)
		sloped = np.gradient(a)
		sloped = np.concatenate((sloped, np.array([sloped[-1]])))

		output.append(sloped)

	return output

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

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
	middle2 = start + datetime.timedelta(2)
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

	plt.xlim([start, end])

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

	plt.xlim([start, end])

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

#-------------------------------------------------------------------------------

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
	middle2 = start + datetime.timedelta(2)
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

################################################################################
################################################################################
