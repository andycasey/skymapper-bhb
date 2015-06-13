import os, re
import random as r
from numpy import *
import pyfits
import datetime
#from pylab import clf, plot, xlim, xlabel, ylabel, draw, savefig

# This is a BHB simulation to get estimations of:
# 	1. How far out can we see BHB stars until we get more than sigmaBHB contamination?
# 	2. From galactic centre to that distance, how many RR Lyrae stars should we find?
# 	3. Assuming RRLyraes ~ xBHB's, how many BHB stars will SkyMapper find?



class extinction():
	
	def __init__(self):

		# Open the maps
		self.sgpMap = pyfits.open('SFD_dust_4096_sgp.fits')
		self.ngpMap = pyfits.open('SFD_dust_4096_ngp.fits')
				
		if (self.ngpMap[0].header['OBJECT'] != 'E(B-V)'): print 'Warning! Northern Galactic Plane extinction map header does not match that expected. Presumed E(B-V) and got %s' % (self.ngpMap[0].header['OBJECT'],)
		if (self.sgpMap[0].header['OBJECT'] != 'E(B-V)'): print 'Warning! Southern Galactic Plane extinction map header does not match that expected. Presumed E(B-V) and got %s' % (self.sgpMap[0].header['OBJECT'],)
		
		return None
	

	def ebv(self, l, b):
		"""
		ebv(l, b)
		
		Returns the E(B-V) value for a given l and b (in radians)
		"""
		
		l, b = [float(l), float(b)]
		
		# Check which plane we are looking in (and what map)
		n = -1. if (b < 0) else +1.
		map = self.sgpMap[0] if (b < 0) else self.ngpMap[0]
		
		# Schlegel et al. Appendix C Eq's C1, C2
		x = 2048*sqrt(1 - n*sin(b))*cos(l) + 2047.5
		y = -2048*n*sqrt(1 - n*sin(b))*sin(l) + 2047.5
		
		# Correct for FITS extension frame (tested July 29th 2010 in all quadrants)
		EBV = map.data[y, x]
		
		return EBV

		
	# SkyMapper u-filter, which is represented as a Sloan u' filter
	def Au(self, l, b):
		"""
		Au(l, b)
		
		Returns the extinction magnitude for the given galactic co-ordinates
		l and b (both in radians) in the SkyMapper u-filter
		"""

		EBV = self.ebv(l, b)
		
		# A/E(B-V) taken from Schlegel et al. 1998 Table 6
		AuoverEBV = 5.155
		
		Au = AuoverEBV * EBV
		return Au
	

	# SkyMapper v-filter, which is represented as a mid-point between the
	# Landolt U and B filters as a 1st-order approximation (Bessell, priv. comm.)
	def Av(self, l, b):
		"""
		Av(l, b [,map=self.ebvMap])
		
		Returns the extinction magnitude for the given galactic co-ordinates
		l and b (both in radians) in the SkyMapper little v-filter
		"""
		
		EBV = self.ebv(l, b)
		
		# A/E(B-V) taken from Schlegel et al. 1998 Table 6
		AUoverEBV = 5.434
		ABoverEBV = 4.315
		
		# Half-way between Landolt U and B (Bessell, priv. comm.)
		AvoverEBV = (AUoverEBV - ABoverEBV)/2 + ABoverEBV
		
		Av = AvoverEBV * EBV
		return Av



class simulation():
	
	def __init__(self):
		"""
		This will specify all the variable quantities for a default simulation.
		If unspecified, the variables will be as follows.
		"""
		self.r_0 = 8.5 # kpc
		self.appBrightLimit = 14 # edible
		self.rStep = 0.1 # edible
		self.elementStep = 2e-4 # edible
		self.population = 10e3 # edible
		self.sigma = 1 # edible
		self.contaminationLimit = 0.05 # edible
		self.rMax = 200 # kpc, edible
		self.colourSep = 0.4 # dex; typical for BHB/BS separation
		self.rMin = 120 # kpc
		self.haloMetallicity = -1.5 # dex, Sirko et al (2003)
		
		
		# SkyMapper technical characteristics
		# Note : Don't change these unless you know what you're doing

		# These values were taken from the SNRcalc_2.xls spreadsheet by Stefan Keller (*responsibility hand-ball*)
		self.filters = 'uvgriz' # Cross-references the following lists
		self.CFHT = [25.77, 26.5, 26.96, 26.47, 26.24, 25.3]
		self.apetures = [0.7, 0.9, 0.9, 0.9, 0.9, 0.9]
		self.Johnson = [-0.98, 0, 0.01, -0.16, -0.40, -0.53]
		self.ZP = [21.035, 22.58, 24.24, 24.16, 23.625, 22.29]
		self.ZP = [self.ZP[i] - self.Johnson[i] for i in range(0, len(self.Johnson))] # Correct ZP with Johnson		
		self.sky = [22.1, 22.0, 22.0, 21.0, 20.0, 19.2]


		# These values are typical for SkyMapper and should not be changed
		self.seeing = 1.5 # arcsecs
		self.obsTime = 110.0 # secs
		self.epoch = 3.0
		self.pxSize = 0.5 # pixel size
		self.RON = 5.0
		self.fov = 5.7 # square degrees
		self.sigmaGlobal = 0.03 # SK, 2010
		
		# Equations
		self.optap = pi*(self.seeing * 0.75/2.0)**2
		self.npix = self.optap/(self.pxSize**2)

		self.contaminationLimitFound = False
		self.distribution = 'spherical'
		
		# Open an extinction call
		self.extinction = extinction()
		
		# Sanity-check	
		lengths = map(len, [self.CFHT, self.apetures, self.Johnson, self.ZP, self.sky])
		j = lengths[0]
		if (max(lengths) != min(lengths)):
			raise ArithmeticError, "Mis-match in SkyMapper photometric dictionaries"

	
			
	def rho(self, d):
		"""
		rho(d)
		
		Returns a spatial density model, depending on what is specified
		in self.distribution. Currently 'spherical' or 'ellipsoidal'
		available.
		
		"""
		if (self.distribution == 'spherical'):
			# log(rho) = (13.86 +/- 0.47) - (3.34 +/- 0.11)*log(R) [Galactocentric distance] (Wetterer 1991)
			rho = 10**(13.86 - 3.34*log10(d*1e3))
			err = 10**(0.47 - 0.11*log10(d*1e3))
			
		elif (self.distribution == 'ellipsoidal'):
			# log(rho) = (15.71 +/- 0.56) - (3.76 +/- 0.13)*log(a) [Galactocentric semimajor distance] (Wetterer 1991)
			rho = 10**(15.71 - 3.76*log10(d*1e3))
			err = 10**(0.56 - 0.13*log10(d*1e3))
		else:	raise TypeError, 'Spatial density distribution unknown, only spherical or ellipsoidal available'

		
		return [rho, err]
		
		
	def cona(self, a):
		"""
		cona(a)
		
		This returns the c/a(a) relationship for an ellipsoidal distribution
		along the semi-major axis a
		"""
		
		if (self.distribution == 'spherical'): raise TypeError, 'Spherical spatial distribution assumed so c/a(a) does not apply here as there is no semi-major axis'
		
		cona = []
		au = 20.0
		cona0 = 0.5
		for ar in a:
			if (ar < au):
				# Preston et al 1991 Eq. 14
				cona.append(cona0 + (1 - cona0)*(ar/au))
			elif (ar > au):
				cona.append(1.0)
			else:
				# Not described in Preston et al what to do here so we will assume;
				cona.append(1.0)
		
		return cona
	

	def fieldFrac(self, l, b, Ri, lim_SNR):
		"""
		fieldFrac(l, b)
		
		Determines the single frame fractional volume of the halo
		surveyed for a field as a function of distance. Does not
		take into account completeness or contamination.
		"""
		absU, absUerr = self.RRLabsU()
		absV, absVerr = self.RRLabsV()

		# At this particular field, we will be constrained by the bright
		# and the faint limits. However, the faint limit includes some
		# extinction, which if is unknown, assumes we can see out further
		# than what we actually can. This is accounted for in self.getContaminationLimit
		

		# Get the minimum distance from our brightness magnitude limitation
		rMin = [10**((self.appBrightLimit - absX + 5.0)/5.0)*1e-3 for absX in [absU, absV]]
		rMin = max(rMin)
		
		# Get the maximum distance from the contamination

		rMax, appU, appV, Au, Av, SNR = self.getContaminationLimit(l, b, Ri, lim_SNR)

		r = arange(rMin, rMax, self.elementStep)
		vol = (1/3.)*((radians(sqrt(self.fov))**2)*(r + self.elementStep)**3 - (radians(sqrt(self.fov))**2)*r**3)
		
		# Galactocentric co-ordinates
		# R = (r^2 + R0^2 + 2*r*R0*cos(b)*cos(l))**0.5
		R = sqrt(r**2 + self.r_0**2 - 2*r*self.r_0*cos(b)*cos(l))

		# Sort the indicies from those in heliocentric to galactocentric
		idx = argsort(R)
		R, vol = [R[idx], vol[idx]]
		
		# Bin them
		
		bR = arange(0, self.rMax, self.rStep)
		bvol = zeros(len(bR))
		bin = around(R/self.rStep)

		for k in arange(0, len(R)):
			bvol[bin[k]] = bvol[bin[k]] + vol[k]
		

		volTotal = (4.0*pi/3.0)*((bR + 0.05)**3 - (bR - 0.05)**3)
		if (self.distribution == 'ellipsoidal'): volTotal = volTotal * self.cona(bR)

		frac = bvol/volTotal
		
		return [rMax, frac, appU, appV, Au, Av, SNR]
	

	def populate(self, sigmaPhot):
		"""
		populate(sigmaPhot)
		
		Generates two samples of size self.population with Gaussian
		error photometry of sigmaPhot. Returns BHB sample, then BS
		sample.
		"""
		
		a = zeros(self.population)
		b = zeros(self.population)
		
		for i in range(0, int(self.population)):
			a[i] = gaussian_rand()*sigmaPhot - sigmaPhot/2.0
			b[i] = gaussian_rand()*sigmaPhot + (self.colourSep - sigmaPhot/2.0)
			
		
		return [a, b]
		
		
		
	def photometry(self, colour):
		"""
		photometry(colour)
		
		Returns a sigmaPhot for a given colour dictionary in the form:
		
		colour = {'u' : 19.5, 'v' : 18.2 }
		
		"""
		
		sigmaPhotColour = []
		for key, filter in enumerate(colour):
			idx = self.filters.index(filter)
			magnitude = colour[filter]

			# This code (and the equations) were checked on July 28th 2010 against SNRCalc_2.xls provided by SK

			# Sky count
			# skycount = 10^((ZP - Sky) * 0.4) * optap * obstime
			skyCount = 10**((self.ZP[idx] - self.sky[idx]) * 0.4) * self.optap * self.obsTime

			# Sky noise
			# skynoise = sqrt(skycount)
			skyNoise = sqrt(skyCount)
			
			# Object count
			# objectcount = 10^((obsmag - ZP + 0.62)/-2.5) * obstime
			objectCount = 10**((magnitude - self.ZP[idx] + 0.62)/-2.5) * self.obsTime

			# Object noise
			# objectnoise = sqrt(objectcount)
			objectNoise = sqrt(objectCount)

			# Total noise
			# totalnoise = sqrt(objectCount + skyCount + npix*ron^2)
			totalNoise = sqrt(objectCount + skyCount + self.npix * self.RON**2)

			# Sigma photometry
			sigmaPhot = totalNoise / objectCount
			
			# Signal-to-Noise ratio
			# SNR = objectCount / totalNoise
			SNR = 1 / sigmaPhot

			# sigmaPhotColor = 1/SNR
			sigmaPhotColour.append(sigmaPhot)
			
		sigmaColours = sum(sigmaPhotColour)
		
		# The following equation has not been double-checked with SK (28.07.2010)
		# sigmaPhot = (sigmaColors**4 + N*sigmaGlobal**2)
		sigmaPhot = sqrt(sigmaColours**4 + len(sigmaPhotColour)*self.sigmaGlobal**2)

		# Account for multiple epochs
		sigmaPhot = sigmaPhot / sqrt(self.epoch)
		
		return sigmaPhot



	def contaminate(self, bhb, bs):
		"""
		contamination(bhb population, bs population)
		
		Returns the fractional contamination of BHB stars by BS stars
		given two populations of photometry.
		
		"""
		
		step = 0.01
		aBHB, aContam, aBS = zeros(3)

		kz = arange(floor(min(bhb)), ceil(max(bhb)), step)
		bhbLimit = min([mean(bhb) + self.sigma*std(bhb), max(bhb)])

		for k in kz:
			u, l = [k + step/2, k - step/2]
			iBHB, iContam, iBS = zeros(3)

			for m in range(0, int(self.population)):
				if (u >= bs[m] > l):
					iBHB += 1
				if (u >= bhb[m] > l):
					if (bhbLimit > l):
						iContam += 1
					iBS += 1
			if (iBHB > iContam):
				iBHB = iContam

			aBHB, aContam, aBS = [x + dx * step for x, dx in zip([aBHB, aContam, aBS], [iBHB, iContam, iBS])]

		contaminated = aBHB/aContam
		
		return contaminated


	def RRLabsV(self, metallicity=None, source='Cacciara&Clementini2003'):
		"""
		RRLabsV(metallicity)
		
		Returns the RRL Absolute V magnitude for a given metallicity and
		paper source (if supplied).
		"""
		
		# Assume a mean halo metallicity if none is specified
		metallicity = self.haloMetallicity if metallicity is None else metallicity
		
		if (source == 'Cacciara&Clementini2003'):
			
			# Cacciara & Clementini (2003)
			absV = 0.23*metallicity + 0.93
			err = 0 # SK
		else: raise TypeError, 'Unknown source for RRL absolute <V> magnitude'
		
		return [absV, err]




	def RRLabsU(self, metallicity=None, source='StefansHead'):
		"""
		RRLAbsU([metallicity = self.metallicity [,source]])
		
		Returns the RRL absolute U magnitude for a given metallicity and
		paper source (if supplied).
		"""
		
		# Assume a mean halo metallicity if none is specified
		metallicity = self.haloMetallicity if metallicity is None else metallicity
		
		if (source == 'StefansHead'):
			# Keller (from a paper, or more likely in his head)
			absV = self.RRLabsV(metallicity)
			absU = absV[0] + 1.1
			err = 0
		else: raise TypeError, 'Unknown source for RRL absolute <U> magnitude'
		
		return [absU, err]



	
	def getContaminationLimit(self, l, b, Ri, lim_SNR, contaminationLimit=None, sigma=None):
		
		sigma = self.sigma if sigma is None else sigma
		contaminationLimit = self.contaminationLimit if contaminationLimit is None else contaminationLimit
		
		absV, absVerr = self.RRLabsV()
		absU, absUerr = self.RRLabsU()
		
		appUerr, appVerr = [absUerr, absVerr]
		
		Au, Av = [self.extinction.Au(l, b), self.extinction.Av(l, b)]
		
		# Even if there is no extinction in this field, we will run the getContaminationLimit because the
		# distance (should) not change if there is a sufficient population of BHB and BS stars, and because
		# transmitting the contamination found at that distance would be cumbersome and require a lot of
		# unnecessary RAM usage.
		
		# Where are we starting from
		dr = Ri

		print '\nField location (l, b): (%3.3f, %3.3f) requires SNR of %2.3f \n' % (degrees(l), degrees(b), lim_SNR,)	
		print '	R	Uapp	Vapp	SNR'
		while True: # No code is complete without the chance of an infinite loop
			
			# Get the colour for this distance (with true extinction)
			appU, appV = [5*log10(dr * 10**3) - 5 + appX + Ax for appX, Ax in zip([absU, absV], [Au, Av])]
			colour = {'u' : appU, 'v' : appV}
			colourErr = {'u' : appU + appUerr, 'v' : appV + appVerr}

			# Establish the SNR
			sigmaPhot = self.photometry(colour)
			SNR = 1/sigmaPhot
			
			# If this is the first run through, set the old* values
			# in case there is no dust in this field.
			
			if (self.rStep > (dr % 1)) or (SNR > lim_SNR):
				print '	%2.1f	%2.2f	%2.2f	%2.3f' % (dr, colour['u'], colour['v'], SNR,)
			
			# We need a minimum SNR, so if this SNR is less than the limit we have reached the contamination zone.	
			if (SNR > lim_SNR):
				# We are stepping outwards with distance, and we have stepped slightly too far.
				return [dr, appU, appV, Au, Av, SNR]

			# Setup for the next loop
			dr -= self.rStep


				
				
			
		
	def approxContaminationLimit(self, contaminationLimit=None, sigma=None):

		sigma = self.sigma if sigma is None else sigma
		contaminationLimit = self.contaminationLimit if contaminationLimit is None else contaminationLimit
		
		# This assumes no metallicity correlation with galactocentric distance
		absV, absVerr = self.RRLabsV()
		absU, absUerr = self.RRLabsU()
		

		appUerr, appVerr = [absUerr, absVerr]
		
		r = arange(self.rMin, self.rMax, self.rStep)
		
		# These are for when the contamination limit is found, we take a step backwards.
		r_old = r[0]
		appU_old, appV_old = [5*log10(r_old * 10**3) - 5 + absX for absX in [absU, absV]]
		sigmaPhot = False
		
		data = open(str(contaminationLimit*100) + 'p' + str(sigma) + 's' + str(self.epoch) + 'e-contaminated.data', 'w')
		data.writelines('# R appU appV SNR contamination\n')
		
		# Boldly go!
		print 'Approximating magnitude limits with no extinction\n'
		print 'We require contamination less than %2.1f percent within %.1f sigma\n' % (contaminationLimit*100, sigma,)
		
		print '	R	Uapp	Vapp	SNR	Contamination	ContaminationLimit'
		for dr in r:

			
			# Apparent magnitude at this distance
			appU, appV = [5*log10(dr * 10**3) - 5 + appX for appX in [absU, absV]]

			colour = {'u' : appU, 'v' : appV}
			colourErr = {'u' : appU + appUerr, 'v' : appV + appVerr}
			# Although the colourErr is known, we do not follow this through into the contamination errors
			
			# Get the photometric uncertainties for this colour
			sigmaPhot = self.photometry(colour)
			if not sigmaPhot:
				sigmaPhot_old = sigmaPhot
			# Populate the BHB and BS samples with these photometric uncertainties
			bhb, bs = self.populate(sigmaPhot)
			
			contamination = self.contaminate(bhb, bs)
			data.writelines(' '.join(map(str, [dr, appU, appV, 1/sigmaPhot, contamination])) + '\n')
			
			print '	%2.1f	%2.2f	%2.2f	%2.3g	%2.3f		%2.3f' % (dr, colour['u'], colour['v'], 1/sigmaPhot, contamination, contaminationLimit,)
			if (contamination > self.contaminationLimit):
				data.close()
				return [r_old, appU_old, appV_old, 1/sigmaPhot_old]
				
			# Goldfish memory
			r_old, appU_old, appV_old, sigmaPhot_old = [dr, appU, appV, sigmaPhot]


	def loadFields(self, file, lcol, bcol):
		"""
		loadFields(file, l column, b column)
		
		Returns a float 2 dimensional array of the l's and b's for each
		field read from the file specified.
		"""
		
		data = loadtxt(file, type('string'), '#', None, None, 0, [lcol, bcol])
		return [map(float, d) for d in data]
		

	def totalFrac(self, fields, lim_Ri, lim_SNR, progress=True):
		"""
		totalFrac(fields, Ri [,progress=True])
		
		Returns the total fractional volume of the halo surveyed of all
		the fields given. Each field will have a limiting magnitude of
		self.appBrightLimit and self.appFaintLimit - this means totalFrac()
		cannot be run until the contamination limit (and subsequently the
		limiting magnitude) has been found.
		
		Ri is the approximate contamination limit point, which is found
		from approxContaminationLimit and required to find the individual
		field contamination limits.
		
		"""
		
		data = open('fields.data', 'w')
		data.writelines("# N l b R fracMedian appU appV Au Av SNR\n")
		data.close()
		
		contamination = []
		fractionSurveyed = []
		if progress: i, numOnScreen, nFields = [0, 50.0, len(fields)]
		
		a = datetime.datetime.now()
		for field in fields:
			
			if progress:
				i += 1
				# Get an estimate time for completion
				b = datetime.datetime.now()
				c = b - a
				
				remaining = c.seconds/float(i) * float(nFields - i)
				est = {}
				est['hours'] = floor(remaining/3600.0)
				est['minutes'] = floor((remaining % 3600.0)/60.0)
				est['seconds'] = floor(remaining % 60.0)
				
				os.system('clear')
				if (remaining > 0):
					print 'At field %s of %s		(approx %1.0f hrs %1.0f min remaining)' % (i, nFields, est['hours'], est['minutes'],)
				else:
					print 'At field %s of %s' % (i, nFields,)
				
				completed = int(floor(100*(float(i)/nFields))/(100/numOnScreen))
				completedStr = '#' * completed
				incompleteStr = '-' * int(numOnScreen - completed)

				print '[' + completedStr + incompleteStr + ']	(' + str(round(100*(float(i)/nFields), 2)) + '%)'

			if (15 > abs(float(field[1]))):
				pass
			else:

				R, frac, appU, appV, Au, Av, SNR = self.fieldFrac(radians(field[0]), radians(field[1]), lim_Ri, lim_SNR)
				
				# Write the data to file
				data = open('fields.data', 'a')
				data.writelines(' '.join(map(str, [i, field[0], field[1], R, median(frac), appU, appV, Au, Av, SNR])) + "\n")
				data.close()
				
				fractionSurveyed.append(frac)
		
		return sum(fractionSurveyed, axis=0)
		
	
	def run(self, contaminationLimit=None, sigma=None, distribution=None, progress=True, loadFrac=None):
		"""
		run([contaminationLimit[, sigma[, progress [,loadFrac=None]]]])
		
		This will run a standard BHB simulation unless other contamination
		limits are given. If loadFrac is supplied then instead of calculating the
		fraction of volume surveyed at a given galactocentric distance it loads
		in the given loadFrac filename.
		"""
		
		sigma = self.sigma if sigma is None else sigma
		contaminationLimit = self.contaminationLimit if contaminationLimit is None else contaminationLimit
		self.distribution = distribution if distribution else self.distribution
		
		
		if not loadFrac:
			# Let's get an approximate level of where to start with
			lim_R, lim_appU, lim_appV, lim_SNR = self.approxContaminationLimit(contaminationLimit, sigma)
	
			# We will use this approximation level as a starting point in each field
			fields = self.loadFields('SkyMapper.fields', 3, 4)
		
			fractionSurveyed = self.totalFrac(fields, lim_R, lim_SNR, True)
			
			# Output the fractional volume of the halo surveyed as a function of galactocentric distance
			data = open('frac-' + self.distribution + '.data', 'w')
			data.writelines('# r frac\n')
			
			R = arange(0, self.rMax, self.rStep)
			for r, frac in zip(R, fractionSurveyed):
				data.writelines(' '.join(map(str, [r, frac])) + "\n")
			
			data.close()

		
		# Variables we need to load in are R, fractionSurveyed		
		else:	R, fractionSurveyed = loadtxt(loadFrac, type(float()), '#', None, None, 0, [0, 1], True)
		
			
		if (self.distribution == 'spherical'):
			
			R = arange(0, self.rMax, self.rStep)
			rhoRRL, rhoRRLerr = self.rho(R)
		
			dNdR = fractionSurveyed * transpose(rhoRRL) * 4 * pi * R**2
			dNdRerr = fractionSurveyed * transpose(rhoRRLerr) * 4 * pi * R**2
			
			# Output the data as normal
			data = open('dist-spherical.data', 'w')
			data.writelines('# r rho dNdR dNdRerr\n')
			
			for r, rho, dndr, dndrerr in zip(R, rhoRRL, dNdR, dNdRerr):
				data.writelines(' '.join(map(str, [r, rho, dndr, dndrerr])) + '\n')
			
			data.close()
			
		
		elif (self.distribution == 'ellipsoidal'):
			
			# This is not tested yet
			a = arange(0, self.rMax, self.rStep)
			
			rhoRRL, rhoRRLerr = self.rho(a)
			cona = self.cona(a)
			
			dNda = fractionSurveyed * transpose(rhoRRL) * transpose(cona) * 4 * pi * a**2
			dNdaerr = fractionSurveyed * transpose(rhoRRLerr) * transpose(cona) * 4 * pi * a**2
		
			# Output the data as normal
			data = open('dist-ellipsoidal.data', 'w')
			data.writelines('# a rho dNda dNdaerr\n')
			
			for a, rho, dnda, dndaerr in zip(a, rhoRRL, dNda, dNdaerr):
				data.writelines(' '.join(map(str, [a, rho, dnda, dndaerr])) + '\n')
			
			data.close()
		
		else: raise TypeError, 'Unknown spatial distribution supplied: ' + self.distribution



# Gaussian population
def gaussian_rand():
	w = 1			     
	while (w >= 1):
		u1, u2 = [2*r.random()-1, 2*r.random()-1]
		w = u1**2 + u2**2
	w = sqrt((-2*log(w))/w)
	g1, g2 = [u1*w, u2*w]
	return g1

# Get uniques from a sequence and keep their index's relative
def unique(s):
	u = []
	for x in s:
		if x not in u:
			u.append(x)
	return u

sim = simulation()
sim.run(None, None, 'ellipsoidal')
