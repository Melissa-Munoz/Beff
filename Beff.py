#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Code to compute the longitudinal and surface magnetic field of an obliquely 
# rotating magnetic dipole.
# 
# Melissa Munoz
# Updated Dec 2020
#
# See also FLDCURV by LANDSTREET, 1986
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#-------------------------------------------------------------------------------
# Library import ---------------------------------------------------------------
#-------------------------------------------------------------------------------


import numpy as np
#import matplotlib.pyplot as plt
from scipy.integrate import simps

#-------------------------------------------------------------------------------
# Longitudinal and surface magnetic field calculation --------------------------
#-------------------------------------------------------------------------------


def Beff(phi,inc,beta,Bd,x0,y0,z0,u):

	#Length of arrays
	PH = len(phi) #Rotational phase
	N = 100 #Spatial grid size

	#Defining polar and azimuthal grids
	MU = np.linspace(0.,1.,N) 
	PHI = np.linspace(0.,1.,N)*2.*np.pi

	#Creating 2D meshgrids
	MU_grid, PHI_grid = np.meshgrid( MU, PHI, indexing='xy')
	X = np.sqrt(1.-MU_grid**2)*np.cos(PHI_grid) 
	Y = np.sqrt(1.-MU_grid**2)*np.sin(PHI_grid) 
	Z = MU_grid 

	#Variable setup
	Bl = np.zeros(PH) #empty array for longitudinal field
	Bs = np.zeros(PH) #empty array for surface field
	inc = np.radians(inc)
	beta = np.radians(beta)
	phi = 2.*np.pi *phi
	for ph in range(0,PH):
		
		#Observer's frame (X,Y,Z)
		#Inclined and rotated frame (XR,YR,ZR)
		#Magnetic frame (XM,YM,ZM)

		#Set of transformations: observer's frame -> magnetic frame
		XR=X*np.cos(inc)-Z*np.sin(inc)
		YR=Y
		ZR=X*np.sin(inc)+Z*np.cos(inc)
		XM=(XR*np.cos(phi[ph])+YR*np.sin(phi[ph]))*np.cos(beta)+ZR*np.sin(beta)
		YM=-XR*np.sin(phi[ph])+YR*np.cos(phi[ph])
		ZM=-(XR*np.cos(phi[ph])+YR*np.sin(phi[ph]))*np.sin(beta)+ZR*np.cos(beta)

		#Dipole offset (x0,y0,z0)
		x = XM-x0
		y = YM-y0
		z = ZM-z0
		r = np.sqrt(x**2 + y**2 + z**2)

		#Components of the dipolar magnetic field vector in cartesian coordinates
		Bx = 3.*Bd/2*x*z/r**5
		By = 3.*Bd/2*y*z/r**5
		Bz = Bd/2*( 3.*z**2 - r**2 )/r**5

		#Set of transformations: magnetic frame -> observer's frame
		BXR = (Bx*np.cos(beta) - Bz*np.sin(beta))*np.cos(phi[ph]) - By*np.sin(phi[ph])
		BYR = (Bx*np.cos(beta) - Bz*np.sin(beta))*np.sin(phi[ph]) + By*np.cos(phi[ph])
		BZR = Bx*np.sin(beta) + Bz*np.cos(beta)
		BX = BXR*np.cos(inc) + BZR*np.sin(inc) 
		BY = BYR
		BZ = -BXR*np.sin(inc)  + BZR*np.cos(inc) 

		#Magnitude of the field
		B = np.sqrt(BX**2 + BY**2 + BZ**2)	

		#Surface integrals over the visible portion of the star
		meanBz = simps(simps(BZ*MU_grid*(1-u+u*MU_grid),MU),PHI)
		meanBs = simps(simps(B*MU_grid*(1-u+u*MU_grid),MU),PHI)
		mean = simps(simps(MU_grid*(1-u+u*MU_grid),MU),PHI) 

		#Resulting longitudinal and surface magentic field 
		Bl[ph] = meanBz/mean 
		Bs[ph] = meanBs/mean 

	return [Bl,Bs]


