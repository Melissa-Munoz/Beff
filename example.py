import numpy as np
import matplotlib.pyplot as plt
from Beff import Beff


#Magnetic properties of the star
#-------------------------------------------------------------------------------
#Bd, dipole magnetic field strenght (input same unit as desired output)
Bd = 2500.0


#Geometric angles
#-------------------------------------------------------------------------------
#Inclination angle, inc, in degrees
#Magnetic obliquity, beta, in degrees
inc = 30.
beta = 60.


#Extra parameters
#-------------------------------------------------------------------------------
#Dipole offset, (x0,y0,z0), in units of stellar radius (between 0 and 1)
#Limbd darkening parameter, u 
x0 = 0. 
y0 = 0.
z0 = 0.
u=0.3

#Calling ADM
#-------------------------------------------------------------------------------
phi = np.linspace(0.,1.,50) #rotational phase
out = Beff(phi,inc,beta,Bd,x0,y0,z0,u)


#Plotting phased longitudinal field curve
#-------------------------------------------------------------------------------
plt.figure(figsize=(9,6))
plt.plot(phi,out[0],'k')
plt.plot(phi+1,out[0],'k')
plt.plot(phi-1,out[0],'k')
plt.xlabel('Rotational phase',fontsize=14)
plt.ylabel('Longitudinal magnetic field',fontsize=14)
plt.xlim([-0.5,1.5])
plt.show()

plt.figure(figsize=(9,6))
plt.plot(phi,out[1],'k')
plt.plot(phi+1,out[1],'k')
plt.plot(phi-1,out[1],'k')
plt.xlabel('Rotational phase',fontsize=14)
plt.ylabel('Surface magnetic field',fontsize=14)
plt.xlim([-0.5,1.5])
plt.show()

