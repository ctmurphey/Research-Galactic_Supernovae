"Calculates the maximum distance an object can be seen based on its location in the sky"
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from math import gamma
import math

#absmag = -18 #dimmest absolute magnitude of a supernova
absmag = -6
minvis = 6 # dimmest apparent magnitude (26 for LSST, 6 for visible)

extinct = absmag - minvis #total exinction along sightline needed

R_sun = 8.5 #kpc

l_max = 180 #max longitude
b_max = 45 #max latitude



def rho_dust(r, z):
    "Dust Density Function"
    R_thin = 2.9 # kpc
    h_thin = 0.095 # kpc

    rho = np.e**(-r/R_thin)*np.e**(-abs(z)/h_thin)
    rho = rho / (R_thin * (1. - np.e**(-R_sun/R_thin)))
    return rho

def dAv_dr(radius, l_rad, b_rad):
    "Extinction Rate due to Dust"
    z_cyl = radius * np.sin(b_rad)
    # solarcentric radius component in plane
    r_par = radius * np.cos(b_rad)
    # Galactocentric cylindrical radius
    R_cyl = np.sqrt(r_par**2 - 2.*R_sun*r_par*np.cos(l_rad) + R_sun**2)

    Av_gc = 30.0
    dAv_dr = Av_gc * rho_dust(R_cyl,z_cyl)
    
    return dAv_dr

def distfinder(l,b): #accepts degrees
    "Finds Distance Based on Magnitudes"
    b_rad = math.radians(b)
    l_rad = math.radians(l)

    startpoint = 0. #where to start integral (kpc)
    endpoint   = 0. #where to end integral (kpc)
    
    newapp = absmag #creating apparent magnitude
    maxvis = 0. #maximum visible distace (pc)
    blocked = 0. # sum of total light lost to extinction, to keep from having to integrate over and over

    while newapp < minvis: #while object is still visible in sky...
        if endpoint > 40: #setting max distance for praticality purposes
                break
        startpoint = endpoint #shift where code will integrate
        endpoint += 0.01 #increment distance by 10pc
        Sigfunct = lambda r: dAv_dr(r, l_rad, b_rad) #so dAv_dr will intgrate correctly
        Sigma, err = integrate.quad(Sigfunct,startpoint,endpoint) #get magnitude loss due to extinction
        blocked += Sigma #add new light extinction to blocked

        distmod = 5*np.log10(endpoint*1000/10) #get magnitude loss due to distance

        newapp = distmod + blocked + absmag #apparent mag is abs mag + extinction + distance modulus
    maxvis = endpoint #i did this for intuitive purposes, it's largely redundant on the next line
    return maxvis

n_b = 9 #number of latitude points on plot
n_l = 15 #number of longitude points

lcoor = np.linspace(0, l_max, n_l) # array of longitude coordinates
bcoor = np.linspace(0, b_max, n_b) # array of latitude coordinates

vals = np.zeros((len(bcoor), len(lcoor))) #set up space for values



for i in range(n_b): #calculate max disance for each point on plot
    for j in range(n_l):
        vals[i][j] = distfinder(lcoor[j], bcoor[i])

#vals = vals * 3262 #convert to light years for presentation

levs = np.linspace(np.min(vals), np.max(vals), 201) #colorbar levels to use
#print(np.max(vals))

#Making a second plot to not have to rerun the code:
plt.contourf(lcoor, bcoor, vals, levs, cmap = 'rainbow') #plot quadrant 1
plt.contourf(-1*lcoor, bcoor, vals, levs, cmap = 'rainbow') #quadrant 2
plt.contourf(lcoor, -1*bcoor, vals, levs, cmap = 'rainbow') #quadrant 4
plt.contourf(-1*lcoor, -1*bcoor, vals, levs, cmap = 'rainbow') #quadrant 3

#making it pretty
plt.title("Maximum Visible Distance for Supergiant Stars", weight = 'bold', fontsize = 20)
plt.xlabel('l (degrees)', weight = 'bold', fontsize = 14)
plt.ylabel('b (degrees)', weight = 'bold', fontsize = 14)
clb = plt.colorbar(orientation = 'vertical')
clb.set_label('Max Distance (kpc)', rotation = 90, weight = 'bold', fontsize = 14)
plt.annotate('Created by Tanner Murphey', (100, -5.5), color = 'w', weight = 'bold') 
#ax = plt.gca()
#ax.set_aspect('equal')
#plt.savefig('lsstful.png') #I'd recommend changing most times the code runs
plt.show()