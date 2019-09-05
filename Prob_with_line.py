##########################################################
##
## SNSky

##
## sky map of Milky Way supernova probability per solid angle

##
## calculated for a grid of points in (l,b) Galactic coordinates
## and plotted as a countour map

##
## method of calculation:

##   given supernova rate density q(R,z) in Galactocentric cylindrical coords

##   total rate is:  R_SN = int q dV

##   rate per solid angle is:  dR/dOmega = int q dV/dOmega = int q r^2 dr

##   and so probability density = probability per solid ange is:

##      dP/dOmega = dR/dOmega / R_SN

##
## Inputs:

##   SN_type:  core collapse, Type Ia

##   SN_dist:  supernova rate density distribution

##   zoom:  show smaller field

##   N_l1q, N_b1q:  number of latitude (l) and longitude (b) points;
##                  this sets the resolution of the map and the runtime

##
## Outputs:

##   contour map of probability per solid angle

############################################################


import numpy as np

import matplotlib.pyplot as plt

import matplotlib.colorbar as colorbar

import scipy.integrate as integrate

from math import gamma


print ("Sky map of MW SN probability")

print ("imposes reflection symmetry in Galactic latitude and longitude")

print( "i.e., calculates for one quadrant and copies to the others")


## model parameters

##
## supernova type
SN_type = 'CC'
#SN_type = 'Ia'

## supernova density distribution model
SN_dist = 'spherical'

SN_dist = 'thindisk'

SN_dist = 'Adams'


# SN_dist = 'Green'

## plot region is zoomed in longitidue

zoom = True

# zoom = False

## resolution: longitude and latitude points in one quadrant

N_l1q = 180

N_b1q = 10



def dPsn_drdOmega(radius):

    # integrand of supernova probability per solid angle on sky

    # dN_sn/dr dOmega = r^2 q_sn(R,z) dr

    # where R = R(r,l,b) and z = z(r,b)

    global l_rad, b_rad

    global R_sun

    z = radius * np.sin(b_rad)

    r_par = radius * np.cos(b_rad)

    R = np.sqrt(r_par**2 - 2.*R_sun*r_par*np.cos(l_rad) + R_sun**2)


    dP_drdOmega = radius**2 * q_SN(R,z)


    return dP_drdOmega




def q_SN(R_gc,z_gc):

    # supernova density in Galactocentric cylindrical coordinates

    global R_disk, h_disk

    global SN_dist


    if (SN_dist == 'Adams'):

        q_0 = 1./(4.*np.pi*R_disk**2*h_disk)

        q = q_0 * np.exp(-R_gc/R_disk) * np.exp(-np.abs(z_gc)/h_disk)

    elif (SN_dist == 'spherical'):

        q_0 = 1./(8.*np.pi*R_disk**3)

        r_gc = np.sqrt(R_gc**2 + z_gc**2)

        q = np.exp(-r_gc/R_disk)

    elif (SN_dist == 'thindisk'):

        q_0 = 1./(8.*np.pi*R_disk*3)

        r_gc = np.sqrt(R_gc**2 + z_gc**2)

        q = np.exp(-r_gc/R_disk)

    elif (SN_dist == 'Green'):
        alpha = 1.09
        beta = 3.87
        R_sun = 8.5 # kpc
        R_sn = 295
        R_0 = R_sun / beta # kpc
        h = 0.095 # kpc
        r_gc = np.sqrt(R_gc**2 + z_gc**2)
        q = (R_sn/(4*np.pi*gamma(alpha+2)*beta**alpha*R_0**2*h))*((abs(r_gc)/R_0)**alpha)*np.exp(-abs(r_gc)/R_0)*np.exp(-abs(z_gc)/h)


    else:

        q = 0


    return q




# geometry parameters

R_sun = 8.5
 # kpc
R_cc = 2.9
 # kpc
R_ia = 2.4
 # kpc
h_cc = 0.05
# kpc
h_ia = 0.8 # kpc



if (zoom):

    l_max_deg = 90.

else:

    l_max_deg = 180.



if (SN_type == 'Ia'):

    h_disk = h_ia

    R_disk = R_ia
 # kpc

    b_max_deg = 10.

    labtext = "Type Ia Supernovae"

    if (SN_dist == 'thindisk'):

        h_disk = 0.350 # kpc
elif (SN_type == 'CC'):

    h_disk = h_cc

    R_disk = R_cc

    b_max_deg = 10.

    labtext = "Core-Collapse Supernovae"

    if (SN_dist == 'thindisk'):

        h_disk = 0.350 # kpc

else:

    print ("bad SN type!")



N_l = 2*N_l1q

N_b = 2*N_b1q


dP_dOmega = np.zeros((N_l,N_b))

lat = np.zeros((N_l,N_b))

long = np.zeros((N_l,N_b))


P_sum = 0.

b_lim = np.zeros(N_l1q)
l_lim = np.zeros(N_l1q)


for i in range(0,N_l1q):

    l_deg = i*l_max_deg/(N_l1q-1)

    l_rad = l_deg*np.pi/180.

    b_lim[i] = 10. * (1 - l_deg / 90)
    l_lim[i] = l_deg

    for j in range(0,N_b1q):


        b_deg = j*b_max_deg/(N_b1q-1)

        b_rad = b_deg*np.pi/180.


        dPdOmega, err = integrate.quad(dPsn_drdOmega,0.,5.*R_sun)


        P_sum = P_sum + dPdOmega * np.cos(b_rad)


        #dP_dOmega[i,j] = dPdOmega


        #long[i,j] = l_deg

        #lat[i,j] = b_deg


        dP_dOmega[N_l1q-i-1,N_b1q+j] = dPdOmega

        dP_dOmega[N_l1q+i,N_b1q+j] = dPdOmega

        dP_dOmega[N_l1q-i-1,N_b1q-j-1] = dPdOmega

        dP_dOmega[N_l1q+i,N_b1q-j-1] = dPdOmega


        long[N_l1q-i-1,N_b1q-j-1] = l_deg

        long[N_l1q-i-1,N_b1q+j] = l_deg

        long[N_l1q+i,N_b1q-j-1] = -l_deg

        long[N_l1q+i,N_b1q+j] = -l_deg

        lat[N_l1q-i-1,N_b1q-j-1] = -b_deg

        lat[N_l1q-i-1,N_b1q+j] = b_deg

        lat[N_l1q+i,N_b1q-j-1] = -b_deg

        lat[N_l1q+i,N_b1q+j] = b_deg



print( "SN type ", SN_type)

P_tot = 4.*P_sum*(b_max_deg/(N_b1q-1.))*(l_max_deg/(N_l1q-1.))*(np.pi/180.)**2


P_max = np.max(dP_dOmega)


print ("estimated P_tot = ",P_tot)


levs = [0.001,0.003,0.01,0.03,0.1,0.3]

levs = [0.0625,0.125,0.25,0.5,0.95]

levs = [0.01,0.03,0.1,0.3,0.99]

levs = np.linspace(0.00, 0.99, 301)

fig = plt.figure(figsize=(15.,8))

plt.title("Supernova Probability Sky Map",weight='bold')


ax1 = plt.subplot()


# cs = ax1.contour(long,lat,dP_dOmega/P_max,levs)
cs = ax1.contourf(long,lat,dP_dOmega/P_max, levs)

#ax1 = plt.subplot(111,projection="aitoff")

#ax1.grid(True)
#ax1.plot(long/(2.*np.pi),lat/(2.*np.pi),dP_dOmega/P_max,'r.')

#cs = ax1.contour(lat/(2.*np.pi),long/(2.*np.pi),dP_dOmega/P_max,levs)

#cs = ax1.contour(long,lat,dP_dOmega/P_max,levs)

#cs = ax1.imshow

ax1.set_xlabel(r"Galactic longitude $\ell$ [deg]",fontsize=20,weight="bold")

ax1.set_ylabel(r"Galactic latitude $b$ [deg]",fontsize=20,weight="bold")

ax1.text(-0.9*l_max_deg,+2.*b_max_deg/3.,labtext,fontsize=20,weight='bold', color = 'white')


# ax1.set_aspect('equal')

cbar = fig.colorbar(cs, format='%.2f')
cbar.set_label(r'Probability $P/P_{\rm max}$', fontsize = 20, weight = 'bold')

plt.plot(l_lim, b_lim, 'r-', linewidth = 3)
plt.plot(l_lim,-1 * b_lim, 'r-', linewidth = 3)
plt.plot(-1*l_lim, b_lim, 'r-', linewidth = 3)
plt.plot(-1*l_lim, -1*b_lim, 'r-', linewidth = 3)
plt.tick_params(axis='both', which='major', labelsize=20)




if zoom:
    plotname = "SNProb_" + SN_type + "_" + SN_dist + "_zoom.pdf"
else:
    plotname = "SNProb_" + SN_type + "_" + SN_dist + ".pdf"

plt.show()

fig.savefig(plotname)
