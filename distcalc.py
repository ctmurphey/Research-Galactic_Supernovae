"Calculates the maximum distance we can see stars of a certain magnitude"

import numpy as np
import matplotlib.pyplot as plt

appmag = 15 #"m" apparant magnitude, visible limit
absmag = np.arange(-8, 15, 0.01)  #"M" absolute magnitude

sgiant = np.arange(-8, -3, 0.01) # Range of supergiants magnitudes

def distvis(M):
    return 10 * 10**((appmag - M)/5) #parsecs


plt.plot(absmag, distvis(absmag)/1000)
plt.xlabel("Absolute Magnitude", weight = 'bold')
plt.xlim(15, -8)
plt.ylabel("Maximum Visible Distance (kpc)", weight = 'bold')
plt.title("Apparent Magnitude 10", weight = 'bold')
plt.grid(True)
plt.fill_between(sgiant, distvis(sgiant)/1000, color = 'purple')
plt.annotate('Tanner Murphey', (14, 2000))
plt.annotate('Supergiants', weight = 'bold', xy =(-5., 100000), xytext = (2, 14000), arrowprops = dict(facecolor='black', shrink=0.05))
plt.savefig('maxdist10.pdf')
plt.show()
