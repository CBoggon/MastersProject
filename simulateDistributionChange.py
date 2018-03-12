#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit
from scipy.stats import moment
import sys

# Schulz distribution
def schulz(self, v, z, z_factorial, v_bar):
    p = [z_factorial, z, v_bar];
    #return p[0]*(v/p[2])**p[1]*np.exp(-(p[1]+1)*(v/p[2]))
    return (v**p[1]/p[0])*((p[1]+1)/p[2])**(p[1]+1)*np.exp(-(v/p[2])*(p[1]+1))


##################################################################

#Create Data set of Schulz distribution:
data = 




