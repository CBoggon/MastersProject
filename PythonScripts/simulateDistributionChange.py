#!python
#
###############################################################################
#In analysing the velocity distributions from particle tracking in sample bulk, I noticed the velocity distributions became more symmetric as the phage infection progressed.
#EDIT: This was not used in report - see report for discussion of why this ocurred. It is due to 2D progression onto 3D image.
    
#Here is code written to demonstrate that when two distributions are iteratively convolved, from the central limit theorem, we expect the resultant distribution to look more and more like a gaussian.
#
###############################################################################


#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.signal import convolve
from scipy.optimize import curve_fit
from scipy.stats import moment, planck
import sys

# Schulz distribution
def schulz(v, z, z_factorial, v_bar):
    p = [z_factorial, z, v_bar];
    #return p[0]*(v/p[2])**p[1]*np.exp(-(p[1]+1)*(v/p[2]))
    return (v**p[1]/p[0])*((p[1]+1)/p[2])**(p[1]+1)*np.exp(-(v/p[2])*(p[1]+1))
    #return np.exp(-(v/p[2])*(p[1]+1))
    #return (v**p[1]/p[0])*((p[1]+1)/p[2])**(p[1]+1)

#Gaussian function
def gaussian(x, mu, sigma):
    #variance = np.square(sigma)
    #return np.exp(-np.square(x-mu)/2*variance)/(np.sqrt(2*np.pi*variance))
    return np.exp(-np.power((x - mu)/sigma, 2.)/2.)

# Exponential distribution
def exponent(v, A, v0):
    return A*np.exp(-(v/v0))

# Plot distribution of velocities
def plotHistogram(data, xlbl='bins', ylbl='Frequency', xlim=np.array(None)):
    plt.figure()
    plt.hist(data[:], bins=200)
    plt.xlabel(xlbl);
    plt.ylabel(ylbl);
    if (xlim.all() != None):
        plt.xlim(xlim[0], xlim[1]);

    #plt.savefig(outputFile);
    #plt.show()
    return

def plotNormalisedHistogram(data, xlbl='bins', ylbl='Frequency', xlim=np.array(None)):
    plt.rc('font', family='serif', size=12); 
    binNum = 100;
    weights = np.ones_like(data)/float(len(data));
    y_data,binEdges = np.histogram(data[:], weights=weights, bins=binNum);
    bincenters_data = 0.5*(binEdges[1:]+binEdges[:-1]);
    plt.plot(bincenters_data, y_data, 'x-')

    return

def plotLine(x, y, lbl, ylbl, xlbl):
    plt.plot(x,y, label=lbl)
    #plt.ylabel('gaussian distribution')
    plt.ylabel('Normalised Frequency')
    plt.xlabel('Normalised Velocity (Arbitrary Units)')

    return


##################################################################
'''
We assume control data points are drawn from a schulz distribution.
But for phages, there is some distribution associated with the infection mechanism so we draw the phages data points from the convolution of the schulz distribution and another distribution.
So here we plot the convolution of schulz with several distributions and show that they all tend towards something more symmetric and normal.
'''

#### Declare a few variables #########
NumIterations = 10;
skewnessArray = np.zeros(NumIterations+1);

### Declare first distribution  ###############
#z = 2.; z_factorial = 2.; v_bar = 10.;
#z = 10.; z_factorial = 20.; v_bar = 40.;
z = 10.; z_factorial = 20.; v_bar = 30.;
M = 1.;
v = np.arange(0,100,.1);
v_norm = v;
#v_norm = v/(v[len(v)-1]);
f = schulz(v, z, z_factorial, v_bar); 
f_norm = f/np.sum(f);
#f_norm = f/len(f);
f_skew = moment(f_norm, moment=3)/(moment(f_norm, moment=2)**(3./2.));

plt.figure(figsize=(10, 6))
plt.rc('font', family='serif', size=11);
plotLine(v_norm,f_norm, 'Original Dist: skew='+str(round(f_skew, 3)), 'Distribution', 'Velocity')
#plt.plot(v,f, label='skew = '+str(round(skewnessArray[0], 3)))
#plt.ylabel('gaussian distribution')
#plt.ylabel('Schulz Distribution')
#plt.xlabel('Velocity')
#plt.legend(loc='upper right')



### Declare second distribution as another schulz distribution  ###################
#z1 = 9.; z1_factorial = 20.; v_bar1 = 10.;
#z1 = 11.; z1_factorial = 40.; v_bar1 = 50.;

z1 = 10.; z1_factorial = 20.; v_bar1 = 60.;
#v1 = np.arange(0, 100, .01);
v1_norm = v;
#v1_norm = v1/np.mean(v1);
f1 = schulz(v, z1, z1_factorial, v_bar1); 
#f1_norm = f1/np.mean(f1);
f1_norm = f1/len(f1);
f1_skew = moment(f1_norm, moment=3)/(moment(f1_norm, moment=2)**(3./2.));
plotLine(v1_norm,f1_norm, 'Schulz: skew='+str(round(f1_skew, 3)), 'Distribution', 'Velocity')


### Declare third distribution as a gaussian distribution  ###################
#v2 = np.arange(0, 100, .01);
#v2_norm = v2/np.mean(v2);
v2_norm = v;
mean = 30.; sigma = 2.0;
f2 = gaussian(v, mean, sigma);
#f2_norm = f2/np.mean(f2);
#f2_norm = f2;
f2_norm = f2/len(f2);
f2_skew = moment(f2_norm, moment=3)/(moment(f2_norm, moment=2)**(3./2.));
plotLine(v2_norm,f2_norm, 'Gaussian: skew='+str(round(f2_skew, 3)), 'Distribution', 'Velocity')



### Declare fourth distribution as an exponential distribution  ###################
#v3 = np.arange(0, 100, .01);
#v3_norm = v3/np.mean(v3);
v3_norm = v;
A = 1000.; v0 = 10.;
f3 = exponent(v, A, v0);
#f3_norm = f3/np.mean(f3);
f3_norm = f3/len(f3);
f3_skew = moment(f3_norm, moment=3)/(moment(f3_norm, moment=2)**(3./2.));
plotLine(v3_norm,f3_norm, 'Exponential: skew='+str(round(f3_skew, 3)), 'Distribution', 'Velocity')


### Declare fourth distribution as an delta function ###################
#v4 = np.arange(0, 100, .01);
#v4 = np.array([35. for i in range(0,len(v))]);
#v4_norm = v4/np.mean(v4);
v4_norm = v_norm;
mean = 35.; sigma = 5.;
f4 = gaussian(v, mean, sigma);
#f4 = np.arange(len(f));
#f4_norm = f4/np.mean(f4);
f4_norm = f4/len(f4);
f4_skew = moment(f4_norm, moment=3)/(moment(f4_norm, moment=2)**(3./2.));
plotLine(v4_norm,f4_norm, 'Delta: skew='+str(round(f4_skew, 3)), 'Distribution', 'Velocity')

plt.legend(loc='upper right')


## Plot delta function
plt.figure(figsize=(10, 6))
plt.rc('font', family='serif', size=11);
plotLine(v_norm,f_norm, 'Original Dist: skew='+str(round(f_skew, 3)), 'Distribution', 'Velocity')
#plotLine(v4_norm,f4, 'Delta: skew='+str(round(f4_skew, 3)), 'Distribution', 'Velocity')
plt.legend(loc='upper right')

#f1_norm = f1/np.mean(f1);
#skewness1 = moment(f1_norm, moment=3)/(moment(f1_norm, moment=2)**(3./2.));

#plt.plot(v1,f1, label='skew 2 = '+str(round(skewness1, 3)))
##plt.ylabel('gaussian distribution')
#plt.ylabel('Schulz Distribution 2')
#plt.xlabel('Velocity')
#plt.legend(loc='upper right')


###### Convolve distributions ####################

#tempArray = np.zeros(len(f));
#f_new = f;
#f1_FT = np.fft.fft(f1);
#for i in range(0, NumIterations):
#    f_new = f1_FT*np.fft.fft(f_new);
#    f_new = np.fft.ifft(f_new);
#    f_new_norm = np.real(f_new)/np.mean(np.real(f_new));
#    skewnessArray[i+1] = moment(f_new_norm, moment=3)/(moment(f_new_norm, moment=2)**(3./2.));
#
#f_new = f;
#for i in range(0, NumIterations):
#    #f_new = convolve(f_new, f1, mode='same')
#    f_new = np.convolve(f_new, f1, mode='same')
#    f_new_norm = f_new/np.mean(f_new);
#    #f_new_norm = np.real(f_new)/np.mean(np.real(f_new));
#    skewnessArray[i+1] = moment(f_new_norm, moment=3)/(moment(f_new_norm, moment=2)**(3./2.));
#

## Draw data from the convolution of data set
f1_new = np.convolve(f, f1, mode='same')
f2_new = np.convolve(f, f2, mode='same')
f3_new = np.convolve(f, f3, mode='same')


f1_new_norm = f1_new/np.mean(f1_new);
f1_new_skew = moment(f1_new_norm, moment=3)/(moment(f1_new_norm, moment=2)**(3./2.));

f2_new_norm = f2_new/np.mean(f2_new);
f2_new_skew = moment(f2_new_norm, moment=3)/(moment(f2_new_norm, moment=2)**(3./2.));

f3_new_norm = f3_new/np.mean(f3_new);
f3_new_skew = moment(f3_new_norm, moment=3)/(moment(f3_new_norm, moment=2)**(3./2.));

# Plot distribution again
#
#plt.figure(figsize=(10, 6))
#plotNormalisedHistogram(f_norm, 'Original: skew='+str(round(f_skew, 3)), 'Distribution', 'Velocity')
#plotNormalisedHistogram(f1_new_norm, 'Schulz Convolution: skew='+str(round(f1_new_skew, 3)), 'Distribution', 'Velocity')
#plotNormalisedHistogram(f2_new_norm, 'Gaussian Convolution: skew='+str(round(f2_new_skew, 3)), 'Distribution', 'Velocity')
#plotNormalisedHistogram(f3_new_norm, 'Exponential Convolution: skew='+str(round(f3_new_skew, 3)), 'Distribution', 'Velocity')
#
plt.figure(figsize=(10, 6))
plotLine(v_norm,f_norm, 'Original: skew='+str(round(f_skew, 3)), 'Distribution', 'Velocity')
plotLine(v1_norm,f1_new_norm, 'Schulz Convolution: skew='+str(round(f1_new_skew, 3)), 'Distribution', 'Velocity')
plotLine(v2_norm,f2_new_norm, 'Gaussian Convolution: skew='+str(round(f2_new_skew, 3)), 'Distribution', 'Velocity')
plotLine(v3_norm,f3_new_norm, 'Exponential Convolution: skew='+str(round(f3_new_skew, 3)), 'Distribution', 'Velocity')
plt.legend(loc='upper right')

# Now calculate skewness and plot again
#f_new_norm = f_new/np.mean(f_new);
#skewness_new = moment(f_new_norm, moment=3)/(moment(f_new_norm, moment=2)**(3./2.));


# Plot distribution again
#plt.figure()
#plt.plot(v, f_new, label='skew = '+str(round(skewnessArray[len(skewnessArray)-1], 3)))
##plt.ylabel('gaussian distribution')
#plt.ylabel('Convolved Schulz Distribution')
#plt.xlabel('Velocity')
#plt.legend(loc='upper right')
#
## Plot skewness vs time
#plt.figure()
#iterations = np.arange(NumIterations+1);
#plt.plot(iterations, skewnessArray, '^-')
##plt.ylabel('gaussian distribution')
#plt.ylabel('Pearsons moment of skewness')
#plt.xlabel('Iteration')
##plt.legend(loc='upper right')


plt.show()


