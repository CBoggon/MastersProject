#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit
import sys
import os

from analyseDDMTracks import analyseTrajectories


### Declare Variables
NumFramesInVideo = 2000;
NumFramesToAverageOver = 3; #Average over a number of frames to reduce the random effects of diffusion on the calculate swimming velocity
minTrajLen = 20*NumFramesToAverageOver;
fps = 50;
timePerFrame = 1./fps;
pixelsToMicrons = 0.702;    # For x20 Mag
#pixelsToMicrons = 0.354;    # For x40 Mag

minStopTimeThreshold = 1*fps;       #minimum length of time a bacteria is expected to stop before lysing. (frames)
minStopVelocityThreshold = 18;      #minimum drop in average velocity that defines a lysis event (micrometers per second)
stoppingVelocityThreshold = 0.2
D = 0.34;    #Diffusion constant micrometers/second
diffusionThreshold = (1/(float(NumFramesToAverageOver)*timePerFrame))*np.sqrt(4*D*(1/pixelsToMicrons)**2*(float(NumFramesToAverageOver)*timePerFrame));     #Above this threshold, bacteria considered to still be swimming.
minSwimmingExponent = 1.5;
minTauVals = 2;

lengthOfVideo = timePerFrame*NumFramesInVideo;  #Total length of video in seconds

### Declare file directory
fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/';
#fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140137-AsImageSequences/';
trackingFile = '/filterTracks2DtOutput/tracks_fixed.dat';

### Initialise Array
A = analyseTrajectories(minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons,  minStopTimeThreshold, minStopVelocityThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold, minSwimmingExponent, minTauVals);

### Loop through all tracking files
calculatedVelocities = [];  #np.zeros(len(os.listdir(fileDir)));
calculatedVelocityError = [];

fileCounter = 0;
fileList = sorted(os.listdir(fileDir));
for folder in fileList:
    #Ensure only the tracking files are read
    if (folder[0:6] != 'Output'):
        continue;
    else:
        #Read in tracking data
        BIGLIST, numberOfFrames = readTrackingFile(folder+'trackingFile');
        
        #### Calculate velocity of particles
        AvVelocityArray, velocityArray, displacementArray = A.calcAverageVelocitiesForAllTraj(A, BIGLIST);
        
        # Calculate average velocity in micrometers/second
        AvVelocityArray_micrometers = pixelsToMicrons*AvVelocityArray;
        
        # Fit Schulz distribution to data:
        fit_params, sigma = A.plotHistogramWithCurveFit(A, AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)', fit='schulz');
        v_bar = fit_params[2];

        calculatedVelocities.append(v_bar);
        calculatedVelocityError.append(sigma);
        
        fileCounter = fileCounter+1;


timeArray = lengthOfVideo*np.arange(fileCounter);

A.plotDataSetsWithErrorBars(time, calculatedVelocities, 'average Velocity', y0_error=calculatedVelocityError, title='Tracked Average Velocities', xlbl='Time (Seconds)', ylbl='Average Velocity (micrometers/second)');


