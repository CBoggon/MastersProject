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
initialFrameNum = 1.0;
NumFramesToAverageOver = 3; #Average over a number of frames to reduce the random effects of diffusion on the calculate swimming velocity
minTrajLen = 10*NumFramesToAverageOver;
fps = 50;
timePerFrame = 1./fps;
pixelsToMicrons = 0.702;    # For x20 Mag
#pixelsToMicrons = 0.354;    # For x40 Mag

minStopTimeThreshold = 1*fps;       #minimum length of time a bacteria is expected to stop before lysing. (frames)
minStopVelocityThreshold = 18;      #minimum drop in average velocity that defines a lysis event (micrometers per second)
stoppingVelocityThreshold = 0.2
D = 0.34;    #Diffusion constant micrometers/second
diffusionThreshold = (1/(float(NumFramesToAverageOver)*timePerFrame))*np.sqrt(4*D*(1/pixelsToMicrons)**2*(float(NumFramesToAverageOver)*timePerFrame));     #Above this threshold, bacteria considered to still be swimming.
minSwimmingExponent = 1.7;
minTauVals = 2;

lengthOfVideo = timePerFrame*NumFramesInVideo;  #Total length of video in seconds
maxVids = 13;   #Define max number of videos to analyse if the last videos have very little to analyse in them.


##### Declare input file directory #####

## MAC ###

#fileDir='../../../../../../../Volumes/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/';
fileDir='../../../../../../../Volumes/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140137-AsImageSequences/';

## UBUNTU ###
#fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/';
#fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140137-AsImageSequences/';

trackingFile = '/filterTracks2DtOutput/tracks_fixed.dat';



##### Create output file directory where tracking plots will be saved #####
outputSaveFileDir = fileDir+'/trackingOutput/';

try:
    os.stat(outputSaveFileDir)
except:
    os.mkdir(outputSaveFileDir);



### Initialise Array
A = analyseTrajectories(minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons,  minStopTimeThreshold, minStopVelocityThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold, minSwimmingExponent, minTauVals);

### Loop through all tracking files
velocityPos00 = [];  #np.zeros(len(os.listdir(fileDir)));
velocityPos01 = [];
velocityPos02 = [];

velocityErrorPos00 = [];
velocityErrorPos01 = [];
velocityErrorPos02 = [];

timePos00 = [];
timePos01 = [];
timePos02 = [];

file00Counter = 0;
file01Counter = 0;
file02Counter = 0;

timeCounterPos00 = 0;
timeCounterPos01 = 0;
timeCounterPos02 = 0;

fileList = sorted(os.listdir(fileDir));
for folder in fileList:
    #Ensure only the tracking files are read
    if (folder[0:6] != 'Output'):
        continue;
    
    else:
        print '\n######## Calculating average v for folder: '+str(folder)+' #########\n'

        #Read in tracking data
        BIGLIST, numberOfFrames = readTrackingFile(fileDir+folder+trackingFile);
        
        #### Calculate velocity of particles
        AvVelocityArray, velocityArray, displacementArray, k_exponentArray = A.calcAverageVelocitiesForAllTraj(A, BIGLIST);
        
        # Calculate average velocity in micrometers/second
        AvVelocityArray_micrometers = pixelsToMicrons*AvVelocityArray;
        
        # Fit Schulz distribution to data:
        fit_params, sigma = A.plotHistogramWithCurveFit(A, AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)', fit='schulz');
	v_bar = fit_params[2];
        
        if (folder[7:12] == 'Pos00'):
            #Skip video if already calculated up to maxVids for pos00.
            if (file00Counter >= maxVids):
                continue; 
            
            velocityPos00.append(v_bar);
            velocityErrorPos00.append(sigma);
            timePos00.append(timeCounterPos00);
            file00Counter = file00Counter + 1;
            timeCounterPos00 = timeCounterPos00+lengthOfVideo;
        
        elif (folder[7:12] == 'Pos01'):
            #Skip video if already calculated up to maxVids for pos01.
            if (file01Counter >= maxVids):
                continue;

            velocityPos01.append(v_bar);
            velocityErrorPos01.append(sigma);
            timePos01.append(timeCounterPos01);
            file01Counter = file01Counter + 1;
            timeCounterPos01 = timeCounterPos01+lengthOfVideo;
        
        elif (folder[7:12] == 'Pos02'):
            #Skip video if already calculated up to maxVids for pos02.
            if (file02Counter >= maxVids):
                continue; 
            
            velocityPos02.append(v_bar);
            velocityErrorPos02.append(sigma);
            timePos02.append(timeCounterPos02);
            file02Counter = file02Counter + 1;
            timeCounterPos02 = timeCounterPos02+lengthOfVideo;
        
        else:
            print 'ERROR: filename not recognised. Dont know where to output velocity array.'

        outputFile = outputSaveFileDir+folder[7:22];
        plt.savefig(outputFile);
        plt.close('all');
        
#Convert lists to numpy arrays:
velocityPos00 = np.asarray(velocityPos00);
velocityErrorPos00 = np.asarray(velocityErrorPos00);
timePos00 = np.asarray(timePos00);

velocityPos01 = np.asarray(velocityPos01);
velocityErrorPos01 = np.asarray(velocityErrorPos01);
timePos01 = np.asarray(timePos01);

velocityPos02 = np.asarray(velocityPos02);
velocityErrorPos02 = np.asarray(velocityErrorPos02);
timePos02 = np.asarray(timePos02);


A.plotDataSetsWithErrorBars(timePos00, velocityPos00, 'average Velocity Pos00', y0_error=velocityErrorPos00, x1=timePos01, y1=velocityPos01, y1_error=velocityErrorPos01, label1='average velocity Pos01', x2=timePos02, y2=velocityPos02, y2_error=velocityErrorPos02, label2='average velocity Pos02', title='Tracked Average Velocities', xlbl='Time (Seconds)', ylbl='Average Velocity (micrometers/second)');

outputFile = outputSaveFileDir+'AvVelocityVsTime';
plt.savefig(outputFile);
plt.show()

