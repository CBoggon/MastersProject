#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit
import sys

from analyseDDMTracks import analyseTrajectories


##########################################################################################################
## Main Program begins here
##########################################################################################################

#Declare Variables
NumFramesToAverageOver = 3; #Average over a number of frames to reduce the random effects of diffusion on the calculate swimming velocity
minTrajLen = 10*NumFramesToAverageOver;
fps = 50;
timePerFrame = 1./fps;
pixelsToMicrons = 0.702;    # For x20 Mag
#pixelsToMicrons = 0.354;    #UNKNOWN FOR X40 MAG

minStopTimeThreshold = 1*fps;       #minimum length of time a bacteria is expected to stop before lysing. (frames)
minStopVelocityThreshold = 18;      #minimum drop in average velocity that defines a lysis event (micrometers per second)
stoppingVelocityThreshold = 0.2
D = 0.34;    #Diffusion constant micrometers/second
diffusionThreshold = (1/(float(NumFramesToAverageOver)*timePerFrame))*np.sqrt(4*D*(1/pixelsToMicrons)**2*(float(NumFramesToAverageOver)*timePerFrame));     #Above this threshold, bacteria considered to still be swimming.
minSwimmingExponent = 1.7;
minTauVals = 2;



####  Mac Filenames ####

#filename = 'testcaseTracks.dat';
#filename = 'OutputPos00_Movie0000-BrightnessContrast/TrackRodsOutput/tracks.dat';
#filename = '../Data/171201-DDM/Output-Pos01_Movie0000/trackRods2DtOutput/tracks.dat';
#filename = '../../../../../../Volumes/CBOGGONUSB/Data/Output-Pos00_Movie0004/trackRods2DtOutput/tracks.dat';

#fileDir = '../../../../../../../Volumes/CBOGGONUSB/Data/20180202-DDM/';
#filename = fileDir+'tracks_fixed.dat';
#outputFilename = fileDir+'DDMTrackingOutputData.dat';

####  Ubuntu Filenames ####

filename = '../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/Output-Pos00_Movie0001/filterTracks2DtOutput/tracks_fixed.dat';
#outputFilename = '../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/Output-Pos00_Movie0001/DDMTrackingOutput';


### Import filename from terminal ###

#filename = sys.argv[1];      #imports as string.
#outputFilename = sys.argv[2];


#Read in tracking data
BIGLIST, numberOfFrames = readTrackingFile(filename)
initialFrameNum = int(BIGLIST[0][0, 0]);
A = analyseTrajectories(minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons,  minStopTimeThreshold, minStopVelocityThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold, minSwimmingExponent, minTauVals);

#BIGLIST, PROPERTIES = readCoordinateFile(filename)
#print BIGLIST
#print numberOfFrames



#### Calculate velocity of particles
AvVelocityArray, velocityArray, displacementArray, k_exponentArray = A.calcAverageVelocitiesForAllTraj(A, BIGLIST);

# Calculate average velocity in micrometers/second
AvVelocityArray_micrometers = pixelsToMicrons*AvVelocityArray;
velocityArray_micrometers = pixelsToMicrons*velocityArray;
displacementArray_micrometers = pixelsToMicrons*displacementArray;

# Plot distribution of velocities
#A.plotHistogram(AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)');

# Plot distribution of k_exponents
A.plotHistogram(k_exponentArray, xlbl='k');
#A.plotHistogram(k_exponentArray, xlbl='k', xlim=np.array([-2., 4.]));

# Fit distribution to data:
fit_params, sigma = A.plotHistogramWithCurveFit(A, AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)', fit='schulz');
if (fit_params == 0 and sigma != 0):
    v_bar = fit_params[2];
    z = fit_params[1];


# Output average velocity
#A.writeToFile(outputFilename, data);




# Plot velocities and displacements throughout trajectories
#A.plotRandomTrajectories(A, velocityArrayByFrame_micrometers, 'time (seconds)', 'Velocity (micrometers/second)', traj2=displacementArrayByFrame_micrometers, xlbl2='time (seconds)', ylbl2='Displacement (micrometers)');

# Plot displacements throughout trajectories
#A.plotRandomTrajectories(A, displacementArrayByFrame_micrometers, 'time (seconds)', 'Displacement Over '+str(NumFramesToAverageOver)+' Frame (micrometers)');


#ID = 295;
#A.plotTrajWithSpecificID(A, BIGLIST, ID);
#A.plotTrajWithSpecificID(A, BIGLIST, ID, plotRodLength=1);

plt.show()


#import sys
#a = sys.argv[1] #imports as string.

