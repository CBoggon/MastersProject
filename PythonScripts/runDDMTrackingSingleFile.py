#!python
#
###############################################################################
#This code was originally designed to do the same as runDDMTrackingAllFiles but only on one file but I haven't touched it in ages so it might be full of buggy.
#
###############################################################################

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

### Declare Variables
NumFramesInVideo = 2000;
initialFrameNum = 1.0;
NumFramesToAverageOver = 3; #Average over a number of frames to reduce the random effects of diffusion on the calculate swimming velocity
minTrajLen = 30;    #10*NumFramesToAverageOver;
fps = 50;
timePerFrame = 1./fps;
pixelsToMicrons = 0.702;    # For x20 Mag
#pixelsToMicrons = 0.354;    # For x40 Mag

minStopTimeThreshold = 1*fps;       #minimum length of time a bacteria is expected to stop before lysing. (frames)
minStopVelocityThreshold = 18;      #minimum drop in average velocity that defines a lysis event (micrometers per second)
stoppingVelocityThreshold = 0.2;
D = 0.34;    #Diffusion constant micrometers/second
diffusionThreshold = (1/(float(NumFramesToAverageOver)*timePerFrame))*np.sqrt(4*D*(1/pixelsToMicrons)**2*(float(NumFramesToAverageOver)*timePerFrame));     #Above this threshold, bacteria considered to still be swimming.
minSwimmingExponent = 1.7;
minTauVals = 1;
BacteriaCounterFrame = 200.;


####  Mac Filenames ####

#filename = 'testcaseTracks.dat';
#filename = 'OutputPos00_Movie0000-BrightnessContrast/TrackRodsOutput/tracks.dat';
#filename = '../Data/171201-DDM/Output-Pos01_Movie0000/trackRods2DtOutput/tracks.dat';
#filename = '../../../../../../Volumes/CBOGGONUSB/Data/Output-Pos00_Movie0004/trackRods2DtOutput/tracks.dat';

#fileDir = '../../../../../../../Volumes/CBOGGONUSB/Data/20180202-DDM/';
#filename = fileDir+'tracks_fixed.dat';
#outputFilename = fileDir+'DDMTrackingOutputData.dat';

####  Ubuntu Filenames ####

filename = '../../../../../../../Volumes/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/Output-Pos00_Movie0000/filterTracks2DtOutput/tracks_fixed.dat';
#outputFilename = '../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/Output-Pos00_Movie0001/DDMTrackingOutput';


### Import filename from terminal ###

#filename = sys.argv[1];      #imports as string.
#outputFilename = sys.argv[2];


#Read in tracking data
BIGLIST, numberOfFrames = readTrackingFile(filename)
initialFrameNum = 1.0;
#initialFrameNum = int(BIGLIST[0][0, 0]);





############################################################################################################################


### Test average speed vs number of frames to average over to check tracking code:

#minTrajLen_Array = np.array([12, 20, 30, 50]);
#minTrajLen_Array = np.array([30, 40, 50, 60]);
#minTrajLen_Array = np.array([20, 30, 40, 60]);

#minK_Array = np.array([1.2, 1.4, 1.6, 1.8]);


minTrajLen_Array = np.array([30]);
#minTrajLen_Array = np.array([12, 20, 30, 40]);
minK_Array = np.array([1.4]);

AvVelocity_SpeedCheck = [];
AvVelocity_SpeedCheck_error = [];

NumFramesToAverageOver_Array = np.arange(12)+1;     # NB: int(minTrajLen/NumFramesToAverageOver) must be > 0 !
AvVelocity_SpeedCheckPlot = np.zeros(len(NumFramesToAverageOver_Array));
AvVelocity_SpeedCheckPlot_error = np.zeros(len(NumFramesToAverageOver_Array));

#for j in range(0, len(minTrajLen_Array)):
for j in range(0, len(minK_Array)):
    print '\n\n\n New j in loop \n\n'
    minTrajLen = minTrajLen_Array[j];
    minSwimmingExponent = minK_Array[j];
    for i in range(0, len(NumFramesToAverageOver_Array)):
        NumFramesToAverageOver = NumFramesToAverageOver_Array[i];
        #print 'NumFramesToAverageOver = '+str(NumFramesToAverageOver)+'\ni = '+str(i)+'\nj = '+str(j)+'\nminTrajLen = '+str(minTrajLen);
        
        from analyseDDMTracks import analyseTrajectories
        A = analyseTrajectories(minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons,  minStopTimeThreshold, minStopVelocityThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold, minSwimmingExponent, minTauVals, BacteriaCounterFrame);
        
        #Calculate velocities
        AvVelocityArray, velocityArray, displacementArray, k_exponentArray, NBacteriaCounter = A.calcAverageVelocitiesForAllTraj(A, BIGLIST);
        
        print AvVelocityArray;

        #Convert to micrometers
        AvVelocityArray_micrometers = pixelsToMicrons*AvVelocityArray;
        
        #Append to array
        AvVelocity_SpeedCheckPlot[i] = np.mean(AvVelocityArray_micrometers);
        AvVelocity_SpeedCheckPlot_error[i] = np.std(AvVelocityArray_micrometers)/2;
        
    AvVelocity_SpeedCheck.append(AvVelocity_SpeedCheckPlot);
    AvVelocity_SpeedCheck_error.append(AvVelocity_SpeedCheckPlot_error);

AvVelocity_SpeedCheck = np.asarray(AvVelocity_SpeedCheck);
AvVelocity_SpeedCheck_error = np.asarray(AvVelocity_SpeedCheck_error);


#A.plotDataSetsWithErrorBars(NumFramesToAverageOver_Array, AvVelocity_SpeedCheckPlot, 'plot', title='tracking speed check', xlbl='Number of frames tracking algorithm averages over', ylbl='Speed (micrometers/second)', plotLegend=False);

A.plotDataSetsWithErrorBars(NumFramesToAverageOver_Array, AvVelocity_SpeedCheckPlot, 'plot', title='tracking speed check', xlbl='Number of frames tracking algorithm averages over', ylbl='Speed (micrometers/second)', plotLegend=False);


#A.plotDataSetsWithErrorBars(NumFramesToAverageOver_Array, AvVelocity_SpeedCheck[0], 'minTrajLen = '+str(minTrajLen_Array[0]), y0_error=AvVelocity_SpeedCheck_error[0], x1=NumFramesToAverageOver_Array, y1=AvVelocity_SpeedCheck[1], y1_error=AvVelocity_SpeedCheck_error[1], label1='minTrajLen = '+str(minTrajLen_Array[1]), x2=NumFramesToAverageOver_Array, y2=AvVelocity_SpeedCheck[2], label2='minTrajLen = '+str(minTrajLen_Array[2]), y2_error=AvVelocity_SpeedCheck_error[2], x3=NumFramesToAverageOver_Array, y3=AvVelocity_SpeedCheck[3], y3_error=AvVelocity_SpeedCheck_error[3], label3='minTrajLen = '+str(minTrajLen_Array[3]), title='tracking speed check', xlbl='Number of frames tracking algorithm averages over', ylbl='Speed (micrometers/second)', plotLegend=False);


#A.plotDataSetsWithErrorBars(NumFramesToAverageOver_Array, AvVelocity_SpeedCheck[0], 'min k = '+str(minK_Array[0]), y0_error=AvVelocity_SpeedCheck_error[0], x1=NumFramesToAverageOver_Array, y1=AvVelocity_SpeedCheck[1], y1_error=AvVelocity_SpeedCheck_error[1], label1='min k = '+str(minK_Array[1]), x2=NumFramesToAverageOver_Array, y2=AvVelocity_SpeedCheck[2], label2='min k = '+str(minK_Array[2]), y2_error=AvVelocity_SpeedCheck_error[2], x3=NumFramesToAverageOver_Array, y3=AvVelocity_SpeedCheck[3], y3_error=AvVelocity_SpeedCheck_error[3], label3='min k = '+str(minK_Array[3]), title='tracking speed check', xlbl='Number of frames tracking algorithm averages over', ylbl='Speed (micrometers/second)', plotLegend=False);




############################################################################################################################



#A = analyseTrajectories(minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons,  minStopTimeThreshold, minStopVelocityThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold, minSwimmingExponent, minTauVals, BacteriaCounterFrame);
#
##BIGLIST, PROPERTIES = readCoordinateFile(filename)
##print BIGLIST
##print numberOfFrames
#
#
##### Calculate velocity of particles
##AvVelocityArray, velocityArray, displacementArray, k_exponentArray, NBacteriaCounter = A.calcAverageVelocitiesForAllTraj(A, BIGLIST);
#
## Calculate average velocity in micrometers/second
##AvVelocityArray_micrometers = pixelsToMicrons*AvVelocityArray;
##velocityArray_micrometers = pixelsToMicrons*velocityArray;
##displacementArray_micrometers = pixelsToMicrons*displacementArray;
#
## Plot distribution of velocities
##A.plotHistogram(AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)');
#
## Plot distribution of k_exponents
##A.plotHistogram(k_exponentArray, xlbl='k');
##A.plotHistogram(k_exponentArray, xlbl='k', xlim=np.array([-2., 4.]));
#
#
#
#
#ID = 295;
#A.plotTrajWithSpecificID(A, BIGLIST, ID);
###A.plotTrajWithSpecificID(A, BIGLIST, ID, plotRodLength=1);
#
#





plt.show()









############################################################################################################################

### OLD STUFF:

# Fit distribution to data:
#fit_params, sigma = A.plotHistogramWithCurveFit(A, AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)', fit='schulz');
#if (fit_params == 0 and sigma != 0):
#    v_bar = fit_params[2];
#    z = fit_params[1];



# Output average velocity
#A.writeToFile(outputFilename, data);

# Plot velocities and displacements throughout trajectories
#A.plotRandomTrajectories(A, velocityArrayByFrame_micrometers, 'time (seconds)', 'Velocity (micrometers/second)', traj2=displacementArrayByFrame_micrometers, xlbl2='time (seconds)', ylbl2='Displacement (micrometers)');

# Plot displacements throughout trajectories
#A.plotRandomTrajectories(A, displacementArrayByFrame_micrometers, 'time (seconds)', 'Displacement Over '+str(NumFramesToAverageOver)+' Frame (micrometers)');


#import sys
#a = sys.argv[1] #imports as string.

