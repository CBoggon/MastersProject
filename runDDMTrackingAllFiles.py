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

plt.close('all');

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
maxVids = 2;   #Define max number of videos to analyse if the last videos have very little to analyse in them.


# Choose whether to use manually added times or calculated times in script.
manualTimes = 1;
t0 = '13:53:28';
#time0 = np.array([0., 41., 491, 607, 733, 859, 985, 1111, 1237, 1363, 1489, 1615, 1741, 1867, 1993, 2119]);
#time1 = np.array([238., 279., 523, 649, 775, 901, 1027, 1153, 1279, 1405, 1531, 1657, 1783, 1909, 2035, 2161]);
#time2 = np.array([392., 433., 565, 691, 817, 943, 1069, 1195, 1321, 1447, 1573, 1699, 1825, 1951, 2077, 2203]);

time0 = np.array([0., 41.]);
time1 = np.array([238., 279.]);
time2 = np.array([392., 433.]);

# Choose whether to plot average velocity from schulz fit to data or from mean of data:
plotCurveFittedData = 0;

##### Declare input file directory #####

## MAC ###

#fileDir='../../../../../../../Volumes/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/';
#fileDir='../../../../../../../Volumes/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140137-AsImageSequences/';

## UBUNTU ###
#fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/';
#fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140137-AsImageSequences/';
fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/';

trackingFile = '/filterTracks2DtOutput/tracks_fixed.dat';



##### Create output file directory where tracking plots will be saved #####
outputSaveFileDir = fileDir+'trackingOutput/';

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

fullTraj_velocityPos00 = [];  #np.zeros(len(os.listdir(fileDir)));
fullTraj_velocityPos01 = [];
fullTraj_velocityPos02 = [];

timePos00 = [];
timePos01 = [];
timePos02 = [];

file00Counter = 0;
file01Counter = 0;
file02Counter = 0;

timeCounterPos00 = 7*60;
timeCounterPos01 = 8*60;
timeCounterPos02 = 9*60;


measurementList = sorted(os.listdir(fileDir));
for measurement in measurementList:
    # Skip folders that do not contain image sequences
    if (measurement[len(measurement)-16:len(measurement)] != 'AsImageSequences'):
        continue;
    
    else:
        fileList = sorted(os.listdir(fileDir+measurement));
        
        for folder in fileList:
            #Ensure only the tracking files are read
            if (folder[0:6] != 'Output'):
                continue;
            
            else:

		if (folder[7:12] == 'Pos00' and file00Counter >= maxVids):
                    #Skip video if already calculated up to maxVids for pos00.
		    continue;
		
		if (folder[7:12] == 'Pos01' and file01Counter >= maxVids):
                    #Skip video if already calculated up to maxVids for pos01.
		    continue;
		
		if (folder[7:12] == 'Pos02' and file02Counter >= maxVids):
                    #Skip video if already calculated up to maxVids for pos02.
		    continue;

                print '\n######## Calculating average v for folder: '+measurement+'/'+folder+' #########\n'

                #Read in tracking data
                BIGLIST, numberOfFrames = readTrackingFile(fileDir+measurement+'/'+folder+trackingFile);
                
                #### Calculate velocity of particles
                AvVelocityArray, velocityArray, displacementArray, k_exponentArray = A.calcAverageVelocitiesForAllTraj(A, BIGLIST);
                
                # Calculate average velocity in micrometers/second
                AvVelocityArray_micrometers = pixelsToMicrons*AvVelocityArray;

                # Fit Schulz distribution to data:
                #fit_params, sigma = A.plotHistogramWithCurveFit(A, AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)', fit='schulz');
                #v_bar = fit_params[2];
                

                if (folder[7:12] == 'Pos00'):
                    
                    if (plotCurveFittedData == 1):
			## Fit Schulz distribution to data:
			fit_params, sigma = A.plotHistogramWithCurveFit(A, AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)', fit='schulz');
			v_bar = fit_params[2];
                        velocityPos00.append(v_bar);
                        velocityErrorPos00.append(sigma/2);     # sigma/2 as plt.errorbar plots magitude of y_error on each side of data point.
                    
                    else:
                        velocityPos00.append(np.mean(AvVelocityArray_micrometers));
                        velocityErrorPos00.append(np.std(AvVelocityArray_micrometers)/2);
			
                    fullTraj_velocityPos00.append(AvVelocityArray_micrometers);
                    timePos00.append(timeCounterPos00);
                    file00Counter = file00Counter + 1;
                    timeCounterPos00 = timeCounterPos00+lengthOfVideo;
                
                elif (folder[7:12] == 'Pos01'):

                    if (plotCurveFittedData == 1):
			## Fit Schulz distribution to data:
			fit_params, sigma = A.plotHistogramWithCurveFit(A, AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)', fit='schulz');
			v_bar = fit_params[2];
                        velocityPos01.append(v_bar);
                        velocityErrorPos01.append(sigma/2);     # sigma/2 as plt.errorbar plots magitude of y_error on each side of data point.
                    
                    else:
                        velocityPos01.append(np.mean(AvVelocityArray_micrometers));
                        velocityErrorPos01.append(np.std(AvVelocityArray_micrometers)/2);
                    
                    fullTraj_velocityPos01.append(AvVelocityArray_micrometers);
                    timePos01.append(timeCounterPos01);
                    file01Counter = file01Counter + 1;
                    timeCounterPos01 = timeCounterPos01+lengthOfVideo;
                
                elif (folder[7:12] == 'Pos02'):
                    
                    if (plotCurveFittedData == 1):
			## Fit Schulz distribution to data:
			fit_params, sigma = A.plotHistogramWithCurveFit(A, AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)', fit='schulz');
			v_bar = fit_params[2];
                        velocityPos02.append(v_bar);
                        velocityErrorPos02.append(sigma/2);     # sigma/2 as plt.errorbar plots magitude of y_error on each side of data point.
                    
                    else:
                        velocityPos02.append(np.mean(AvVelocityArray_micrometers));
                        velocityErrorPos02.append(np.std(AvVelocityArray_micrometers)/2);
                    
		    fullTraj_velocityPos02.append(AvVelocityArray_micrometers);
                    timePos02.append(timeCounterPos02);
                    file02Counter = file02Counter + 1;
                    timeCounterPos02 = timeCounterPos02+lengthOfVideo;
                
                else:
                    print 'ERROR: filename not recognised. Dont know where to output velocity array.'

                #outputFile = outputSaveFileDir+folder[7:22];
                #plt.savefig(outputFile);
                #plt.close('all');
                


#Convert lists to numpy arrays:
velocityPos00 = np.asarray(velocityPos00);
velocityErrorPos00 = np.asarray(velocityErrorPos00);
fullTraj_velocityPos00 = np.asarray(fullTraj_velocityPos00);
timePos00 = np.asarray(timePos00);

velocityPos01 = np.asarray(velocityPos01);
velocityErrorPos01 = np.asarray(velocityErrorPos01);
fullTraj_velocityPos01 = np.asarray(fullTraj_velocityPos01);
timePos01 = np.asarray(timePos01);

velocityPos02 = np.asarray(velocityPos02);
velocityErrorPos02 = np.asarray(velocityErrorPos02);
fullTraj_velocityPos02 = np.asarray(fullTraj_velocityPos02);
timePos02 = np.asarray(timePos02);


if (manualTimes == 1):
    timePos00 = time0;
    timePos01 = time1;
    timePos02 = time2;

#Plot overlapping, normalised histograms for each position at different times to see if there is a difference in shape.
#timesToPlotHist = [0, 4, 8, 12, 16];
timesToPlotHist = [0, 1, 0, 1, 0];

A.plotNormalisedHistograms(fullTraj_velocityPos00[timesToPlotHist[0]], str(timePos00[timesToPlotHist[0]])+' s', fullTraj_velocityPos00[timesToPlotHist[1]], str(timePos00[timesToPlotHist[1]])+' s', fullTraj_velocityPos00[timesToPlotHist[2]], str(timePos00[timesToPlotHist[2]])+' s', fullTraj_velocityPos00[timesToPlotHist[3]], str(timePos00[timesToPlotHist[3]])+' s', fullTraj_velocityPos00[timesToPlotHist[4]], str(timePos00[timesToPlotHist[4]])+' s', xlbl='Average Velocity (micrometers)', saveFilename=outputSaveFileDir+'Pos00Histograms');

A.plotNormalisedHistograms(fullTraj_velocityPos01[timesToPlotHist[0]], str(timePos01[timesToPlotHist[0]])+' s', fullTraj_velocityPos01[timesToPlotHist[1]], str(timePos01[timesToPlotHist[1]])+' s', fullTraj_velocityPos01[timesToPlotHist[2]], str(timePos01[timesToPlotHist[2]])+' s', fullTraj_velocityPos01[timesToPlotHist[3]], str(timePos01[timesToPlotHist[3]])+' s', fullTraj_velocityPos01[timesToPlotHist[4]], str(timePos01[timesToPlotHist[4]])+' s', xlbl='Average Velocity (micrometers)', saveFilename=outputSaveFileDir+'Pos01Histograms');

A.plotNormalisedHistograms(fullTraj_velocityPos02[timesToPlotHist[0]], str(timePos02[timesToPlotHist[0]])+' s', fullTraj_velocityPos02[timesToPlotHist[1]], str(timePos02[timesToPlotHist[1]])+' s', fullTraj_velocityPos02[timesToPlotHist[2]], str(timePos02[timesToPlotHist[2]])+' s', fullTraj_velocityPos02[timesToPlotHist[3]], str(timePos02[timesToPlotHist[3]])+' s', fullTraj_velocityPos02[timesToPlotHist[4]], str(timePos02[timesToPlotHist[4]])+' s', xlbl='Average Velocity (micrometers)', saveFilename=outputSaveFileDir+'Pos02Histograms');


#Plot average velocity vs time with all positions (control and both pahge videos) on one plot.
A.plotDataSetsWithErrorBars(timePos00, velocityPos00, 'average Velocity Pos00', y0_error=velocityErrorPos00, x1=timePos01, y1=velocityPos01, y1_error=velocityErrorPos01, label1='average velocity Pos01', x2=timePos02, y2=velocityPos02, y2_error=velocityErrorPos02, label2='average velocity Pos02', title='Tracked Average Velocities', xlbl='Time (Seconds)', ylbl='Average Velocity (micrometers/second)');

outputFile = outputSaveFileDir+'AvVelocityVsTime';
plt.savefig(outputFile);
plt.show()

