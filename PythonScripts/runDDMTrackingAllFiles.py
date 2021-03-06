#!python
#
###############################################################################
#Master program for running tracking analysis on videos taken in bulk of sample (same videos as used for DDM analysis). After running 'runTracking.sh', this program reads in all the files called 'tracks_fixed.dat', which are saved in a folder within the folder containing video as an image sequence, calculates all trajectories, filters out non swimmers and plots results as; histograms of velocity distribution for each sample, average tracked velocity vs time and number of tracked bacteria vs time.
#
#My experiment involved sequentially videoing 3 different samples for 2000 frames at 50 fps (so each video is 40s a part and each subsequent video on the same sample is 3x40s=2 mins a part). The program identifies which file belongs to which video by reading the name of the folder it is in. Pos00 = control, Pos01 = phage1, pos02 = phage2.
#
#
###############################################################################


#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
#from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit
import datetime
import sys
import os

from analyseDDMTracks import analyseTrajectories

plt.close('all');

### Declare Variables
NumFramesInVideo = 2000;
initialFrameNum = 1.0;
NumFramesToAverageOver = 3; #Average over a number of frames to reduce the random effects of diffusion on the calculate swimming velocity
minTrajLen = 30;    #3*NumFramesToAverageOver+1;  
fps = 50;
timePerFrame = 1./fps;
pixelsToMicrons = 0.702;    # For x20 Mag
#pixelsToMicrons = 0.354;    # For x40 Mag

minStopTimeThreshold = 1*fps;       #minimum length of time a bacteria is expected to stop before lysing. (frames)
minStopVelocityThreshold = 18;      #minimum drop in average velocity that defines a lysis event (micrometers per second)
stoppingVelocityThreshold = 0.2;
D = 0.34;    #Diffusion constant micrometers/second
diffusionThreshold = (1/(float(NumFramesToAverageOver)*timePerFrame))*np.sqrt(4*D*(1/pixelsToMicrons)**2*(float(NumFramesToAverageOver)*timePerFrame));     #Above this threshold, bacteria considered to still be swimming.
minSwimmingExponent = 1.5;
minTauVals = 1;
BacteriaCounterFrame = 200.;

lengthOfVideo = timePerFrame*NumFramesInVideo;  #Total length of video in seconds

## Code takes ages to run so if I want to shorten program during debugging, I just change the value of maxVids.
maxVids = 16;   #Define max number of videos to analyse if the last videos have very little to analyse in them.


# Choose whether to use manually added times or calculated times in script.
manualTimes = 1;    # 1 = use manual times, 0 use calculated times.
t0 = '13:53:28';
time0 = np.array([0., 41., 491, 607, 733, 859, 985, 1111, 1237, 1363, 1489, 1615, 1741, 1867, 1993, 2119]);
time1 = np.array([238., 279., 523, 649, 775, 901, 1027, 1153, 1279, 1405, 1531, 1657, 1783, 1909, 2035, 2161]);
time2 = np.array([392., 433., 565, 691, 817, 943, 1069, 1195, 1321, 1447, 1573, 1699, 1825, 1951, 2077, 2203]);

#time0 = np.array([0., 41., 491, 607, 733, 859, 985, 1111, 1237, 1363, 1489, 1615, 1741, 1867, 1993]);
#time1 = np.array([238., 279., 523, 649, 775, 901, 1027, 1153, 1279, 1405, 1531, 1657, 1783, 1909, 2035]);
#time2 = np.array([392., 433., 565, 691, 817, 943, 1069, 1195, 1321, 1447, 1573, 1699, 1825, 1951, 2077]);

#time0 = np.array([0., 41.]);
#time1 = np.array([238., 279.]);
#time2 = np.array([392., 433.]);

#time0 = np.array([0.]);
#time1 = np.array([238.]);
#time2 = np.array([392.]);


# Choose whether to plot average velocity from schulz fit to data or from mean of data:
plotCurveFittedData = 0;    #The schulz fit wasn't very good, so set this to 0 to plot mean of data.

##### Declare input file directory #####
# (Mac and Ubuntu have slightly different paths to external hard drives)

## MAC ###

#fileDir='../../../../../../../Volumes/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/';
#fileDir='../../../../../../../Volumes/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140137-AsImageSequences/';
fileDir='../../../../../../../Volumes/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/';

## UBUNTU ###
#fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/';
#fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140137-AsImageSequences/';
#fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/';
#fileDir='../../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/';

trackingFile = '/filterTracks2DtOutput/tracks_fixed.dat';


##### Create output file directory where tracking plots will be saved #####

#outputSaveFileDir = fileDir+'trackingOutput/';
#outputSaveFileDir = '../../Data/Results/DDM/';
outputSaveFileDir = '../../../../../../../../Volumes/CBOGGONUSB/Data/DDMResults/';
#outputSaveFileDir = '../../../../../../../media/cameron/CBOGGONUSB/Data/DDMResults/';

try:
    os.stat(outputSaveFileDir)
except:
    os.mkdir(outputSaveFileDir);



### Initialise Array
A = analyseTrajectories(minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons,  minStopTimeThreshold, minStopVelocityThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold, minSwimmingExponent, minTauVals, BacteriaCounterFrame);

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

NArrayPos00 = [];   #Number of bacteria
NArrayPos01 = [];
NArrayPos02 = [];

NArrayErrorPos00 = [];
NArrayErrorPos01 = [];
NArrayErrorPos02 = [];

timePos00 = [];
timePos01 = [];
timePos02 = [];

file00Counter = 0;
file01Counter = 0;
file02Counter = 0;

timeCounterPos00 = 7*60;
timeCounterPos01 = 8*60;
timeCounterPos02 = 9*60;

#begin searching through all files in directory.
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
                AvVelocityArray, velocityArray, displacementArray, k_exponentArray, NBacteriaCounter = A.calcAverageVelocitiesForAllTraj(A, BIGLIST);
                
                # Calculate average velocity in micrometers/second
                AvVelocityArray_micrometers = pixelsToMicrons*AvVelocityArray;
                #AvVelocityArray_micrometers = AvVelocityArray;

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
                        velocityErrorPos00.append(np.std(AvVelocityArray_micrometers)/2);   # standard deviation
		    
                    fullTraj_velocityPos00.append(AvVelocityArray_micrometers);
		    NArrayPos00.append(NBacteriaCounter);
		    NArrayErrorPos00.append(np.sqrt(NBacteriaCounter)/2);     #Assume n taken from poisson distribution
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
		    NArrayPos01.append(NBacteriaCounter);
		    NArrayErrorPos01.append(np.sqrt(NBacteriaCounter)/2);
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
		    NArrayPos02.append(NBacteriaCounter);
		    NArrayErrorPos02.append(np.sqrt(NBacteriaCounter)/2);
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
NArrayPos00 = np.asarray(NArrayPos00);
NArrayErrorPos00 = np.asarray(NArrayErrorPos00);
timePos00 = np.asarray(timePos00);

velocityPos01 = np.asarray(velocityPos01);
velocityErrorPos01 = np.asarray(velocityErrorPos01);
fullTraj_velocityPos01 = np.asarray(fullTraj_velocityPos01);
NArrayPos01 = np.asarray(NArrayPos01);
NArrayErrorPos01 = np.asarray(NArrayErrorPos01);
timePos01 = np.asarray(timePos01);

velocityPos02 = np.asarray(velocityPos02);
velocityErrorPos02 = np.asarray(velocityErrorPos02);
fullTraj_velocityPos02 = np.asarray(fullTraj_velocityPos02);
NArrayPos02 = np.asarray(NArrayPos02);
NArrayErrorPos02 = np.asarray(NArrayErrorPos02);
timePos02 = np.asarray(timePos02);


if (manualTimes == 1):
    timePos00 = time0[0:maxVids];
    timePos01 = time1[0:maxVids];
    timePos02 = time2[0:maxVids];

#Plot overlapping, normalised histograms for each position at different times to see if there is a difference in shape.
timesToPlotHist = [0, 3, 6, 9];
#timesToPlotHist = [2, 4, 6, 7, 8];
#timesToPlotHist = [0, 1, 0, 1, 0];
#timesToPlotHist = [0, 0, 0, 0, 0];

#Define labels for plots
labelPos00_0 = str(datetime.timedelta(seconds=round(timePos00[timesToPlotHist[0]],0)))[2:10];
labelPos00_1 = str(datetime.timedelta(seconds=round(timePos00[timesToPlotHist[1]],0)))[2:10];
labelPos00_2 = str(datetime.timedelta(seconds=round(timePos00[timesToPlotHist[2]],0)))[2:10];
labelPos00_3 = str(datetime.timedelta(seconds=round(timePos00[timesToPlotHist[3]],0)))[2:10];

labelPos01_0 = str(datetime.timedelta(seconds=round(timePos01[timesToPlotHist[0]],0)))[2:10];
labelPos01_1 = str(datetime.timedelta(seconds=round(timePos01[timesToPlotHist[1]],0)))[2:10];
labelPos01_2 = str(datetime.timedelta(seconds=round(timePos01[timesToPlotHist[2]],0)))[2:10];
labelPos01_3 = str(datetime.timedelta(seconds=round(timePos01[timesToPlotHist[3]],0)))[2:10];

labelPos02_0 = str(datetime.timedelta(seconds=round(timePos02[timesToPlotHist[0]],0)))[2:10];
labelPos02_1 = str(datetime.timedelta(seconds=round(timePos02[timesToPlotHist[1]],0)))[2:10];
labelPos02_2 = str(datetime.timedelta(seconds=round(timePos02[timesToPlotHist[2]],0)))[2:10];
labelPos02_3 = str(datetime.timedelta(seconds=round(timePos02[timesToPlotHist[3]],0)))[2:10];

#Convert time to minutes
timePos00 = timePos00/60;
timePos01 = timePos01/60;
timePos02 = timePos02/60;


# Plot normalised histrograms for Pos00, Pos01, pos02.
A.plotNormalisedHistograms(fullTraj_velocityPos00[timesToPlotHist[0]], labelPos00_0, fullTraj_velocityPos00[timesToPlotHist[1]], labelPos00_1, fullTraj_velocityPos00[timesToPlotHist[2]], labelPos00_2, fullTraj_velocityPos00[timesToPlotHist[3]], labelPos00_3, xlbl='Normalised Velocity (v/<v>)', plotAsLines=True, saveFilename=outputSaveFileDir+'Pos00Histograms');

A.plotNormalisedHistograms(fullTraj_velocityPos01[timesToPlotHist[0]], labelPos01_0, fullTraj_velocityPos01[timesToPlotHist[1]], labelPos01_1, fullTraj_velocityPos01[timesToPlotHist[2]], labelPos01_2, fullTraj_velocityPos01[timesToPlotHist[3]], labelPos01_3, xlbl='Normalised Velocity (v/<v>)', plotAsLines=True, saveFilename=outputSaveFileDir+'Pos01Histograms');

A.plotNormalisedHistograms(fullTraj_velocityPos02[timesToPlotHist[0]], labelPos02_0, fullTraj_velocityPos02[timesToPlotHist[1]], labelPos02_1, fullTraj_velocityPos02[timesToPlotHist[2]], labelPos02_2, fullTraj_velocityPos02[timesToPlotHist[3]], labelPos02_3, xlbl='Normalised Velocity (v/<v>)', plotAsLines=True, saveFilename=outputSaveFileDir+'Pos02Histograms');


#Plot histograms of all positions on same figure, shortly after infection 
# This is a fancy plot with all three plots on one figure.

#A.plotHistogramsInSameFig(fullTraj_velocityPos00[0], labelPos00_0, data0_1=fullTraj_velocityPos00[1], label0_1=labelPos00_1, data1_0=fullTraj_velocityPos01[0], label1_0=labelPos01_0, data1_1=fullTraj_velocityPos01[1], label1_1=labelPos01_1, data2_0=fullTraj_velocityPos02[0], label2_0=labelPos02_0, data2_1=fullTraj_velocityPos02[1], label2_1=labelPos02_1, xlbl='Normalised Velocity (v/<v>)', saveFilename=outputSaveFileDir+'HistogramsShortlyAfterInfection');
A.plotHistogramsInSameFig(fullTraj_velocityPos00[timesToPlotHist[0]], labelPos00_0, data0_1=fullTraj_velocityPos00[timesToPlotHist[1]], label0_1=labelPos00_1, data0_2=fullTraj_velocityPos00[timesToPlotHist[2]], label0_2=labelPos00_2, data1_0=fullTraj_velocityPos01[timesToPlotHist[0]], label1_0=labelPos01_0, data1_1=fullTraj_velocityPos01[timesToPlotHist[1]], label1_1=labelPos01_1, data1_2=fullTraj_velocityPos01[timesToPlotHist[2]], label1_2=labelPos01_2, data2_0=fullTraj_velocityPos02[timesToPlotHist[0]], label2_0=labelPos02_0, data2_1=fullTraj_velocityPos02[timesToPlotHist[1]], label2_1=labelPos02_1, data2_2=fullTraj_velocityPos02[timesToPlotHist[2]], label2_2=labelPos02_2, xlbl='Normalised Velocity (v/<v>)', saveFilename=outputSaveFileDir+'HistogramsInSameFig');

A.plotHistogramsInSameFig(fullTraj_velocityPos00[timesToPlotHist[0]], labelPos00_0, data0_1=fullTraj_velocityPos00[timesToPlotHist[2]], label0_1=labelPos00_2, data0_2=fullTraj_velocityPos00[timesToPlotHist[3]], label0_2=labelPos00_3, data1_0=fullTraj_velocityPos01[timesToPlotHist[0]], label1_0=labelPos01_0, data1_1=fullTraj_velocityPos01[timesToPlotHist[2]], label1_1=labelPos01_2, data1_2=fullTraj_velocityPos01[timesToPlotHist[3]], label1_2=labelPos01_3, data2_0=fullTraj_velocityPos02[timesToPlotHist[0]], label2_0=labelPos02_0, data2_1=fullTraj_velocityPos02[timesToPlotHist[2]], label2_1=labelPos02_2, data2_2=fullTraj_velocityPos02[timesToPlotHist[3]], label2_2=labelPos02_3, xlbl='Normalised Velocity (v/<v>)', saveFilename=outputSaveFileDir+'HistogramsInSameFig');


# Plot lysis line -- this is a terrible name. This line indicates when the phages were added.
lysisLine_x = np.array([0.5 for i in range(0,9)]);
lysisLine_y = np.array([5.0*i for i in range(0,9)]);
N_lysisLine_y = np.array([200.0*i for i in range(0,9)]);

#Plot skewness of histograms vs time
skewPos00, skewPos01, skewPos02 = A.calculateSkewness(fullTraj_velocityPos00, fullTraj_velocityPos01, fullTraj_velocityPos02);
A.plotDataSetsWithErrorBars(timePos00, skewPos00, 'Control', x1=timePos01, y1=skewPos01, label1='Phage 1', x2=timePos02, y2=skewPos02, label2='Phage 2', xlbl='Time (Minutes)', ylbl='Pearsons Moment Coefficient of Skewness');
outputFile = outputSaveFileDir+'SkewVsTime';
plt.savefig(outputFile);
#title='Skewness of velocity distribution over time'


#Plot number of motile bacteria vs time
A.plotDataSetsWithErrorBars(timePos00, NArrayPos00, 'Control', y0_error=NArrayErrorPos00, x1=timePos01, y1=NArrayPos01, y1_error=NArrayErrorPos01, label1='Phage 1', x2=timePos02, y2=NArrayPos02, y2_error=NArrayErrorPos02, label2='Phage 2', x3=lysisLine_x, y3=N_lysisLine_y, label3='Phage Infection', xlbl='Time (Minutes)', ylbl='Number of Tracked Bacteria');
outputFile = outputSaveFileDir+'NVsTime';
plt.savefig(outputFile);
#title='Tracked Number of Motile Bacteria'

#Plot average velocity vs time with all positions (control and both pahge videos) on one plot.
A.plotDataSetsWithErrorBars(timePos00, velocityPos00, 'Control', y0_error=velocityErrorPos00, x1=timePos01, y1=velocityPos01, y1_error=velocityErrorPos01, label1='Phage 1', x2=timePos02, y2=velocityPos02, y2_error=velocityErrorPos02, label2='Phage 2', x3=lysisLine_x, y3=lysisLine_y, label3='Phage Infection', xlbl='Time (Minutes)', ylbl='Average Velocity (micrometers/second)', plotLegend=False);
plt.legend(loc='lower center')
outputFile = outputSaveFileDir+'AvVelocityVsTime';
plt.savefig(outputFile);
#title='Tracked Average Velocities'
plt.show()

