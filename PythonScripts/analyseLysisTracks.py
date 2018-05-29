#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit
import sys
import os
import fileinput

##########################################################################################################
## Define functions
##########################################################################################################

class analyseTrajectories:
 
    def findLysisEvent(self, A, BIGLIST_traj):
        xPositions = BIGLIST_traj[:, 2]
        yPositions = BIGLIST_traj[:, 3]
        traj_frames = BIGLIST_traj[:, 0];

        #Calculate distance travelled between frames
        displacementArray = A.calcDistBetweenPoints(xPositions, yPositions);
        
        #Calculate velocity at each frame
        velocityArray = displacementArray/(self.NumFramesToAverageOver*self.timePerFrame);
        
        #Calculate normalised direction vector to find correlation between points
        directionCorrelationArray = A.calcNormalisedDirectionCorrelation(displacementArray);

        #Calculate Running average Velocity
        runningVelocityAverage, runningVelocityAverage_error = A.calcRunningAverageVelocity(velocityArray);        

        stopFrame = -1;
        timeStopped = -1.;
        noLysisMarker = -1;
        counter = 0;
        counterJump = 1;
        while (counter < len(runningVelocityAverage)-3*counterJump):
        #for i in range (self.minStopTimeThreshold, len(velocityArray)-self.minStopTimeThreshold):   #Set self.minStopTimeThreshold as minimum amount of time a bacteria needs to have stopped to be considered a stopping event.

            if ((runningVelocityAverage[counter + 2*counterJump] < self.stoppingVelocityThreshold*runningVelocityAverage[counter]) and (runningVelocityAverage[counter + 3*counterJump] <= self.stoppingVelocityThreshold*runningVelocityAverage[counter]) and (np.mean(runningVelocityAverage[(counter+2*counterJump):len(runningVelocityAverage)]) <= self.diffusionThreshold)):
                
                # Detected possible stopping time.
                print 'Detected possible stopping time for ID = '+str(BIGLIST_traj[0, 12]);

                for i in range(counter+2*counterJump, len(runningVelocityAverage)-2):
                    if (runningVelocityAverage[i] >= self.diffusionThreshold and runningVelocityAverage[i+1] >= self.diffusionThreshold and runningVelocityAverage[i+2] >= self.diffusionThreshold):
                        #particle got stuck on glass and started moving again. Not lysis event.
                        noLysisMarker = 1;
                        print 'In first branch. No stop time.'
                        break;

                if (np.mean(runningVelocityAverage) < self.diffusionThreshold):
                    #Particle could possibly be either stuck to the glass or diffusing throughout whole trajectory.
                    swimmerArray = [];
                    for i in range(0, len(runningVelocityAverage)-1):
                        if ((runningVelocityAverage[i-1] > self.diffusionThreshold) and (runningVelocityAverage[i] > self.diffusionThreshold) and (runningVelocityAverage[i+1] > self.diffusionThreshold)):
                            #Found swimmer.
                            swimmerArray.append(i);

                    if (len(swimmerArray) < self.minSwimTimeThreshold):
                        #No swimming for a long period of time. No lysis event.
                        print 'No swimming for long period of time.'
                        noLysisMarker = 1;

                # Continue searching through trajectory
                if(noLysisMarker == 1):
                    noLysisMarker = -1;    #Reset noLysisMarker in case there is lysis event later on in the trajectory.
                    counter = counter + counterJump;
                    print 'noLysisMarker == 1. Continuing through trajectory.'
                    continue;
                
                else:
                    #Detected Lysis event. Return trajectory details.
                    stopFrame = (counter+3)*self.RunningAvTrajLen;
                    stopTime = stopFrame*self.NumFramesToAverageOver+int(traj_frames[0]);
                    timeStopped = float((len(velocityArray) - len(velocityArray[0:stopFrame]))*self.NumFramesToAverageOver);
                    print 'Frame at which stopping event occurs = '+str(stopTime);
                    #A.plotDataSetsWithErrorBars(traj_frames[1:len(traj_frames)], velocityArray, 'Velocity', None, traj_frames[2:len(traj_frames)], runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Average', title='Lysis Event', xlbl='frame', ylbl='Velocity');
                    
                    return stopTime, timeStopped, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray
            
            else:
                counter = counter + counterJump;

        #StopFrame not updated, no lysis event detected.
        return stopFrame, timeStopped, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray
    
    def calcDistBetweenPoints(self, xPositions, yPositions):
        '''
        Calls: nothing

        Returns: an array containing the radial distance a particle has travelled between successive frames separated by self.NumFramesToAverageOver (units = pixels).
        '''

        #check positions arrays are the same length to be sure
        if (len(xPositions) != len(yPositions)):
            print '\nERROR: len(xPositions) != len(yPositions)\n'

        r = np.zeros(int(len(xPositions)/(self.NumFramesToAverageOver)-1));
        counter = 0;
        for i in range(0, len(r)):
            #r[i] = abs(np.sqrt(xPositions[counter+self.NumFramesToAverageOver]**2 + yPositions[counter+self.NumFramesToAverageOver]**2) - np.sqrt(xPositions[counter]**2 + yPositions[counter]**2));
	    r[i] = np.sqrt(abs((xPositions[counter+self.NumFramesToAverageOver] - xPositions[counter])**2 + (yPositions[counter+self.NumFramesToAverageOver] - yPositions[counter])**2)); 
            counter = counter + self.NumFramesToAverageOver;

        return r

    def calcNormalisedDirectionCorrelation(self, displacementArray):
        '''
        n = normalised direction vector. Should be correlated for swimmers and uncorrelated for non-swimmers.
        '''
        n = np.zeros(len(displacementArray)-1);
        for i in range(0,len(n)):
            n[i] = (displacementArray[i+1]-displacementArray[i])/(np.abs(displacementArray[i+1]-displacementArray[i]));
        
        return n

    def calcRunningAverageVelocity(self, velocityArray):
        '''
        Note: first average calculate with only 2 data sets.
        '''
        runningVelocityAverage = np.zeros(int(len(velocityArray)/self.RunningAvTrajLen));
        runningVelocityAverage_error = np.zeros(int(len(velocityArray)/self.RunningAvTrajLen));
        
        counter = self.RunningAvTrajLen;
        for i in range(0, len(runningVelocityAverage)):
            runningVelocityAverage[i] = np.mean(velocityArray[counter-self.RunningAvTrajLen:counter+1]);
            runningVelocityAverage_error[i] = np.std(velocityArray[counter-self.RunningAvTrajLen:counter+1])/np.sqrt(len(velocityArray[counter-self.RunningAvTrajLen:counter+1]));   #error = error on the mean = sigma/sqrt(N)
            counter = counter + self.RunningAvTrajLen;
        
        return runningVelocityAverage, runningVelocityAverage_error

    def plotLysisTrajectory(self, A, lysisVelocityArray, stopTimeArray):
        indexToPlot = 0;
        time = self.timePerFrame*self.NumFramesToAverageOver*np.arange(len(lysisVelocityArray[indexToPlot]));
        
        #Create vertical line to mark where the bacterium has stopped swimming.
        yAxisLine = np.arange(len(np.max(lysisVelocityArray[indexToPlot])));
        xAxisLine = np.array([stopTimeArray[indexToPlot] for i in range(0, len(yAxisLine))]);

        A.plotDataSetsWithErrorBars(time, lysisVelocityArray[indexToPlot], 'Velocity', x1=xAxisLine, y1=yAxisLine, label1='stopTime', title='Lysis Event', xlbl='(time seconds)', ylbl='Velocity (micrometers/second)');
        
        return

    def plotTrajWithSpecificID(self, A, BIGLIST_traj, ID, plotRodLength=0):
        
        # Ensure BIGLIST[ID] exists.
        try:
            print 'Found trajectory beginning at frame '+str(BIGLIST_traj[0, 0]);
            
        # Throw error if BIGLIST_traj does not exist.
        except:
            print 'ERROR: No trajectory identified. Trajectory ID '+str(ID)+' does not exist in frame';
            return

        stopTime_frame, timeStopped, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray = A.findLysisEvent(A, BIGLIST_traj);
        
        #Convert data to units of micrometers and seconds.
        velocityArray = self.pixelsToMicrons*velocityArray;
        runningVelocityAverage = self.pixelsToMicrons*runningVelocityAverage;
        runningVelocityAverage_error = self.pixelsToMicrons*(runningVelocityAverage_error/runningVelocityAverage);
        #traj_frames = self.timePerFrame*traj_frames;

        #Define array to go on x axis of graph
        velocityArray_frames = np.zeros(len(velocityArray));
        counter = self.NumFramesToAverageOver-1;
        for i in range(0, len(velocityArray)):
            velocityArray_frames[i] = traj_frames[counter];
            counter = counter + self.NumFramesToAverageOver

        
        runningVelocityAverage_frames = np.zeros(int(len(velocityArray)/self.RunningAvTrajLen));
        counter = self.RunningAvTrajLen*self.NumFramesToAverageOver-1;
        for i in range(0, len(runningVelocityAverage)):
            runningVelocityAverage_frames[i] = traj_frames[counter];
            counter = counter + self.RunningAvTrajLen*self.NumFramesToAverageOver;

        
        #Covert frames to seconds for plot:
        velocityArray_seconds = self.timePerFrame*(velocityArray_frames - velocityArray_frames[0]);
        runningVelocityAverage_seconds = self.timePerFrame*(runningVelocityAverage_frames - runningVelocityAverage_frames[0]);

        if (plotRodLength == 0):

            #A.plotDataSetsWithErrorBars(traj_frames[1:len(traj_frames)], velocityArray, 'Velocity', np.array(None), traj_frames[2:len(traj_frames)], runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', title='Lysis Event', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Velocity (micrometers/second)');
            #A.plotDataSetsWithErrorBars(velocityArray_frames, velocityArray, 'Velocity', np.array(None), runningVelocityAverage_frames, runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', title='Lysis Event', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Velocity (frames/second)');
            
            A.plotDataSetsWithErrorBars(velocityArray_seconds, velocityArray, 'Velocity', np.array(None), runningVelocityAverage_seconds, runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', xlbl='time (seconds)', ylbl='Velocity (micrometers/second)');

            return stopTime_frame, timeStopped, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray

        else:
            rodLengthArray = self.pixelsToMicrons*BIGLIST[ID][:, 6];
            
            #Fit straight line to these lengths
            M,C = curve_fit(A.fitStraightline, traj_frames[0:len(traj_frames)], rodLengthArray)[0]
            rodLengthFit = M*traj_frames[0:len(traj_frames)] + C;

            #Plot velocity array and running average velocity with error
            #A.plotDataSetsWithErrorBars(velocityArray_frames, velocityArray, 'Velocity', np.array(None), runningVelocityAverage_frames, runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', title='Lysis Event', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Velocity (micrometers/second)');
            
            A.plotDataSetsWithErrorBars(velocityArray_seconds, velocityArray, 'Velocity', np.array(None), runningVelocityAverage_seconds, runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', title='Lysis Event', xlbl='time (seconds)', ylbl='Velocity (micrometers/second)');
            
            #Plot rod length
            #A.plotDataSetsWithErrorBars(traj_frames[0:len(traj_frames)], rodLengthArray, 'Length', np.array(None), traj_frames[0:len(traj_frames)], rodLengthFit, label1='Fit', title='Bacteria Lengths through trajectory', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Rod Length (micrometers)');
            A.plotDataSetsWithErrorBars(velocityArray_frames, rodLengthArray, 'Length', np.array(None), runningVelocityAverage_frames, rodLengthFit, label1='Fit', title='Bacteria Lengths through trajectory', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Rod Length (micrometers)');
                
            return stopTime_frame, timeStopped, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray


    def fitStraightline(self, x, M, C): # this is your 'straight line' y=f(x)
        return M*x + C


    def calcAverageVelocity(self, A, velocityArray):
        '''
        Calls: Nothing

        Returns: Array containing mean velocity for all trajectories and standard error on mean (frames per second)
        '''
        averageVelocityArray = np.zeros(len(velocityArray));
        averageVelocityArray_error = np.zeros(len(velocityArray));
        for i in range(0, len(velocityArray)):
            averageVelocityArray[i] = np.mean(velocityArray[i]);
            averageVelocityArray_error = np.std(velocityArray[i])/np.sqrt(len(velocityArray[i]));

        return averageVelocityArray, averageVelocityArray_error

    def addBIGLISTTrajectories(self, A, BIGLIST_traj1, BIGLIST_traj2):
        NEWLIST = BIGLIST_traj1.tolist();
        APPENDLIST = BIGLIST_traj2.tolist();
        for i in range(0, len(APPENDLIST)):
            NEWLIST.append(APPENDLIST[i]);

        NEWLIST = np.asarray(NEWLIST);
        return NEWLIST

    def __init__(self, minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons, RunningAvTrajLen, minStopTimeThreshold, minSwimTimeThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold):
        self.timePerFrame = timePerFrame;   #For calculating the average velocity over one frame
        self.NumFramesToAverageOver = NumFramesToAverageOver;   #number of frames over which to calculate the displacement over before taking average.
        self.minTrajLen = minTrajLen;
        self.pixelsToMicrons = pixelsToMicrons;
        self.RunningAvTrajLen = RunningAvTrajLen;     #number of frames over which to calculate the running average
        self.minStopTimeThreshold = minStopTimeThreshold;       #minimum length of stopping time.
        self.minSwimTimeThreshold = minSwimTimeThreshold  # min number of frames over which a particle has v > diffusionThreshold to define as swimming. N_frames = minSwimTimeThreshold*NumFramesToAverageOver*RunningAvTrajLen  
        self.initialFrameNum = initialFrameNum;
        self.stoppingVelocityThreshold = stoppingVelocityThreshold;     #fraction drop in average running velocity required to define stopping event.
        self.diffusionThreshold = diffusionThreshold;   #Stopped bacteria still diffuse. Above this threshold, bacteria considered still swimming.

    def appendTrajectory(self, A, BIGLIST_traj, IDs, stopTime_frame, timeStopped, lysisEventDataFileDir, lysisEventDataFile, trajArrayName, AvVelocityBeforeStopArray, AvVelocityBeforeStopArray_error):
        
        # Get the filename as string
        fn = str(input("Append trajectories to Lysis Events File as string ['Y', 'n'] : "))
        #a = sys.argv[1] #imports as string.

        if (fn == 'n'):
            return

        # Function just outputs random errors of [+10, -10] on stop time. Need to alter this by hand later.
        elif (fn == 'Y'):
            A.writeToFile(lysisEventDataFile, [IDs, stopTime_frame, timeStopped, stopTime_frame-10, stopTime_frame+10, AvVelocityBeforeStopArray, AvVelocityBeforeStopArray_error, trajArrayName]);
            np.save(lysisEventDataFileDir+trajArrayName, BIGLIST_traj);
            return

        else:
            print 'Unknown response. Please answer either Y for yes or n for no.';
            return

    # write values to file
    def writeToFile(self, filename, data):
            myFile = open(filename,'a')
            #myFile.write('i Pdf\n\n');
            myFile.write('%s    %d      %f      %f      %f      %f       %f       %s\n' % tuple(data));
            myFile.close()

    # Plot distribution of velocities
    def plotHistogram(self, data, xlbl='bins', ylbl='Frequency'):
	plt.hist(data[:], bins=100)
	plt.xlabel(xlbl);
	plt.ylabel(ylbl);
	#plt.savefig('outputHistogram.pdf')
	#plt.show()
        return

    # Plot up to five sets of data with y error bars on same graph 
    def plotDataSetsWithErrorBars(self, x0, y0, label0, y0_error=np.array(None), x1=np.array(None), y1=np.array(None), y1_error=np.array(None), label1=None, x2=np.array(None), y2=np.array(None), y2_error=np.array(None), label2=None, x3=np.array(None), y3=np.array(None), y3_error=np.array(None), label3=None, x4=np.array(None), y4=np.array(None), y4_error=np.array(None), label4=None, title=None, xlbl=None, ylbl=None):
        
        plt.figure(figsize=(14, 6))
        plt.rc('font', family='serif', size=15);

        print 'y0_error:'
        print y0_error

        if (y0_error.all() == None):
            plt.plot(x0, y0, 'x-', label=label0)
        else:
            plt.errorbar(x0, y0, yerr=y0_error, fmt='o-', label = label0)

        if(y1.all() != None): 
            if (y1_error.all() == None):
                plt.plot(x1, y1, 'x-', label=label1)
            else:
                plt.errorbar(x1, y1, yerr=y1_error, fmt='x-', label = label1)

        if(y2.all() != None): 
            if (y2_error.all() == None):
                plt.plot(x2, y2, '^-', label=label2)
            else:
                plt.errorbar(x2, y2, yerr=y2_error, fmt='^-', label = label2)
        
        if(y3.all() != None): 
            if (y3_error.all() == None):
                plt.plot(x3, y3, 's-', label=label3)
            else:
                plt.errorbar(x3, y3, yerr=y3_error, fmt='s-', label = label3)
        
        if(y4.all() != None): 
            if (y4_error.all() == None):
                plt.plot(x4, y4, 'd-', label=label4)
            else:
                plt.errorbar(x4, y4, yerr=y4_error, fmt='d-', label = label4)
        
        if(title != None): plt.suptitle(title);
        if(xlbl != None): plt.xlabel(xlbl);
        if(ylbl != None): plt.ylabel(ylbl);
        if(y0_error.size != 0 and y1_error.size != 0 and y2_error.size != 0 and y3_error.size != 0 and y4_error.size != 0): plt.legend(loc='upper right');
        
        #plt.rc('font', family='serif', size=15);
        #plt.savefig(outputPlotName)
        #plt.show()

##########################################################################################################
## Main Program begins here
##########################################################################################################

#Declare Variables
NumFramesToAverageOver = 3;
minTrajLen = 30*NumFramesToAverageOver;
RunningAvTrajLen = int(minTrajLen/20);  #Number of points in velocity array to average over in running average.
minSwimTimeThreshold = 2    # min number of frames over which a particle has v > diffusionThreshold to define as swimming. N_frames = minSwimTimeThreshold*NumFramesToAverageOver*RunningAvTrajLen 
fps = 25;
timePerFrame = 1./fps;
pixelsToMicrons = 0.702;    # For x20 Mag
#pixelsToMicrons = 0.354;    #UNKNOWN FOR X40 MAG
minStopTimeThreshold = 1*fps;       #minimum length of time a bacteria is expected to stop before lysing. (frames)
#jumpsToFindLysis = 2*RunningAvTrajLen;  #Jumps in runningVelocityAverage array used to define whether a lysis event is detected.
stoppingVelocityThreshold = 0.6;
D = 0.5;    #Diffusion constant micrometers squared/second
D_ps = D/(pixelsToMicrons**2);    #Diffusion constant micrometers squared/second
#diffusionThreshold = (1/(float(NumFramesToAverageOver)*timePerFrame))*np.sqrt(4*D*(1/pixelsToMicrons)**2*(float(NumFramesToAverageOver)*timePerFrame));     #Above this threshold, bacteria considered to still be swimming.
tau = float(NumFramesToAverageOver*timePerFrame);
diffusionThreshold = (1/tau)*np.sqrt(6*D_ps*tau);   #Above this threshold, bacteria considered to still be swimming.
diffusionThreshold = diffusionThreshold+1.;

#### Declare Filename #######

### MAC ######
#filename = '../../../../../../../Volumes/CBOGGONUSB/Data/20180206-SurfaceVid/trackRods2DtOutput/tracks.dat';

#filename = '../../../../../../../Volumes/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-151805-AsImageSequences/Output-Pos02_Movie0001/filterTracks2DtOutput/tracks_fixed.dat';
#fileDir = '../../../../../../../../Volumes/MyBook/MastersProject/Data/20180227/20180227SurfaceVid1-25fpsx20Mag2000Frames/DDMmovies180227-143931-AsImageSequences/Output-Pos02_Movie0010';
#fileDir = '../../../../../../../../Volumes/MyBook/MastersProject/Data/20180227/20180227SurfaceVid1-25fpsx20Mag2000Frames/DDMmovies180227-143931-AsImageSequences/Output-Pos03_Movie0007';
#fileDir = '../../../../../../../../Volumes/MyBook/MastersProject/Data/20180227/20180227SurfaceVid1-25fpsx20Mag2000Frames/DDMmovies180227-143931-AsImageSequences/Output-Pos01_Movie0009';
#fileDir = '../../../../../../../../Volumes/MyBook/MastersProject/Data/20180227/20180227SurfaceVid2-25fpsx20Mag2000Frames/DDMmovies180227-155825-AsImageSequences/Output-Pos01_Movie0009';
fileDir = '../../../../../../../../Volumes/MyBook/MastersProject/Data/20180227/20180227SurfaceVid2-25fpsx20Mag2000Frames/DDMmovies180227-155825-AsImageSequences/Output-Pos01_Movie0013';
#fileDir = '../../../../../../../../Volumes/MyBook/MastersProject/Data/20180227/20180227SurfaceVid2-25fpsx20Mag2000Frames/DDMmovies180227-155825-AsImageSequences/Output-Pos02_Movie0010';

### UBUNTU ######
#fileDir = '../../../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-152022-AsImageSequences/Output-Pos03_Movie0027';
#fileDir = '../../../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-154132-AsImageSequences/Output-Pos03_Movie0002';

lysisEventDataFileDir = '../../../../../../../../Volumes/CBOGGONUSB/Data/20180227-Surface'+'/lysisTrackingOutput/';

trackFilename = fileDir+'/filterTracks2DtOutput/tracks_fixed.dat';
#lysisEventDataFile = fileDir+'/lysisTrackingOutput/lysisEvents.dat';
lysisEventDataFile = lysisEventDataFileDir+'lysisEvents.dat';


## Create directory if lysisEventDataFileDir does not exist.
try:
    os.stat(lysisEventDataFileDir);
except:
    os.mkdir(lysisEventDataFileDir);


#Read in tracking data
BIGLIST, numberOfFrames = readTrackingFile(trackFilename)
initialFrameNum = int(BIGLIST[0][0, 0]);
A = analyseTrajectories(minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons, RunningAvTrajLen,  minStopTimeThreshold, minSwimTimeThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold);


### Get trajectories, then calculate cumulative velocity of trajectories and apply running average.
# Define threshold for which if average velocity drops significantly. Save these trajectories and cumulative velocities of them. Velocities should fall to zero, within a certain error.
# Calculate time between drop and end of trajectory. If velocity recovers again, don't save the trajectory. Save also the ID and the frame so can check them by eye.
# calc2362ate number of frames/time between stop and end of trajectory and plot stop time.


# Find lysis trajectories
#lysisTraj, lysisVelocityArray, stopTimeArray, lysisDirCorrelationArray, nonLysisTraj, nonLysisVelocityArray = A.getLysisTrajectories(A, BIGLIST);

# Calculate average velocity of non-lysis trajectories
#averageVelocityArray, averageVelocityArray_error = A.calcAverageVelocity(A, nonLysisVelocityArray);

#Convert average velocity to micrometers per second and plot velocity distribution:
#averageVelocityArray_micrometers = averageVelocityArray*pixelsToMicrons;
#A.plotHistogram(averageVelocityArray_micrometers, xlbl='Average Velocity (micrometers)');


#Convert velocity to micrometers per second and plot lysis trajectory:
#lysisVelocityArray = pixelsToMicrons*lysisVelocityArray;
#stopTimeArray = timePerFrame*stopTimeArray;

#Plot lysis trajectory
#A.plotLysisTrajectory(A, lysisVelocityArray, stopTimeArray);


## Plot displacements throughout trajectories
#A.plotRandomTrajectories(A, lysisVelocityArray, 'time (seconds)', 'Velocity (pixels)');


if (len(sys.argv) > 3):
    ## Add two BIGLIST trajectories together:
    print 'sys.argv[1] = '+str(sys.argv[1])
    ID = int(sys.argv[1]);
    ID2 = int(sys.argv[2]);
    ID3 = int(sys.argv[2]);
    #ID = 342;
    #ID2 = 4779;
    #ID = 8720;
    #ID2 = 8721;
    Temp_BIGLIST_traj = A.addBIGLISTTrajectories(A, BIGLIST[ID], BIGLIST[ID2]);
    BIGLIST_traj = A.addBIGLISTTrajectories(A, Temp_BIGLIST_traj, BIGLIST[ID3]);

    stopTime_frame, timeStopped, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray = A.plotTrajWithSpecificID(A, BIGLIST_traj, ID);
    
    # Calculate velociy before stopping time:
    stopIndex = int(stopTime_frame - traj_frames[0] - 10);
    velocityBeforeStopArray = velocityArray[0:stopIndex];
    AvVelocityBeforeStopArray = np.mean(velocityBeforeStopArray);
    AvVelocityBeforeStopArray_error = np.std(velocityBeforeStopArray)/2;

    IDs = str(ID)+'+'+str(ID2)+'+'+str(ID3);
    trajArrayName = fileDir[len(fileDir)-15:len(fileDir)]+'_ID'+str(ID);

    plt.savefig(lysisEventDataFileDir+trajArrayName);
    plt.show()

    A.appendTrajectory(A, BIGLIST_traj, IDs, stopTime_frame, timeStopped, lysisEventDataFileDir, lysisEventDataFile, trajArrayName, AvVelocityBeforeStopArray, AvVelocityBeforeStopArray_error);

elif (len(sys.argv) == 3):
    ## Add two BIGLIST trajectories together:
    print 'sys.argv[1] = '+str(sys.argv[1])
    ID = int(sys.argv[1]);
    ID2 = int(sys.argv[2]);
    #ID = 342;
    #ID2 = 4779;
    #ID = 8720;
    #ID2 = 8721;
    BIGLIST_traj = A.addBIGLISTTrajectories(A, BIGLIST[ID], BIGLIST[ID2]);

    stopTime_frame, timeStopped, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray = A.plotTrajWithSpecificID(A, BIGLIST_traj, ID);
    
    # Calculate velociy before stopping time:
    stopIndex = int(stopTime_frame - traj_frames[0] - 10);
    velocityBeforeStopArray = velocityArray[0:stopIndex];
    AvVelocityBeforeStopArray = np.mean(velocityBeforeStopArray);
    AvVelocityBeforeStopArray_error = np.std(velocityBeforeStopArray)/2;

    IDs = str(ID)+'+'+str(ID2)
    trajArrayName = fileDir[len(fileDir)-15:len(fileDir)]+'_ID'+str(ID);

    plt.savefig(lysisEventDataFileDir+trajArrayName);
    plt.show()

    A.appendTrajectory(A, BIGLIST_traj, IDs, stopTime_frame, timeStopped, lysisEventDataFileDir, lysisEventDataFile, trajArrayName, AvVelocityBeforeStopArray, AvVelocityBeforeStopArray_error);

elif (len(sys.argv) > 1):
    ID = int(sys.argv[1]);
    #ID = 2362;
    BIGLIST_traj = BIGLIST[ID];
    stopTime_frame, timeStopped, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray = A.plotTrajWithSpecificID(A, BIGLIST_traj, ID);
    ##A.plotTrajWithSpecificID(A, BIGLIST, ID, plotRodLength=1);
    
    # Calculate velociy before stopping time:
    stopIndex = int(stopTime_frame - traj_frames[0] - 10);
    velocityBeforeStopArray = velocityArray[0:stopIndex];
    AvVelocityBeforeStopArray = np.mean(velocityBeforeStopArray);
    AvVelocityBeforeStopArray_error = np.std(velocityBeforeStopArray)/2;

    IDs = ID;
    trajArrayName = fileDir[len(fileDir)-15:len(fileDir)]+'_ID'+str(ID);
    
    plt.savefig(lysisEventDataFileDir+trajArrayName);
    plt.show()

    A.appendTrajectory(A, BIGLIST_traj, IDs, stopTime_frame, timeStopped, lysisEventDataFileDir, lysisEventDataFile, trajArrayName, AvVelocityBeforeStopArray, AvVelocityBeforeStopArray_error);


else:
    print 'unrecognised input arguments'




#np.load(trajArrayName)

'''
Notes on things not done:

    defines lysis event using running avergage velocity. So does not calculate an exact time for lysis. --- this may not actually be true anymore

    lysis event definition includes minimum number of frames over which v > diffusionThreshold. But no contingency that number of detected swimming events in trajectory are all one after another. So a swimming type that moves very close to diffusion limit could come up as lysis event.
    
'''
