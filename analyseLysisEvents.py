#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit

##########################################################################################################
## Define functions
##########################################################################################################

class analyseTrajectories:

    def getLysisTrajectories(self, A, BIGLIST):
        #Get xPositions and yPositions from BIGLIST
        lysisTraj = [];
        lysisVelocityArray = [];
        nonLysisTraj = [];
        nonLysisVelocityArray = [];

        #AvVelocityArray = np.zeros(len(BIGLIST));
        #displacementArrayByFrame = [];
        #velocityArrayByFrame = [];   #Array of velocities averaged over NumFramesToAverageOver
        i = 0;
        for i in range(0, len(BIGLIST)):   # iterate through trajectories
            if(len(BIGLIST[i]) < self.minTrajLen):
                continue;
            else:
                stopTime_frame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray = A.findLysisEvent(A, BIGLIST[i]);
                #Might need to remove zeroes from end of list.
                
                if (stopTime_frame != -1):
                    lysisTraj.append(i);
                    lysisVelocityArray.append(velocityArray);

                else:
                    nonLysisTraj.append(i);
                    nonLysisVelocityArray.append(velocityArray);
        
        lysisTraj = np.asarray(lysisTraj);
        lysisVelocityArray = np.asarray(lysisVelocityArray);
        nonLysisVelocityArray = np.asarray(nonLysisVelocityArray);
        nonLysisVelocityArray = np.asarray(nonLysisVelocityArray);

        return lysisTraj, lysisVelocityArray, nonLysisTraj, nonLysisVelocityArray 

    
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
        
####  OLD Version
#        stopFrame = -1;
#        for i in range (self.minStopTimeThreshold, len(velocityArray)-self.minStopTimeThreshold):   #Set self.minStopTimeThreshold as minimum amount of time a bacteria needs to have stopped to be considered a stopping event.
#            if ((velocityArray[i] < self.stoppingVelocityThreshold*runningVelocityAverage[i-1]) and (velocityArray[i+1] < 1.5*self.stoppingVelocityThreshold*runningVelocityAverage[i-1]) and (np.mean(runningVelocityAverage[i:len(runningVelocityAverage)]) <= self.diffusionThreshold)):
#                # Detected lysis event.
#                stopFrame = i;
#                timeStopped = len(velocityArray) - len(velocityArray[0:stopFrame]);
#                #A.plotDataSetsWithErrorBars(traj_frames[1:len(traj_frames)], velocityArray, 'Velocity', None, traj_frames[2:len(traj_frames)], runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Average', title='Lysis Event', xlbl='frame', ylbl='Velocity');
#            
#                return stopFrame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray
#
#
#        #StopFrame not updated, no lysis event detected.
#        return stopFrame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray

        stopFrame = -1;
        noLysisMarker = -1;
        counter = 0;
        counterJump = 1;
        while (counter < len(runningVelocityAverage)-2*counterJump):
        #for i in range (self.minStopTimeThreshold, len(velocityArray)-self.minStopTimeThreshold):   #Set self.minStopTimeThreshold as minimum amount of time a bacteria needs to have stopped to be considered a stopping event.

            if ((runningVelocityAverage[counter + counterJump] < self.stoppingVelocityThreshold*runningVelocityAverage[counter]) and (runningVelocityAverage[counter + 2*counterJump] <= self.stoppingVelocityThreshold*runningVelocityAverage[counter]) and (np.mean(runningVelocityAverage[counter:len(runningVelocityAverage)]) <= self.diffusionThreshold)):
                
                # Detected possible lysis event.
                
                for i in range(counter, len(runningVelocityAverage)-1):
                    if (runningVelocityAverage[i] >= self.diffusionThreshold and runningVelocityAverage[i+1] >= self.diffusionThreshold):
                        #particle got stuck on glass and started moving again. Not lysis event.
                        noLysisMarker = 1;
                        break;

                if (np.mean(runningVelocityAverage) < self.diffusionThreshold):
                    #Particle could possibly be stuck to the glass or diffusing throughout whole trajectory.
                    swimmerArray = [];
                    for i in range(0, len(runningVelocityAverage)):
                        if ((runningVelocityAverage[i-1] > self.diffusionThreshold) and (runningVelocityAverage[i] > self.diffusionThreshold) and (runningVelocityAverage[i+1] > self.diffusionThreshold)):
                            #Found swimmer.
                            swimmerArray.append(i);

                    if (len(swimmerArray) < self.minSwimTimeThreshold):
                        #No swimming for a long period of time. No lysis event.
                        noLysisMarker = 1;

                # Continue searching through trajectory
                if(noLysisMarker == 1):
                    noLysisMarker = -1;    #Reset noLysisMarker in case there is lysis event later on in the trajectory.
                    counter = counter + counterJump;
                    continue;
                
                else:
                    #Detected Lysis event. Return trajectory details.
                    stopFrame = counter*self.RunningAvTrajLen;
                    timeStopped = len(velocityArray) - len(velocityArray[0:stopFrame]);
                    #A.plotDataSetsWithErrorBars(traj_frames[1:len(traj_frames)], velocityArray, 'Velocity', None, traj_frames[2:len(traj_frames)], runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Average', title='Lysis Event', xlbl='frame', ylbl='Velocity');
            
                    return stopFrame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray
            
            else:
                counter = counter + counterJump;

        #StopFrame not updated, no lysis event detected.
        return stopFrame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray
    
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

    
    def plotTrajWithSpecificID(self, A, BIGLIST, ID, plotRodLength=0):
        
        # Ensure BIGLIST[ID] exists.
        try:
            print 'Found trajectory beginning at frame '+str(BIGLIST[ID][0, 0]);
            
        # Throw error if BIGLIST[ID] does not exist.
        except:
            print 'ERROR: No trajectory identified. Trajectory ID '+str(ID)+' does not exist in frame '+str(frame);
            return

        stopTime_frame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames, directionCorrelationArray = A.findLysisEvent(A, BIGLIST[ID]);
        
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


        if (plotRodLength == 0):

            #A.plotDataSetsWithErrorBars(traj_frames[1:len(traj_frames)], velocityArray, 'Velocity', np.array(None), traj_frames[2:len(traj_frames)], runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', title='Lysis Event', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Velocity (micrometers/second)');
            A.plotDataSetsWithErrorBars(velocityArray_frames, velocityArray, 'Velocity', np.array(None), runningVelocityAverage_frames, runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', title='Lysis Event', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Velocity (micrometers/second)');

            return

        else:
            rodLengthArray = self.pixelsToMicrons*BIGLIST[ID][:, 6];
            
            #Fit straight line to these lengths
            M,C = curve_fit(A.fitStraightline, traj_frames[0:len(traj_frames)], rodLengthArray)[0]
            rodLengthFit = M*traj_frames[0:len(traj_frames)] + C;

            #Plot velocity array and running average velocity with error
            A.plotDataSetsWithErrorBars(velocityArray_frames, velocityArray, 'Velocity', np.array(None), runningVelocityAverage_frames, runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', title='Lysis Event', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Velocity (micrometers/second)');
            
            #Plot rod length
            #A.plotDataSetsWithErrorBars(traj_frames[0:len(traj_frames)], rodLengthArray, 'Length', np.array(None), traj_frames[0:len(traj_frames)], rodLengthFit, label1='Fit', title='Bacteria Lengths through trajectory', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Rod Length (micrometers)');
            A.plotDataSetsWithErrorBars(velocityArray_frames, rodLengthArray, 'Length', np.array(None), runningVelocityAverage_frames, rodLengthFit, label1='Fit', title='Bacteria Lengths through trajectory', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Rod Length (micrometers)');
                
            return


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


#    def calcMSD(self, F, y, x, n, N):
#	r = F.calcDisplacement(y, x);
#	print('N = '+str(N));
#	print('len(r) = '+str(len(r)));
#	MSDCumulative = np.zeros((N-1,len(x[0])));
#	MSD = np.zeros(len(r[0]));
#	for i in range(0,len(x[0])):
#	    for j in range(0,N-1):
#		MSD[i] = MSD[i] + (n/N)*(r[j+1][i] - r[j][i])**2;
#		MSDCumulative[j][i] = MSD[i];
#	return MSD, MSDCumulative;
#
#    # Calculate the displacement of a particle in microns given cartesian data points.
#    def calcDisplacement(self, y, x):
#	r = np.zeros((len(x), len(x[0])));  #r = [no. points in trajectory][particle number]
#	for j in range(0,len(x)):
#	    for i in range(0,len(x[0])):
#		r[j][i] = np.sqrt(x[j][i]**2 + y[j][i]**2);
#	return r
#
#    # To find the average mean squared displacement
#    def averageMSD(self, MSDData):
#	meanMSD = np.zeros(len(MSDData));
#	for j in range(0, len(MSDData)):
#	    print(j)
#	    accumulator = 0;
#	    for i in range(0, len(MSDData[0])-1):
#		accumulator = accumulator + MSDData[j][i];
#		#array[j] = np.nansum([array[j], im[im.columns[j]][i]])
#		#array[j] = array[j]+im[im.columns[j]][i];
#		#print(array)
#	    meanMSD[j] = accumulator/len(MSDData[0]);    #Calculate average of array
#	return meanMSD
#
#    def checkTrajUniformity(self, A, BIGLIST):
#        '''
#        To check trajectories are roughly uniform so that it is reasonable to calculate average velocity from average displacements over entire length of trajectory.
#        Is there is a significant gradient in the data, this would not be reasonable.
#        This might happen if the particles are slowing down fast enough during the trajectory.
#
#        Calls: 
#            - plotDataSets
#
#        Returns: Nothing
#        '''
#        trajList = [];
#        i = 0;
#        while(i < 5):
#            temp = BIGLIST[int(len(BIGLIST)*np.random.uniform())];
#            if(len(temp) < self.minTrajLen):
#                continue;
#            else:
#                trajList.append(temp);
#                i = i+1;
#
##        traj0 = np.asarray(trajList[0]);
##        traj1 = np.asarray(trajList[1]);
##        traj2 = np.asarray(trajList[2]);
##        traj3 = np.asarray(trajList[3]);
##        traj4 = np.asarray(trajList[4]);
#        
#        AvVelocityArray_traj0, displacementArray_SeveralFrames_traj0 = A.calcAverageParticleVelocity(A, np.asarray(trajList[0]));
#        AvVelocityArray_traj1, displacementArray_SeveralFrames_traj1 = A.calcAverageParticleVelocity(A, np.asarray(trajList[1]));
#        AvVelocityArray_traj2, displacementArray_SeveralFrames_traj2 = A.calcAverageParticleVelocity(A, np.asarray(trajList[2]));
#        AvVelocityArray_traj3, displacementArray_SeveralFrames_traj3 = A.calcAverageParticleVelocity(A, np.asarray(trajList[3]));
#        AvVelocityArray_traj4, displacementArray_SeveralFrames_traj4 = A.calcAverageParticleVelocity(A, np.asarray(trajList[4]));
#
#        A.plotDataSets(np.arange(len(displacementArray_SeveralFrames_traj0)), displacementArray_SeveralFrames_traj0, 'traj0', x1=np.arange(len(displacementArray_SeveralFrames_traj1)), y1=displacementArray_SeveralFrames_traj1, label1='traj1', x2=np.arange(len(displacementArray_SeveralFrames_traj2)), y2=displacementArray_SeveralFrames_traj2, label2='traj2', x3=np.arange(len(displacementArray_SeveralFrames_traj3)), y3=displacementArray_SeveralFrames_traj3, label3='traj3', x4=np.arange(len(displacementArray_SeveralFrames_traj4)), y4=displacementArray_SeveralFrames_traj4, label4='traj4', xlbl='Iteration (iterates over '+str(self.NumFramesToAverageOver)+' Frame(s))', ylbl='Displacement Over '+str(self.NumFramesToAverageOver)+' Frame (pixels)');
#
#        return
#

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

#    def calcAverageVelocitiesForAllTraj(self, A, BIGLIST):        
#        '''
#        Note: This assumes the particles behave uniformly throughout the course of the trajectory i.e they do not slow down.
#        
#        Calls:
#            - calcAverageParticleVelocity
#        
#        Returns:
#            - average velocity of all particle trajectories in pixels per seconds. 
#
#        '''
#        #Get xPositions and yPositions from BIGLIST
#        AvVelocityArray = np.zeros(len(BIGLIST));
#        displacementArrayByFrame = [];
#        velocityArrayByFrame = [];   #Array of velocities averaged over NumFramesToAverageOver
#        i = 0;
#        for i in range(0, len(BIGLIST)):   # iterate through trajectories
#            if(len(BIGLIST[i]) < self.minTrajLen):
#                continue;
#            else:
#                AvVelocityArray[i], temp_velocityArrayByFrame, displacementArray_SeveralFrames = A.calcAverageParticleVelocity(A, BIGLIST[i]);
#                velocityArrayByFrame.append(temp_velocityArrayByFrame[np.nonzero(temp_velocityArrayByFrame)]);
#                displacementArrayByFrame.append(displacementArray_SeveralFrames/self.NumFramesToAverageOver);
#    
#        AvVelocityArray = AvVelocityArray[np.nonzero(AvVelocityArray)];  # return non zero elements of array (elements with long enough trajectories).
#        displacementArrayByFrame = np.asarray(displacementArrayByFrame);
#        velocityArrayByFrame = np.asarray(velocityArrayByFrame);
#
#        return AvVelocityArray, velocityArrayByFrame, displacementArrayByFrame
#
#    def calcAverageParticleVelocity(self, A, BIGLIST_traj):
#        '''
#        Calls:
#            - calcAverageDisplacement
#
#        Returns:
#            - average velocity of a single particle (single trajectory) in pixels per seconds
#            - array containing total displacement over 'NumFramesToAverageOver' frames.
#        '''
#        xPositions = np.zeros(len(BIGLIST_traj));
#        yPositions = np.zeros(len(BIGLIST_traj));
#        for j in range(0, len(xPositions)):     # iterate through all tracked frames in this trajectory
#            xPositions[j] = BIGLIST_traj[j][2];
#            yPositions[j] = BIGLIST_traj[j][3];
#
#        #Calculate distance between points
#        displacementArray = A.calcDistBetweenPoints(xPositions, yPositions);
#        AvDisplacement, displacementArray_SeveralFrames = A.calcAverageDisplacement(displacementArray);
#        
#        AvVelocity = AvDisplacement/self.timePerFrame;
#        velocityArrayByFrame = displacementArray_SeveralFrames/(self.NumFramesToAverageOver*self.timePerFrame);
#
#        return AvVelocity, velocityArrayByFrame, displacementArray_SeveralFrames
#
#    
#
#    def calcAverageDisplacement(self, displacementArray):
#        '''
#        Calculate the the distance travelled over 'NumFramesToAverageOver' frames. 
#        Then average all these displacements to find the average displacement of particle per frame in video.
#
#        Calls: nothing.
#        
#        Returns: 
#            - float of average displacement per frame of a particle in pixels.
#            - array containing total displacement over 'NumFramesToAverageOver' frames.
#        '''
#        displacementArray_SeveralFrames = np.zeros(int(len(displacementArray)/self.NumFramesToAverageOver));
#        tempArray = np.zeros(self.NumFramesToAverageOver);
#        counter = 0;
#        for i in range(0, len(displacementArray_SeveralFrames)):
#            for j in range(0, self.NumFramesToAverageOver):
#                tempArray[j] = displacementArray[counter+j];
#
#            displacementArray_SeveralFrames[i] = np.sum(tempArray);
#            counter = counter + self.NumFramesToAverageOver;
#        
#        #Calculate average displacement per frame
#        AvDisplacement = np.sum(displacementArray_SeveralFrames)/(self.NumFramesToAverageOver*len(displacementArray_SeveralFrames));
#        
#        #Set displacementArray_severalFrames to global
#        #self.displacementArrayByFrame = displacementArray_severalFrames/self.NumFramesToAverageOver;
#        #self.velocityArrayByFrame = displacementArray_severalFrames/(self.NumFramesToAverageOver*self.timePerFrame);    #Array of velocities averaged over NumFramesToAverageOver 
#
#        return AvDisplacement, displacementArray_SeveralFrames

    
    #def plotRandomTrajectories(self, A, traj, xlbl='Iteration (iterates over '+str(self.NumFramesToAverageOver)+' Frame(s))', ylbl='Displacement Over '+str(self.NumFramesToAverageOver)+' Frame (pixels)'):
    def plotRandomTrajectories(self, A, traj, xlbl, ylbl, traj2=None, xlbl2=None, ylbl2=None):
        '''
        Plots the displacement of the particle for each iteration of NumFramesToAverageOver for 5 random trajectories.

        Inputs: 
            - traj = 2D numpy array of all trajectories
            - xlbl, ylbl = strings for plot

        Calls: 
            - plotDataSets

        Returns: Nothing
        '''
        velTrajList = [];
        if(xlbl2 != None): 
            dispTrajList = [];
        
        i = 0;
        for i in range(0,5):
            rand = np.random.uniform();
            #temp = traj[int(len(traj)*rand)];
            velTrajList.append(traj[int(len(traj)*rand)]);
            if(xlbl2 != None): 
                dispTrajList.append(traj2[int(len(traj2)*rand)]);
            i = i+1;

        time_traj0 = np.arange(len(velTrajList[0]))*self.timePerFrame*self.NumFramesToAverageOver;
        time_traj1 = np.arange(len(velTrajList[1]))*self.timePerFrame*self.NumFramesToAverageOver;
        time_traj2 = np.arange(len(velTrajList[2]))*self.timePerFrame*self.NumFramesToAverageOver;
        time_traj3 = np.arange(len(velTrajList[3]))*self.timePerFrame*self.NumFramesToAverageOver;
        time_traj4 = np.arange(len(velTrajList[4]))*self.timePerFrame*self.NumFramesToAverageOver;
        
        print 'xlbl2 = '+str(xlbl2);
        #print 'len(traj0) = '+str(len(trajList[0]));
        #print 'len(traj1) = '+str(len(trajList[1]));
        #print 'len(traj2) = '+str(len(trajList[2]));
        #print 'len(traj3) = '+str(len(trajList[3]));
        #print 'len(traj4) = '+str(len(trajList[4]));

#        traj0 = np.asarray(trajList[0]);
#        traj1 = np.asarray(trajList[1]);
#        traj2 = np.asarray(trajList[2]);
#        traj3 = np.asarray(trajList[3]);
#        traj4 = np.asarray(trajList[4]);
        
#        AvVelocityArray_traj0, displacementArray_SeveralFrames_traj0 = A.calcAverageParticleVelocity(A, np.asarray(trajList[0]));
#        AvVelocityArray_traj1, displacementArray_SeveralFrames_traj1 = A.calcAverageParticleVelocity(A, np.asarray(trajList[1]));
#        AvVelocityArray_traj2, displacementArray_SeveralFrames_traj2 = A.calcAverageParticleVelocity(A, np.asarray(trajList[2]));
#        AvVelocityArray_traj3, displacementArray_SeveralFrames_traj3 = A.calcAverageParticleVelocity(A, np.asarray(trajList[3]));
#        AvVelocityArray_traj4, displacementArray_SeveralFrames_traj4 = A.calcAverageParticleVelocity(A, np.asarray(trajList[4]));


        A.plotDataSets(time_traj0, velTrajList[0], 'traj0', x1=time_traj1, y1=velTrajList[1], label1='traj1', x2=time_traj2, y2=velTrajList[2], label2='traj2', x3=time_traj3, y3=velTrajList[3], label3='traj3', x4=time_traj4, y4=velTrajList[4], label4='traj4', xlbl=xlbl, ylbl=ylbl);
         
        if(xlbl2 != None): 
            A.plotDataSets(time_traj0, dispTrajList[0], 'traj0', x1=time_traj1, y1=dispTrajList[1], label1='traj1', x2=time_traj2, y2=dispTrajList[2], label2='traj2', x3=time_traj3, y3=dispTrajList[3], label3='traj3', x4=time_traj4, y4=dispTrajList[4], label4='traj4', xlbl=xlbl2, ylbl=ylbl2);

        return


    # Plot distribution of velocities
    def plotHistogram(self, data, xlbl='bins', ylbl='Frequency'):
	plt.hist(data[:], bins=100)
	plt.xlabel(xlbl);
	plt.ylabel(ylbl);
	#plt.savefig('outputHistogram.pdf')
	#plt.show()
        return

    # Plot up to five sets of data on same graph 
    def plotDataSets(self, x0, y0, label0, x1=None, y1=None, label1=None, x2=None, y2=None, label2=None, x3=None, y3=None, label3=None, x4=None, y4=None, label4=None, title=None, xlbl=None, ylbl=None):
        plt.figure()
        plt.plot(x0, y0, 'o-', label=label0)
        
        if(y1.any() != None): 
            plt.plot(x1, y1, 'x-', label=label1)
        
        if(y2.any() != None): 
            plt.plot(x2, y2, '^-', label=label2)
        
        if(y3.any() != None): 
            plt.plot(x3, y3, 's-', label=label3)
        
        if(y4.any() != None): 
            plt.plot(x4, y4, 'd-', label=label4)
        
        if(title != None): plt.suptitle(title);
        if(xlbl != None): plt.xlabel(xlbl);
        if(ylbl != None): plt.ylabel(ylbl)
        plt.legend(loc='upper right');
        #plt.rc('font', family='serif', size=15);
        #plt.savefig(outputPlotName)
        #plt.show()


    # Plot up to five sets of data with y error bars on same graph 
    def plotDataSetsWithErrorBars(self, x0, y0, label0, y0_error=np.array(None), x1=np.array(None), y1=np.array(None), y1_error=np.array(None), label1=None, x2=np.array(None), y2=np.array(None), y2_error=np.array(None), label2=None, x3=np.array(None), y3=np.array(None), y3_error=np.array(None), label3=None, x4=np.array(None), y4=np.array(None), y4_error=np.array(None), label4=None, title=None, xlbl=None, ylbl=None):
        
        plt.figure()
       
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
        if(ylbl != None): plt.ylabel(ylbl)
        plt.legend(loc='upper right');
        #plt.rc('font', family='serif', size=15);
        #plt.savefig(outputPlotName)
        #plt.show()

##########################################################################################################
## Main Program begins here
##########################################################################################################

#Declare Variables
NumFramesToAverageOver = 3;
minTrajLen = 30*NumFramesToAverageOver;
RunningAvTrajLen = int(minTrajLen/10);  #Number of points in velocity array to average over in running average.
minSwimTimeThreshold = 6    # min number of frames over which a particle has v > diffusionThreshold to define as swimming. N_frames = minSwimTimeThreshold*NumFramesToAverageOver*RunningAvTrajLen 
fps = 50;
timePerFrame = 1./fps;
pixelsToMicrons = 0.702;    # For x20 Mag
#pixelsToMicrons = 0.354;    #UNKNOWN FOR X40 MAG
minStopTimeThreshold = 1*fps;       #minimum length of time a bacteria is expected to stop before lysing. (frames)
#jumpsToFindLysis = 2*RunningAvTrajLen;  #Jumps in runningVelocityAverage array used to define whether a lysis event is detected.
stoppingVelocityThreshold = 0.2;
D = 0.34;    #Diffusion constant micrometers/second
diffusionThreshold = (1/(float(NumFramesToAverageOver)*timePerFrame))*np.sqrt(4*D*(1/pixelsToMicrons)**2*(float(NumFramesToAverageOver)*timePerFrame));     #Above this threshold, bacteria considered to still be swimming.
#B = analyseLysisEvents(A)


#filename = 'testcaseTracks.dat';
#filename = 'OutputPos00_Movie0000-BrightnessContrast/TrackRodsOutput/tracks.dat';
#filename = '../Data/171201-DDM/Output-Pos01_Movie0000/trackRods2DtOutput/tracks.dat';
filename = '../../../../../../../Volumes/CBOGGONUSB/Data/20180206-SurfaceVid/trackRods2DtOutput/tracks.dat';


#Read in tracking data
BIGLIST, numberOfFrames = readTrackingFile(filename)
initialFrameNum = int(BIGLIST[0][0, 0]);
A = analyseTrajectories(minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons, RunningAvTrajLen,  minStopTimeThreshold, minSwimTimeThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold);

#BIGLIST, PROPERTIES = readCoordinateFile(filename)
#print BIGLIST
#print numberOfFrames


#### Calculate velocity of particles
#AvVelocityArray, velocityArrayByFrame, displacementArrayByFrame = A.calcAverageVelocitiesForAllTraj(A, BIGLIST);

## Plot distribution of velocities
##A.plotHistogram(AvVelocityArray, xlbl='Average Velocity (pixels)');
#
## Check particles behave uniformly throughout trajectories
##A.checkTrajUniformity(A, BIGLIST); #Not relevant anymore!
#
## Calculate average velocity in micrometers/second
#AvVelocityArray_micrometers = pixelsToMicrons*AvVelocityArray;
#velocityArrayByFrame_micrometers = pixelsToMicrons*velocityArrayByFrame;
#displacementArrayByFrame_micrometers = pixelsToMicrons*displacementArrayByFrame;
#
## Plot distribution of velocities
#A.plotHistogram(AvVelocityArray_micrometers, xlbl='Average Velocity (micrometers)');
#
## Plot velocities and displacements throughout trajectories
#A.plotRandomTrajectories(A, velocityArrayByFrame_micrometers, 'time (seconds)', 'Velocity (micrometers/second)', traj2=displacementArrayByFrame_micrometers, xlbl2='time (seconds)', ylbl2='Displacement (micrometers)');
#
## Plot displacements throughout trajectories
##A.plotRandomTrajectories(A, displacementArrayByFrame_micrometers, 'time (seconds)', 'Displacement Over '+str(NumFramesToAverageOver)+' Frame (micrometers)');




### Get trajectories, then calculate cumulative velocity of trajectories and apply running average.
# Define threshold for which if average velocity drops significantly. Save these trajectories and cumulative velocities of them. Velocities should fall to zero, within a certain error.
# Calculate time between drop and end of trajectory. If velocity recovers again, don't save the trajectory. Save also the ID and the frame so can check them by eye.
# calculate number of frames/time between stop and end of trajectory and plot stop time.


# Find lysis trajectories
lysisTraj, lysisVelocityArray, nonLysisTraj, nonLysisVelocityArray = A.getLysisTrajectories(A, BIGLIST);

# Calculate average velocity of non-lysis trajectories
averageVelocityArray, averageVelocityArray_error = A.calcAverageVelocity(A, nonLysisVelocityArray);


#Convert average velocity to micrometers per second and plot velocity distribution:
averageVelocityArray_micrometers = averageVelocityArray*pixelsToMicrons;
A.plotHistogram(averageVelocityArray_micrometers, xlbl='Average Velocity (micrometers)');

## Plot displacements throughout trajectories
#A.plotRandomTrajectories(A, lysisVelocityArray, 'time (seconds)', 'Velocity (pixels)');


ID = 14;
A.plotTrajWithSpecificID(A, BIGLIST, ID);
#A.plotTrajWithSpecificID(A, BIGLIST, ID, plotRodLength=1);

plt.show()


'''
Notes on things not done:

    defines lysis event using running avergage velocity. So does not calculate an exact time for lysis. --- this may not actually be true anymore

    lysis event definition includes minimum number of frames over which v > diffusionThreshold. But no contingency that number of detected swimming events in trajectory are all one after another. So a swimming type that moves very close to diffusion limit could come up as lysis event.
    
'''
