#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit
import sys

##########################################################################################################
## Define functions
##########################################################################################################

class analyseTrajectories:

#    def getLysisTrajectories(self, A, BIGLIST):
#        #Get xPositions and yPositions from BIGLIST
#        lysisTraj = [];
#        lysisVelocityArray = [];
#        nonLysisTraj = [];
#        nonLysisVelocityArray = [];
#
#        #AvVelocityArray = np.zeros(len(BIGLIST));
#        #displacementArrayByFrame = [];
#        #velocityArrayByFrame = [];   #Array of velocities averaged over NumFramesToAverageOver
#        i = 0;
#        for i in range(0, len(BIGLIST)):   # iterate through trajectories
#            if(len(BIGLIST[i]) < self.minTrajLen):
#                continue;
#            else:
#                stopTime_frame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames = A.findLysisEvent(A, BIGLIST[i]);
#                #Might need to remove zeroes from end of list.
#                
#                if (stopTime_frame != -1):
#                    lysisTraj.append(i);
#                    lysisVelocityArray.append(velocityArray);
#
#                else:
#                    nonLysisTraj.append(i);
#                    nonLysisVelocityArray.append(velocityArray);
#        
#        lysisTraj = np.asarray(lysisTraj);
#        lysisVelocityArray = np.asarray(lysisVelocityArray);
#        nonLysisVelocityArray = np.asarray(nonLysisVelocityArray);
#        nonLysisVelocityArray = np.asarray(nonLysisVelocityArray);
#
#        return lysisTraj, lysisVelocityArray, nonLysisTraj, nonLysisVelocityArray 
#
#    
#    def findLysisEvent(self, A, BIGLIST_traj):
#        xPositions = BIGLIST_traj[:, 2]
#        yPositions = BIGLIST_traj[:, 3]
#        traj_frames = BIGLIST_traj[:, 0];
#
#        #Calculate distance travelled between frames
#        displacementArray = A.calcDistBetweenPoints(xPositions, yPositions);
#        
#        #Calculate velocity at each frame
#        velocityArray = displacementArray/self.timePerFrame;
#
#        #Calculate Running average Velocity
#        runningVelocityAverage, runningVelocityAverage_error, runningVelocityBackwardsAverage, runningVelocityBackwardsAverage_error = A.calcRunningAverageVelocity(velocityArray);
#        
#        stopFrame = -1;
#        for i in range (self.minStopTimeThreshold, len(velocityArray)-self.minStopTimeThreshold):   #Set self.minStopTimeThreshold as minimum amount of time a bacteria needs to have stopped to be considered a stopping event.
#            if ((velocityArray[i] < self.stoppingVelocityThreshold*runningVelocityAverage[i-1]) and (velocityArray[i+1] < 1.5*self.stoppingVelocityThreshold*runningVelocityAverage[i-1]) and (np.mean(runningVelocityAverage[i:len(runningVelocityAverage)]) <= self.diffusionThreshold)):
#                # Detected lysis event.
#                stopFrame = i;
#                timeStopped = len(velocityArray) - len(velocityArray[0:stopFrame]);
#                print 'timeStopped = '+str(timeStopped)
#                print 'stopFrame = '+str(stopFrame)
#                #A.plotDataSetsWithErrorBars(traj_frames[1:len(traj_frames)], velocityArray, 'Velocity', None, traj_frames[2:len(traj_frames)], runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Average', title='Lysis Event', xlbl='frame', ylbl='Velocity');
#            
#                return stopFrame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames
#
#
#        #StopFrame not updated, no lysis event detected.
#        return stopFrame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames
#
#    
#    def calcRunningAverageVelocity(self, velocityArray):
#        '''
#        Note: first average calculate with only 2 data sets.
#        '''
#        runningVelocityAverage = np.zeros(len(velocityArray)-1);
#        runningVelocityAverage_error = np.zeros(len(velocityArray)-1);
#        runningVelocityBackwardsAverage = np.zeros(len(velocityArray)-1);
#        runningVelocityBackwardsAverage_error = np.zeros(len(velocityArray)-1);
#        for i in range(0,len(runningVelocityAverage)):
#            runningVelocityAverage[i] = np.mean(velocityArray[0:i+1]);
#            runningVelocityAverage_error[i] = np.std(velocityArray[0:i+1])/np.sqrt(len(velocityArray[0:i+1]))   #error = error on the mean = sigma/sqrt(N)
#            
#            runningVelocityBackwardsAverage[i] = np.mean(velocityArray[len(runningVelocityAverage)-i-1:len(runningVelocityAverage)]);
#            runningVelocityBackwardsAverage_error[i] = np.std(velocityArray[len(runningVelocityAverage)-i-1:len(runningVelocityAverage)])/np.sqrt(len(velocityArray[len(runningVelocityAverage)-i-1:len(runningVelocityAverage)]))   #error = error on the mean = sigma/sqrt(N)
#
#        return runningVelocityAverage, runningVelocityAverage_error, runningVelocityBackwardsAverage, runningVelocityBackwardsAverage_error
#
#
#    
#    def plotTrajWithSpecificID(self, A, BIGLIST, ID, plotRodLength=0):
#        
#        # Ensure BIGLIST[ID] exists.
#        try:
#            print 'Found trajectory beginning at frame '+str(BIGLIST[ID][0, 0]);
#            
#        # Throw error if BIGLIST[ID] does not exist.
#        except:
#            print 'ERROR: No trajectory identified. Trajectory ID '+str(ID)+' does not exist in frame '+str(frame);
#            return
#
#        stopTime_frame, velocityArray, runningVelocityAverage, runningVelocityAverage_error, traj_frames = A.findLysisEvent(A, BIGLIST[ID]);
#        
#        #Convert data to units of micrometers and seconds.
#        velocityArray = self.pixelsToMicrons*velocityArray;
#        runningVelocityAverage = self.pixelsToMicrons*runningVelocityAverage;
#        runningVelocityAverage_error = self.pixelsToMicrons*(runningVelocityAverage_error/runningVelocityAverage);
#        #traj_frames = self.timePerFrame*traj_frames;
#
#
#        if (plotRodLength == 0):
#
#            A.plotDataSetsWithErrorBars(traj_frames[1:len(traj_frames)], velocityArray, 'Velocity', np.array(None), traj_frames[2:len(traj_frames)], runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', title='Lysis Event', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Velocity (micrometers/second)');
#
#            return
#
#        else:
#            rodLengthArray = self.pixelsToMicrons*BIGLIST[ID][:, 6];
#            
#            #Fit straight line to these lengths
#            M,C = curve_fit(A.fitStraightline, traj_frames[0:len(traj_frames)], rodLengthArray)[0]
#            rodLengthFit = M*traj_frames[0:len(traj_frames)] + C;
#
#            #Plot velocity array and running average velocity with error
#            A.plotDataSetsWithErrorBars(traj_frames[1:len(traj_frames)], velocityArray, 'Velocity', np.array(None), traj_frames[2:len(traj_frames)], runningVelocityAverage, y1_error=runningVelocityAverage_error, label1='Running Average', title='Lysis Event', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Velocity (micrometers/second)');
#            
#            #Plot rod length
#            A.plotDataSetsWithErrorBars(traj_frames[0:len(traj_frames)], rodLengthArray, 'Length', np.array(None), traj_frames[0:len(traj_frames)], rodLengthFit, label1='Fit', title='Bacteria Lengths through trajectory', xlbl='frame ('+str(self.timePerFrame)+' seconds/frame)', ylbl='Rod Length (micrometers)');
#                
#            return
#
#
#    def fitStraightline(self, x, M, C): # this is your 'straight line' y=f(x)
#        return M*x + C
#        

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

    def __init__(self, minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons, minStopTimeThreshold, minStopVelocityThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold, minSwimmingExponent, minTauVals):
        self.timePerFrame = timePerFrame;   #For calculating the average velocity over one frame
        self.NumFramesToAverageOver = NumFramesToAverageOver;   #number of frames over which to calculate the displacement over before taking average.
        self.minTrajLen = minTrajLen;
        self.pixelsToMicrons = pixelsToMicrons;
        self.minStopTimeThreshold = minStopTimeThreshold;       #minimum length of stopping time.
        self.minStopVelocityThreshold = minStopVelocityThreshold    #minimum drop in velocity required to define a stopping/lysis event
        self.initialFrameNum = initialFrameNum;
        self.stoppingVelocityThreshold = stoppingVelocityThreshold;     #fraction drop in average running velocity required to define stopping event.
        self.diffusionThreshold = diffusionThreshold;   #Stopped bacteria still diffuse. Above this threshold, bacteria considered still swimming.
        self.minSwimmingExponent = minSwimmingExponent;     #min value of k in function <r**2> = Ct^k to define swimmer. Ideally, k = 2 for swimmer, k = 1 for diffuser and k = 0 for adherer.
        self.minTauVals = minTauVals;    #minimum number of of tau points used to calculate k for separating swimmers from diffusers.

    def calcAverageVelocitiesForAllTraj(self, A, BIGLIST):        
        '''
        Note: This assumes the particles behave uniformly throughout the course of the trajectory i.e they do not slow down.
        
        Calls:
            - calcAverageParticleVelocity
        
        Returns:
            - average velocity of all particle trajectories in pixels per seconds. 
            - array containing velocities of all trajectories frame by frame
            - array containing displacements of all trajectories frame by frame
        '''
        #Get xPositions and yPositions from BIGLIST
        avVelocityAllTraj = []; #Array of average velocity throughout each trajectory
        displacementArray = [];     #Array of displacements throughout trajectories
        velocityArray = [];   #Array of velocities throughout trajectories
        trajID = [];  
        k_exponentArray = [];
        i = 0;
        for i in range(0, len(BIGLIST)):   # iterate through trajectories
            if(len(BIGLIST[i]) < self.minTrajLen):
                continue;
            else:
                temp_avVelocityAllTraj, temp_displacementArray, temp_velocityArray, k_exponent = A.calcAverageParticleVelocity(A, BIGLIST[i]);
                
                if (k_exponent == -10):
                    print 'curve fit failed for trajID = '+str(i);
                    continue;
                
                elif (k_exponent > self.minSwimmingExponent):
                    avVelocityAllTraj.append(temp_avVelocityAllTraj);
                    displacementArray.append(temp_displacementArray);
                    velocityArray.append(temp_velocityArray);
                    trajID.append(i);
                    k_exponentArray.append(k_exponent);

                else:
                    k_exponentArray.append(k_exponent);
                    continue;

        #Convert lists to numpy arrays
        avVelocityAllTraj = np.asarray(avVelocityAllTraj);
        displacementArray = np.asarray(displacementArray);
        velocityArray = np.asarray(velocityArray);
        k_exponentArray = np.asarray(k_exponentArray);
        
        return avVelocityAllTraj, velocityArray, displacementArray, k_exponentArray 
    

    def calcDistBetweenPoints(self, xPositions, yPositions):
        '''
        Calls: nothing

        Returns: an array containing the radial distance a particle has travelled between successive frames separated by self.NumFramesToAverageOver (units = pixels).
        '''

        #check positions arrays are the same length to be sure
        if (len(xPositions) != len(yPositions)):
            print '\nERROR: len(xPositions) != len(yPositions)\n'
        
        r = np.zeros(int(len(xPositions)/self.NumFramesToAverageOver)-1);
        counter = 0;
        for i in range(0, len(r)):
            r[i] = np.sqrt(abs((xPositions[counter+self.NumFramesToAverageOver] - xPositions[counter])**2 + (yPositions[counter+self.NumFramesToAverageOver] - yPositions[counter])**2)); 
            counter = counter + self.NumFramesToAverageOver;
        
        return r

    
    def calcAverageParticleVelocity(self, A, BIGLIST_traj):
        '''
        Calls:
            - calcDistBetweenPoints

        Returns:
            - average velocity of a single particle (single trajectory) in pixels per seconds
        '''
        xPositions = np.zeros(len(BIGLIST_traj));
        yPositions = np.zeros(len(BIGLIST_traj));
        for j in range(0, len(xPositions)):     # iterate through all tracked frames in this trajectory
            xPositions[j] = BIGLIST_traj[j][2];
            yPositions[j] = BIGLIST_traj[j][3];

        #Calculate distance between points
        displacementArray = A.calcDistBetweenPoints(xPositions, yPositions);
            
        #Calculate average velocity and distance travelled by particle
        AvDisplacement = np.mean(displacementArray);
        AvDisplacement_error = np.std(displacementArray)/np.sqrt(len(displacementArray));   #Standard error on mean
        
        velocityArray = displacementArray/(self.NumFramesToAverageOver*self.timePerFrame);
        AvVelocity = np.mean(velocityArray);
        #AvVelocity = AvDisplacement/(self.NumFramesToAverageOver*self.timePerFrame);
        AvVelocity_error = np.std(velocityArray)/np.sqrt(len(velocityArray));   #Standard deviation (NB: Not standard error on mean)

        #Calculate average mean squared displacement as a function of tau
        meanSquaredDispArray, tauArray = A.calcMeanSquaredDisplacement(xPositions, yPositions);

        #Determine if trajectory is for swimmers, diffusers or adherers
        k_exponent = A.separateDiffusersAndSwimmers(A, meanSquaredDispArray, tauArray);
        
        return AvVelocity, displacementArray, velocityArray, k_exponent


    def calcMeanSquaredDisplacement(self, xPositions, yPositions):
        '''
        Calculated mean squared displacement as a function of tau.
        
        Returns: array of average mean squared velocities with increasing values of tau.
        '''
        #tau = np.arange(len()/self.minTrajLen);
        #self.NumFramesToAverageOver
        minTauRange = self.minTauVals*self.NumFramesToAverageOver;

        meanSquaredDispArray = np.zeros(int(len(xPositions)/minTauRange));
        tauArray = np.zeros(int(len(xPositions)/minTauRange));

        tau = self.NumFramesToAverageOver;
        for j in range(0, int(len(xPositions)/minTauRange)):
            counter = 0;
            temp_array = np.zeros(int(len(xPositions)/tau));

            if (len(xPositions) % tau == 0):
                #tau divides len(xPositions) exactly so will get an index error if don't account for this in defintion of i in for loop.
                for i in range(0, int(len(xPositions)/tau)-1):
                    temp_array[i] = abs((xPositions[counter+tau] - xPositions[counter])**2 + (yPositions[counter+tau] - yPositions[counter])**2);
                    counter = counter + tau;
                
            else:
                #don't need to worry about index error so can get extra tau point in data with this branch.
                for i in range(0, int(len(xPositions)/tau)):
                    temp_array[i] = abs((xPositions[counter+tau] - xPositions[counter])**2 + (yPositions[counter+tau] - yPositions[counter])**2);
                    counter = counter + tau;
            
            meanSquaredDispArray[j] = np.mean(temp_array);
            tauArray[j] = tau;
            tau = tau + self.NumFramesToAverageOver;

        return meanSquaredDispArray, tauArray

#    
#        squaredDisplacementArray = displacementArray**2;
#        cumulativeMeanSquaredDispArray = np.zeros(len(displacementArray));
#        counter=3;
#        cumulativeMeanSquaredDispArray[0] = 0;
#        for i in range(0, len(displacementArray)):
#            
#            cumulativeMeanSquaredDispArray[i] = (1/len(squaredDisplacementArray[0:counter]))*np.sum(squaredDisplacementArray[0:counter]);
#            counter = counter+3;
#            #cumulativeMeanSquaredDispArray[i] = (1/len(squaredDisplacementArray[0:i+1]))*np.sum(squaredDisplacementArray[0:i+1]);
#            #cumulativeMeanSquaredDispArray[i] = (np.sum(displacementArray[0:i+1]))**2;
#
#        return cumulativeMeanSquaredDispArray


    def separateDiffusersAndSwimmers(self, A, meanSquaredDisplacementArray, tauArray):
        
        t = tauArray*self.timePerFrame;
        #t = self.timePerFrame*3*self.NumFramesToAverageOver*np.arange(len(cumulativeMeanSquaredDisplacementArray));
        init_vals = [1, 1];     # for [amp, cen, wid]

        #print tauArray
        #print meanSquaredDisplacementArray

        try:
            best_vals, covar = curve_fit(A.fitDiffusersAndSwimmers, t, meanSquaredDisplacementArray, p0=init_vals);
            #print(best_vals);
            k = best_vals[1];
            C = best_vals[0];
        
        except RuntimeError:
            print("Optimal parameters not found: Number of calls to function has reached maxfev = 600.\n");
            k = -10; 
        
        return k
        
#        from numpy import random
#
#        i = 0;
#        k = [];
#        C = [];
#        while (i < len(im.columns)):
#            try:
#                x = im[im.columns[i]];
#            except KeyError:
#                print("im[im.columns[i]] is not defined. i = "+str(i));
#                i = i+1;
#                continue;
#            else:
#                #print("sure, it was defined. i = "+str(i));
#                x = x[np.logical_not(np.isnan(x))];  #remove NaN from array
#                y = fitDiffusersAndSwimmers(x, 2.33, 0.21) + random.normal(0, 0.2, len(x));
#                init_vals = [1, 1];     # for [amp, cen, wid]
#                best_vals, covar = curve_fit(fitHist, x, y, p0=init_vals);
#                print(best_vals);
#                k.append(best_vals[1])
#                C.append(best_vals[0])
#                i = i+1;
        return
 
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


    ### Calculate historgram to separate diffusers from bacteria stuck to surface
    def fitDiffusersAndSwimmers(self, x, C, k):
        "Equation to fit data: y = Cx^k"
        return (C*(x**k))

    
    # Plot distribution of velocities
    def plotHistogram(self, data, xlbl='bins', ylbl='Frequency', xlim=np.array(None)):
        plt.figure()
	plt.hist(data[:], bins=60)
	plt.xlabel(xlbl);
	plt.ylabel(ylbl);
        if (xlim.all() != None):
            plt.xlim(xlim[0], xlim[1]);

        plt.savefig(outputFile);
	#plt.show()
        return

    # Plot distribution of velocities
    def plotNormalisedHistograms(self, data00, label00=np.array(None), data01=np.array(None), label01=np.array(None), data02=np.array(None), label02=np.array(None), data03=np.array(None), label03=np.array(None), data04=np.array(None), label04=np.array(None), data05=np.array(None), label05=np.array(None), data06=np.array(None), label06=np.array(None), xlbl='bins', ylbl='Normalised Frequency', xlim=np.array(None), saveFilename=None):
        
        binNum = 60;
        alphaVal = 0.4;

        plt.figure()
        data00 = data00/len(data00);
        weights = np.ones_like(data00)/float(len(data00));
        plt.hist(data00[:], weights=weights, label=label00, bins=binNum, alpha=alphaVal);
        
        if (data01.all() == None):
            data01 = data01/len(data01);
            weights = np.ones_like(data01)/float(len(data01));
	    plt.hist(data01[:], weights=weights, label=label01, bins=binNum, alpha=alphaVal);
	
        if (data02.all() == None):
            data02 = data02/len(data02);
            weights = np.ones_like(data02)/float(len(data02));
	    plt.hist(data02[:], weights=weights, label=label02, bins=binNum, alpha=alphaVal);
        
        if (data03.all() == None):
            data03 = data03/len(data03);
            weights = np.ones_like(data03)/float(len(data03));
	    plt.hist(data03[:], weights=weights, label=label03, bins=binNum, alpha=alphaVal);
        
        if (data04.all() == None):
            data04 = data04/len(data04);
            weights = np.ones_like(data04)/float(len(data04));
	    plt.hist(data04[:], weights=weights, label=label04, bins=binNum, alpha=alphaVal);
        
        if (data05.all() == None):
            data05 = data05/len(data05);
            weights = np.ones_like(data05)/float(len(data05));
	    plt.hist(data05[:], weights=weights, label=label05, bins=binNum, alpha=alphaVal);
        
        if (data06.all() == None):
            data06 = data06/len(data06);
            weights = np.ones_like(data01)/float(len(data06));
	    plt.hist(data06[:], weights=weights, label=label06, bins=binNum, alpha=alphaVal);
        
	plt.xlabel(xlbl);
	plt.ylabel(ylbl);
        if (xlim.all() != None):
            plt.xlim(xlim[0], xlim[1]);

        plt.legend(loc='upper right');
        
        if (saveFilename != None):
            plt.savefig(saveFilename);
	
        #plt.show()
        
        return
    
    def plotHistogramWithCurveFit(self, A, data, xlbl='bins', ylbl='Frequency', fit='schulz'):
        popt_RuntimeError = [-10, -10, -10];
	sigma_RuntimeError = -10;

	plt.figure()
        bin_heights, bin_borders, _ = plt.hist(data, bins=60, label='histogram')
        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
        x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000);
        
        if (fit == 'gaussian'):
	    try:
            	popt, _ = curve_fit(A.gaussian, bin_centers, bin_heights, p0=[30., 10., 10.])
            	lbl = 'gaussian fit';
            	y_for_plot = A.gaussian(x_interval_for_fit, *popt); 
            	sigma = popt[2];
           	print 'Gaussian fit with equation: P(v) = %f exp((v-%f/%f)**2) \n' % tuple(popt);
		
	    except RuntimeError:
		print 'RuntimeError: Optimal parameters not found: Number of calls to function has reached maxfev = 800. \n';
		return popt_RuntimeError, sigma_RuntimeError
            

        elif (fit == 'schulz'):
	    try:
            	popt, _ = curve_fit(A.schulz, bin_centers, bin_heights, p0=[1., 2100, 40.])
            	lbl = 'schulz fit';
            	y_for_plot = A.schulz(x_interval_for_fit, *popt); 
            	v_bar = popt[2];
            	z = popt[1];
            	variance = np.sqrt(v_bar**2/(z+1));
            	#sigma = np.sqrt(variance)/np.sqrt(len(y_for_plot));     #Standard error on the mean.
            	sigma = np.sqrt(variance);     #Standard deviation (not standard error on mean).
                print 'Schulz fit with equation: P(v) = (v**z/z!)*((z+1)/v_bar)**(z+1)*np.exp(-(v/v_bar)*(z+1)). Optimised constants: z_factorial = %f, z = %f, v_bar = %f \n' % tuple(popt); 
		
	    except RuntimeError:
		print 'RuntimeError: Optimal parameters not found: Number of calls to function has reached maxfev = 800.\n';
        	return popt_RuntimeError, sigma_RuntimeError

        else:
            print 'ERROR: Undefined distribution for fit\n';
            return

	# Curve fit succeeded so plot and return fit_params, popt and sigma
        plt.plot(x_interval_for_fit, y_for_plot, label=lbl);
        plt.legend();
        plt.xlabel(xlbl);
        plt.ylabel(ylbl);
        
        return popt, sigma

    #Gaussian function
    def gaussian(self, x, x0, y0, sigma):
        p = [x0, y0, sigma]
        return p[1]* np.exp(-((x-p[0])/p[2])**2)

    # Schulz distribution
    def schulz(self, v, z, z_factorial, v_bar):
        p = [z_factorial, z, v_bar];        
        #return p[0]*(v/p[2])**p[1]*np.exp(-(p[1]+1)*(v/p[2]))
        return (v**p[1]/p[0])*((p[1]+1)/p[2])**(p[1]+1)*np.exp(-(v/p[2])*(p[1]+1))

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
#    def plotDataSetsWithErrorBars(self, x0, y0, label0, y0_error=np.array(None), x1=np.array(None), y1=np.array(None), label1=None, y1_error=np.array(None), x2=np.array(None), y2=np.array(None), label2=None, y2_error=np.array(None), x3=np.array(None), y3=np.array(None), label3=None, y3_error=np.array(None), x4=np.array(None), y4=np.array(None), label4=None, y4_error=np.array(None), title=None, xlbl=None, ylbl=None):
    def plotDataSetsWithErrorBars(self, x0, y0, label0, y0_error=np.array(None), x1=np.array(None), y1=np.array(None), y1_error=np.array(None), label1=None, x2=np.array(None), y2=np.array(None), y2_error=np.array(None), label2=None, x3=np.array(None), y3=np.array(None), y3_error=np.array(None), label3=None, x4=np.array(None), y4=np.array(None), y4_error=np.array(None), label4=None, title=None, xlbl=None, ylbl=None):
        plt.figure()
       
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
#        if(y1.all() != None): plt.legend(loc='upper right');
        #plt.legend(loc='upper right');
        #plt.rc('font', family='serif', size=15);
        #plt.savefig(outputPlotName)
        #plt.show()


    # write values to file
    def writeToFile(self, filename, data):
        # ADD IF STATEMENT, CHECK IF FILE EXISTS AND IF NOT CREATE ONE.
        myFile = open(filename,'w')
        myFile.write('i Pdf\n\n');    
        for i in range(0,len(data)): myFile.write('%d           %f\n' % (i, data[i]));
        myFile.close()

