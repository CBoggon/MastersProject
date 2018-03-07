#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit


# Plot distribution of velocities
def plotHistogram(data, xlbl='bins', ylbl='Frequency', title=None):
    plt.hist(data[:], bins=20)
    plt.xlabel(xlbl);
    plt.ylabel(ylbl);
    plt.title(title);
    #plt.savefig('outputHistogram.pdf')
    #plt.show()
    return

# write values to file
def readToFile(filename):
            myFile = open(filename,'r')
            
            #Read data in from file
            with open(filename) as f:
                content = f.readlines()
                # you may also want to remove whitespace characters like `\n` at the end of each line
            
            myFile.close()
            content = [x.strip() for x in content]
            
            #append data to 
            IDArray = [];
            stopFramearray = [];
            timeStoppedArray = [];
            veloctyFileArray = [];
            for i in range(0,len(content)):
                data = content[i].split();
                IDArray.append(data[0]);
                stopFramearray.append(float(data[1]));
                timeStoppedArray.append(float(data[2]));
                veloctyFileArray.append(data[3]);

            stopFramearray = np.asarray(stopFramearray);
            timeStoppedArray = np.asarray(timeStoppedArray);
            return IDArray, stopFramearray, timeStoppedArray, veloctyFileArray


fps = 25;
timePerFrame = 1./fps;

fileDir = '../../../../../../../../Volumes/CBOGGONUSB/Data/20180227-Surface/lysisTrackingOutput/lysisEvents.dat';

IDArray, stopFramearray, timeStoppedArray, veloctyFileArray = readToFile(fileDir);

# Plot timeStopped histogram
timeStoppedArray = timePerFrame*timeStoppedArray;
plotHistogram(timeStoppedArray, xlbl='time (seconds)', ylbl='Frequency', title='Time Bacteria Stopped Before Lysis (9 data points)');

plt.show()
