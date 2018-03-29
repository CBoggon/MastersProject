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

def plotDataSetsWithUpperAndLowerErrorBars(x0, y0, label0, y0_error_lower=np.array(None), y0_error_upper=np.array(None), x1=np.array(None), y1=np.array(None), y1_error_lower=np.array(None), y1_error_upper=np.array(None), label1=None, x2=np.array(None), y2=np.array(None), y2_error_lower=np.array(None), y2_error_upper=np.array(None), label2=None, x3=np.array(None), y3=np.array(None), y3_error_lower=np.array(None), y3_error_upper=np.array(None), label3=None, x4=np.array(None), y4=np.array(None), y4_error_lower=np.array(None), y4_error_upper=np.array(None), label4=None, title=None, xlbl=None, ylbl=None, saveFileName=None, plotLegend=True):

    plt.figure(figsize=(10, 6))
    plt.rc('font', family='serif', size=15);

    if (y0_error_lower.all() == None):
        plt.plot(x0, y0, 'x', label=label0)
    else:
        plt.errorbar(x0, y0, yerr=(y0_error_lower, y0_error_upper), fmt='x', label = label0)

    if(y1.all() != None):
        if (y1_error_lower.all() == None):
            plt.plot(x1, y1, 'o', label=label1)
        else:
            plt.errorbar(x1, y1, yerr=(y1_error_lower, y0_error_upper), fmt='o', label = label1)

    if(y2.all() != None):
        if (y2_error_lower.all() == None):
            plt.plot(x2, y2, '^', label=label2)
        else:
            plt.errorbar(x2, y2, yerr=(y2_error_lower, y0_error_upper), fmt='^', label = label2)

    if(y3.all() != None):
        if (y3_error_lower.all() == None):
            plt.plot(x3, y3, 's', label=label3)     #Removed 's-' so does not show data marker on this plot.
        else:
            plt.errorbar(x3, y3, yerr=(y3_error_lower, y0_error_upper), fmt='s', label = label3)

    if(y4.all() != None):
        if (y4_error_lower.all() == None):
            plt.plot(x4, y4, 'd', label=label4)
        else:
            plt.errorbar(x4, y4, yerr=(y4_error_lower, y0_error_upper), fmt='d', label = label4)

    if(title != None): plt.suptitle(title);
    if(xlbl != None): plt.xlabel(xlbl);
    if(ylbl != None): plt.ylabel(ylbl);
    if(plotLegend == True): plt.legend(loc='upper left');
    #if (xlim.all() != None):
    plt.ylim(0, 26);
    if(saveFileName != None): plt.savefig(saveFileName);
    plt.xticks(np.linspace(0,1,0.5), np.array(['0' for i in range(0,1)]));
    plt.rc('font', family='serif', size=15);
    #plt.show()

    return

def plotDataSetsWithErrorBars(x0, y0, label0, y0_error=np.array(None), x1=np.array(None), y1=np.array(None), y1_error=np.array(None), label1=None, x2=np.array(None), y2=np.array(None), y2_error=np.array(None), label2=None, x3=np.array(None), y3=np.array(None), y3_error=np.array(None), label3=None, x4=np.array(None), y4=np.array(None), y4_error=np.array(None), label4=None, title=None, xlbl=None, ylbl=None, plotLegend=True):

    plt.figure(figsize=(10, 6))
    plt.rc('font', family='serif', size=15);

    if (y0_error.all() == None):
        plt.plot(x0, y0, 'o', label=label0)
    else:
        plt.errorbar(x0, y0, yerr=y0_error, fmt='o', label = label0)

    if(y1.all() != None):
        if (y1_error.all() == None):
            plt.plot(x1, y1, 'x', label=label1)
        else:
            plt.errorbar(x1, y1, yerr=y1_error, fmt='x', label = label1)

    if(y2.all() != None):
        if (y2_error.all() == None):
            plt.plot(x2, y2, '^', label=label2)
        else:
            plt.errorbar(x2, y2, yerr=y2_error, fmt='^', label = label2)

    if(y3.all() != None):
        if (y3_error.all() == None):
            plt.plot(x3, y3, 's', label=label3)     #Removed 's-' so does not show data marker on this plot.
        else:
            plt.errorbar(x3, y3, yerr=y3_error, fmt='s', label = label3)

    if(y4.all() != None):
        if (y4_error.all() == None):
            plt.plot(x4, y4, 'd', label=label4)
        else:
            plt.errorbar(x4, y4, yerr=y4_error, fmt='d', label = label4)

    if(title != None): plt.suptitle(title);
    if(xlbl != None): plt.xlabel(xlbl);
    if(ylbl != None): plt.ylabel(ylbl);
    if(plotLegend == True): plt.legend(loc='upper left');
    if(saveFileName != None): plt.savefig(saveFileName);
    
#        if(y1.all() != None): plt.legend(loc='upper right');
    #plt.legend(loc='upper right');
    #plt.rc('font', family='serif', size=13);
    #plt.savefig(outputPlotName)
    #plt.show()

    return


# write values to file
def readFile(filename):
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
            stopFramearray_Lerror = [];
            stopFramearray_Uerror = [];
            veloctyFileArray = [];
            for i in range(0,len(content)):
                data = content[i].split();
                IDArray.append(data[0]);
                stopFramearray.append(float(data[1]));
                timeStoppedArray.append(float(data[2]));
                stopFramearray_Lerror.append(float(data[3]));
                stopFramearray_Uerror.append(float(data[4]));
                veloctyFileArray.append(data[5]);

            stopFramearray = np.asarray(stopFramearray);
            timeStoppedArray = np.asarray(timeStoppedArray);
            stopFramearray_Lerror = np.asarray(stopFramearray_Lerror);
            stopFramearray_Uerror = np.asarray(stopFramearray_Uerror);
            
            # Calc position of upper and lower error in time
            timeStoppedArray_Uerror = (1./25)*(stopFramearray - stopFramearray_Lerror);
            timeStoppedArray_Lerror = (1./25)*(stopFramearray_Uerror - stopFramearray);
            
            # Calc error wrt timeStopped
            #timeStoppedArray_Uerror = timeStoppedArray_Uerror - timeStoppedArray;
            #timeStoppedArray_Lerror = timeStoppedArray - timeStoppedArray_Lerror;
            
            return IDArray, stopFramearray, timeStoppedArray, timeStoppedArray_Lerror, timeStoppedArray_Uerror, veloctyFileArray, stopFramearray_Lerror, stopFramearray_Uerror


fps = 25;
timePerFrame = 1./fps;

fileDir = '../../../../../../../../Volumes/CBOGGONUSB/Data/20180227-Surface/lysisTrackingOutput/lysisEvents.dat';

IDArray, stopFramearray, timeStoppedArray, timeStoppedArray_Lerror, timeStoppedArray_Uerror, veloctyFileArray, stopFramearray_Lerror, stopFramearray_Uerror = readFile(fileDir);

# Plot timeStopped histogram
timeStoppedArray = timePerFrame*timeStoppedArray;
#timeStoppedError = NA;
plotHistogram(timeStoppedArray, xlbl='time (seconds)', ylbl='Frequency', title='Time Bacteria Stopped Before Lysis (9 data points)');

x_onGraph = 1./len(timeStoppedArray)*np.arange(len(timeStoppedArray));
#plotDataSetsWithErrorBars(x_onGraph, timeStoppedArray, None, xlbl='Lysis Event', ylbl='Time Stopped (seconds)', title='Time Bacteria Stopped Before Lysis (9 data points)');
plotDataSetsWithUpperAndLowerErrorBars(x_onGraph, timeStoppedArray, None, timeStoppedArray_Lerror, timeStoppedArray_Uerror, xlbl='Lysis Event', ylbl='Time Stopped (seconds)', saveFileName='../../../../../../../../Volumes/CBOGGONUSB/Data/20180227-Surface/lysisTrackingOutput/LysisStopTimes');

#title='Time Bacteria Stopped Before Lysis'

plt.show()
