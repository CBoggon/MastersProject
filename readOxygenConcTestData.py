#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
#from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit
import os

#from analyseDDMTracks import analyseTrajectories


def readImageJFile(filename, skipLines=1):
    myFile = open(filename,'r')

    #Read data in from file
    with open(filename) as f:
        content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line

    myFile.close()

    # Declare arrays to depend data to
    pixels = np.zeros(len(content));
    intensity = np.zeros(len(content));

    for i in range(skipLines, len(content)):
        data = content[i].rstrip().split(',')
        pixels[i] = float(data[0]);
        intensity[i] = float(data[1]);

    return pixels, intensity


def fitStraightLine(x, M, C): # this is your 'straight line' y=f(x)
    return M*x + C

def calcAvIntensity(pixelArray, intensityArray):
    
    AvIntensity = np.zeros(len(intensityArray));
    AvIntensity_error = np.zeros(len(intensityArray));
    for i in range(0, len(pixelArray)):
        #M,C = curve_fit(fitStraightLine, pixelArray[i], intensityArray[i])[0]
        #AvIntensity[i] = C;
        AvIntensity[i] = np.mean(intensityArray[i]);
        AvIntensity_error[i] = np.std(intensityArray[i])/2;  #Standard deviation

    return AvIntensity, AvIntensity_error


def plotDataSetsWithErrorBars(x0, y0, label0, y0_error=np.array(None), x1=np.array(None), y1=np.array(None), y1_error=np.array(None), label1=None, x2=np.array(None), y2=np.array(None), y2_error=np.array(None), label2=None, x3=np.array(None), y3=np.array(None), y3_error=np.array(None), label3=None, x4=np.array(None), y4=np.array(None), y4_error=np.array(None), label4=None, title=None, xlbl=None, ylbl=None, plotLegend=True, saveFileName=None):

    plt.figure(figsize=(10, 6))

    if (y0_error.all() == None):
        plt.plot(x0, y0, 'o-', label=label0)
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
            plt.plot(x3, y3, 's-', label=label3)     #Removed 's-' so does not show data marker on this plot.
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
    if(plotLegend == True): plt.legend(loc='upper left');
    if(saveFileName != None): plt.savefig(saveFileName);

#   if(y1.all() != None): plt.legend(loc='upper right');
    #plt.legend(loc='upper right');
    #plt.rc('font', family='serif', size=13);
    #plt.savefig(outputPlotName)
    #plt.show()

    return


#######################################################################################################
#####    MAIN
#######################################################################################################


#Define filenames
FileDir = '../../Data/171122-Oxygen/ImageJ-IntensityFiles/';

outputSaveFileDir = FileDir+'AnalysisOutput/';


# initialise arrays for data
ControlPixelsArray = [];
ControlIntensityArray = [];
GlassPixelsArray = [];
GlassIntensityArray = [];
PDMSPixelsArray = [];
PDMSIntensityArray = [];

folderList = sorted(os.listdir(FileDir))
for folder in folderList:
    #Ensure only the images beginning with found are outputted
    if (folder[0:7] == 'Control'):
        fileList = sorted(os.listdir(FileDir+'Control/'))
        for file in fileList:
            temp_pixels, temp_intensity = readImageJFile(FileDir+'Control/'+file);
            ControlPixelsArray.append(temp_pixels);
            ControlIntensityArray.append(temp_intensity);
            
    elif (folder[0:5] == 'Glass'):
        fileList = sorted(os.listdir(FileDir+'Glass/'))
        for file in fileList:
            temp_pixels, temp_intensity = readImageJFile(FileDir+'Glass/'+file);
            GlassPixelsArray.append(temp_pixels);
            GlassIntensityArray.append(temp_intensity);

    elif (folder[0:4] == 'PDMS'):
        fileList = sorted(os.listdir(FileDir+'PDMS/'))
        for file in fileList:
            temp_pixels, temp_intensity = readImageJFile(FileDir+'PDMS/'+file);
            PDMSPixelsArray.append(temp_pixels);
            PDMSIntensityArray.append(temp_intensity);
    
    else:
        continue;


ControlPixelsArray = np.asarray(ControlPixelsArray);
ControlIntensityArray = np.asarray(ControlIntensityArray);
GlassPixelsArray = np.asarray(GlassPixelsArray);
GlassIntensityArray = np.asarray(GlassIntensityArray);
PDMSPixelsArray = np.asarray(PDMSPixelsArray);
PDMSIntensityArray = np.asarray(PDMSIntensityArray);

## Now perfom straight line fit to find average intensity at each data point
time_Control = np.array([30, 33, 37]);
time_Glass = np.array([19, 22, 25, 29, 32, 43]);
time_PDMS = np.array([50, 53, 57, 60, 63, 66, 68, 71, 74, 77, 80, 82]);

time_Control = time_Control - time_Control[0];
time_Glass = time_Glass - time_Glass[0];
time_PDMS = time_PDMS - time_PDMS[0];

#time_Control = np.arange(len(ControlPixelsArray));
#time_Glass = np.arange(len(GlassPixelsArray));
#time_PDMS = np.arange(len(PDMSPixelsArray));

AvIntensity_Control, AvIntensity_Control_error = calcAvIntensity(ControlPixelsArray, ControlIntensityArray);
AvIntensity_Glass, AvIntensity_Glass_error = calcAvIntensity(GlassPixelsArray, GlassIntensityArray);
AvIntensity_PDMS, AvIntensity_PDMS_error = calcAvIntensity(PDMSPixelsArray, PDMSIntensityArray);


#plotDataSetsWithErrorBars(time_Glass, AvIntensity_Glass, None, AvIntensity_Glass_error);

plotDataSetsWithErrorBars(time_Glass, AvIntensity_Glass, 'Glass', AvIntensity_Glass_error, x1=time_Control, y1=AvIntensity_Control, y1_error=AvIntensity_Control_error, label1='PDMS-Control', x2=time_PDMS, y2=AvIntensity_PDMS, y2_error=AvIntensity_PDMS_error, label2='PDMS', xlbl='Time (minutes)', ylbl='Intensity (arbitrary units)', plotLegend=False);
plt.legend(loc='upper right')

#
## Plot DICF vs tau:
#tau, gArray0, gArray1, gArray2
#gPlotTimes = [4, 23, 26];
#A.plotDataSetsWithErrorBars(tau, gArray0, 't = '+str(gPlotTimes[0])+' mins', x1=tau, y1=gArray1, label1='t = '+str(gPlotTimes[1])+' mins', x2=tau, y2=gArray2, label2='t = '+str(gPlotTimes[2])+' mins', xlbl='Delay Time, tau (seconds)', ylbl='g(q, tau)', plotLegend=False);
#plt.legend(loc='upper left')
#plt.xscale('log')
##plt.ylim([-0,1.1])
#outputFile = outputSaveFileDir+'DDMGVsQ';
#plt.savefig(outputFile);
#

plt.show()
