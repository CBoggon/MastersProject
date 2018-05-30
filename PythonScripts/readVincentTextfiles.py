#!python
#
###############################################################################
#Program reads in text files containing data from DDM (given to me by Vincent Martinez) and plots data.
#
###############################################################################

#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
#from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat
from scipy.optimize import curve_fit

from analyseDDMTracks import analyseTrajectories

# Schulz distribution
#def schulz(v, z, z_factorial, v_bar):
#    p = [z_factorial, z, v_bar];
#    #return p[0]*(v/p[2])**p[1]*np.exp(-(p[1]+1)*(v/p[2]))
#    return (v**p[1]/p[0])*((p[1]+1)/p[2])**(p[1]+1)*np.exp(-(v/p[2])*(p[1]+1))

def schulzSwimmer(z, z_factorial, v_bar, tau, q):
    delta = (q*v_bar*tau)/(z+1);
    f = ((z+1)/(q*z*v_bar*tau))*(np.sin(z*np.arctan(delta))/((1+delta**2)**(z/2)));
    return f

def produceISFExamples(A):
    D = 0.3;
    q = 1.;
    v_bar = 15.;
    sigma = 7.5;
    beta = 0.3;
    z = 2;
    z_factorial = 2;

    tau = np.arange(0, 40, 0.01);
    
    f_diff = np.exp(-D*tau*q**2);
    f_swim = np.sin(q*v_bar*tau)/(q*v_bar*tau);
    
    #v = np.arange(0,40, 1);
#    #f_temp = [];
#    f_temp = np.zeros(len(f_swim));
#    iterations = 20
#    for i in range(0,iterations):
#        v = 40*np.random.normal();
#        P = schulz(v, z, z_factorial, v_bar)
#        #f_temp.append(P*f_swim);
#        f_temp = f_temp + P*f_swim;
#    
#    #f_temp = np.asarray(f_temp);
#    f_swimAv = f_temp/float(iterations);
#    f_swimAv = f_swimAv/f_swimAv[1];
#    print f_swimAv
    
    f_swimAv = schulzSwimmer(z, z_factorial, v_bar, tau, q);
    f_mix = f_diff*((1-beta) + beta*f_swimAv);

    g_diff = 1 - f_diff;
    g_swim = 1 - f_swim;
    g_swimAv = 1 - f_swimAv;
    g_mix = 1 - f_mix;

    #A.plotDataSetsWithErrorBars(tau, f_diff, 'Diffusers', x1=tau, y1=f_swim, label1='Swimmers', x2=tau, y2=f_mix, label2=str(int(100*beta))+':'+str(int(100*(1-beta)))+' swimmers:diffusers', xlbl='Delay time (seconds)', ylbl='f(q, tau)', plotLegend=True);

    plt.figure(figsize=(12, 7))
    plt.plot(tau, f_diff, '-', label='Diffusers')
    plt.plot(tau, f_swim, '-', label='Single Swimmer')
    plt.plot(tau, f_swimAv, '-', label='Average Swimmers')
    plt.plot(tau, f_mix, '-', label=str(int(100*beta))+':'+str(int(100*(1-beta)))+' swimmers:diffusers')
    plt.xlabel('Delay time (seconds)')
    plt.ylabel('f(q, tau)')
    plt.xscale('log')
    plt.rc('font', family='serif', size=15);
    plt.legend(loc='upper right')


    plt.figure(figsize=(12, 7))
    plt.plot(tau, g_diff, '-', label='Diffusers')
    plt.plot(tau, g_swim, '-', label='Single Swimmer')
    plt.plot(tau, g_swimAv, '-', label='Average Swimmers')
    plt.plot(tau, g_mix, '-', label=str(int(100*beta))+':'+str(int(100*(1-beta)))+' swimmers:diffusers')
    plt.xlabel('Delay time (seconds)')
    plt.ylabel('g(q, tau)')
    plt.xscale('log')
    plt.rc('font', family='serif', size=15);
    plt.legend(loc='lower right')
    
    return

def readGVsQFile(filename, skipLines=1):
    myFile = open(filename,'r')

    #Read data in from file
    with open(filename) as f:
        content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line

    myFile.close()

    # Declare arrays to depend data to
    tau = [];
    gArray0 = [];
    gArray1 = [];
    gArray2 = [];

    i = skipLines*4;
    data = content[0].split();
    while (i < len(data)):
    #for i in range(skipLines,len(content)):
        tau.append(float(data[i]));
        gArray0.append(float(data[i+1]));
        gArray1.append(float(data[i+2]));
        gArray2.append(float(data[i+3]));
        i = i+4;

    tau = np.asarray(tau);
    gArray0 = np.asarray(gArray0);
    gArray1 = np.asarray(gArray1);
    gArray2 = np.asarray(gArray2);

    return tau, gArray0, gArray1, gArray2


def readFile(filename, skipLines=1):
    myFile = open(filename,'r')

    #Read data in from file
    with open(filename) as f:
        content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line

    myFile.close()
    #content = [x.strip() for x in content]

    #append data to arrays
    AvVelocityArray = [];
    AvVelocityArray_error = [];
    NonMotileFractionArray = [];
    NonMotileFractionArray_error = [];
    Time = [];
    NormalisedAmplitudeArray = [];
    NormalisedAmplitudeArray_error = [];
    
    i = skipLines*11;
    data = content[0].split();
    while (i < len(data)):
    #for i in range(skipLines,len(content)):
        AvVelocityArray.append(float(data[i+1]));
        AvVelocityArray_error.append(float(data[i+2]));
        NonMotileFractionArray.append(float(data[i+3]));
        NonMotileFractionArray_error.append(float(data[i+4]));
        Time.append(float(data[i+6]));
        NormalisedAmplitudeArray.append(float(data[i+8]));
        NormalisedAmplitudeArray_error.append(float(data[i+9]));

        i = i+11;

    AvVelocityArray = np.asarray(AvVelocityArray);
    AvVelocityArray_error = np.asarray(AvVelocityArray_error);
    NonMotileFractionArray = np.asarray(NonMotileFractionArray);
    NonMotileFractionArray_error = np.asarray(NonMotileFractionArray_error);
    Time = np.asarray(Time);
    NormalisedAmplitudeArray = np.asarray(NormalisedAmplitudeArray);
    NormalisedAmplitudeArray_error = np.asarray(NormalisedAmplitudeArray_error);

    return AvVelocityArray, AvVelocityArray_error, NonMotileFractionArray, NonMotileFractionArray_error, Time, NormalisedAmplitudeArray, NormalisedAmplitudeArray_error



#######################################################################################################
#####    MAIN
#######################################################################################################

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
minSwimmingExponent = 1.3;
minTauVals = 2;
BacteriaCounterFrame = 200.;


### Initialise Array
A = analyseTrajectories(minTrajLen, NumFramesToAverageOver, timePerFrame, pixelsToMicrons,  minStopTimeThreshold, minStopVelocityThreshold, initialFrameNum, stoppingVelocityThreshold, diffusionThreshold, minSwimmingExponent, minTauVals, BacteriaCounterFrame);



#Define filenames
outputSaveFileDir = '../../Data/180202DDMVincentData/more_data/';
Filename_Pos00 = outputSaveFileDir+'Table_Pos00.txt';
Filename_Pos01 = outputSaveFileDir+'Table_Pos01.txt';
Filename_Pos02 = outputSaveFileDir+'Table_Pos02.txt';
Filename_gVsq = outputSaveFileDir+'gqtau_data_q0p716.txt';


### Plot example DDM plots of ISF
produceISFExamples(A);


## Read data 
AvVelocityArray_Pos00, AvVelocityArray_error_Pos00, NonMotileFractionArray_Pos00, NonMotileFractionArray_error_Pos00, Time_Pos00, NormalisedAmplitudeArray_Pos00, NormalisedAmplitudeArray_error_Pos00 = readFile(Filename_Pos00);

AvVelocityArray_Pos01, AvVelocityArray_error_Pos01, NonMotileFractionArray_Pos01, NonMotileFractionArray_error_Pos01, Time_Pos01, NormalisedAmplitudeArray_Pos01, NormalisedAmplitudeArray_error_Pos01 = readFile(Filename_Pos01);

AvVelocityArray_Pos02, AvVelocityArray_error_Pos02, NonMotileFractionArray_Pos02, NonMotileFractionArray_error_Pos02, Time_Pos02, NormalisedAmplitudeArray_Pos02, NormalisedAmplitudeArray_error_Pos02 = readFile(Filename_Pos02);

tau, gArray0, gArray1, gArray2 = readGVsQFile(Filename_gVsq);

# Calculate motile fraction
MotileFractionArray_Pos00 = 1 - NonMotileFractionArray_Pos00;
MotileFractionArray_Pos01 = 1 - NonMotileFractionArray_Pos01;
MotileFractionArray_Pos02 = 1 - NonMotileFractionArray_Pos02;

# Plot lysis line
lysisLine_x = np.array([0.5 for i in range(0,10)]);
lysisLine_y = np.array([5.0*i for i in range(0,10)]);
Normalised_lysisLine_y = np.linspace(0, 1, 10);


# Plot average speed:
A.plotDataSetsWithErrorBars(Time_Pos00, AvVelocityArray_Pos00, 'Control', y0_error=AvVelocityArray_error_Pos00, x1=Time_Pos01, y1=AvVelocityArray_Pos01, y1_error=AvVelocityArray_error_Pos01, label1='Phage 1', x2=Time_Pos02, y2=AvVelocityArray_Pos02, y2_error=AvVelocityArray_error_Pos02, label2='Phage 2', x3=lysisLine_x, y3=lysisLine_y, label3='Phage Added', xlbl='Time (Minutes)', ylbl='Average Velocity (micrometers/second)', plotLegend=False);
plt.legend(loc='lower left')
plt.ylim([-2,50])
outputFile = outputSaveFileDir+'DDMAvVelocityVsTime';
plt.savefig(outputFile);


# Plot normalised amplitude vs time:
A.plotDataSetsWithErrorBars(Time_Pos00, NormalisedAmplitudeArray_Pos00, 'Control', y0_error=NormalisedAmplitudeArray_error_Pos00, x1=Time_Pos01, y1=NormalisedAmplitudeArray_Pos01, y1_error=NormalisedAmplitudeArray_error_Pos01, label1='Phage 1', x2=Time_Pos02, y2=NormalisedAmplitudeArray_Pos02, y2_error=NormalisedAmplitudeArray_error_Pos02, label2='Phage 2', x3=lysisLine_x, y3=4*Normalised_lysisLine_y, label3='Phage Added', xlbl='Time (Minutes)', ylbl='Normalised Amplitude, A(q)', plotLegend=False);
plt.legend(loc='upper center')
#plt.ylim([-2,50])
outputFile = outputSaveFileDir+'DDMNAqVsTime';
plt.savefig(outputFile);



# Remove specific data points that are very noisy:
deletePosition = 13;
NonMotileFractionArray_Pos02_del = np.delete(NonMotileFractionArray_Pos02, deletePosition)
NonMotileFractionArray_error_Pos02_del = np.delete(NonMotileFractionArray_error_Pos02, deletePosition)
Time_Pos02_del = np.delete(Time_Pos02, deletePosition)

MotileFractionArray_Pos02_del = np.delete(MotileFractionArray_Pos02, deletePosition)



# Plot motile fraction vs time:
#A.plotDataSetsWithErrorBars(Time_Pos00, NonMotileFractionArray_Pos00, 'Control', y0_error=NonMotileFractionArray_error_Pos00, x1=Time_Pos01, y1=NormalisedAmplitudeArray_Pos01, y1_error=NormalisedAmplitudeArray_error_Pos01, label1='A(q)', x2=Time_Pos01, y2=NonMotileFractionArray_Pos02_del, y2_error=NonMotileFractionArray_error_Pos02_del, label2='Phage 2', x3=lysisLine_x, y3=Normalised_lysisLine_y, label3='Phage Infection', xlbl='Time (Minutes)', ylbl='Non-Motile Fraction', plotLegend=False);
#plt.legend(loc='lower center')
#plt.ylim([-0,1.1])
#outputFile = outputSaveFileDir+'DDMNonMotileFractionVsTime';
#plt.savefig(outputFile);


# Plot motile fraction vs time:
A.plotDataSetsWithErrorBars(Time_Pos00, MotileFractionArray_Pos00, 'Control', y0_error=NonMotileFractionArray_error_Pos00, x1=Time_Pos01, y1=MotileFractionArray_Pos01, y1_error=NonMotileFractionArray_error_Pos01, label1='Phage 1', x2=Time_Pos02_del, y2=MotileFractionArray_Pos02_del, y2_error=NonMotileFractionArray_error_Pos02_del, label2='Phage 2', x3=lysisLine_x, y3=Normalised_lysisLine_y, label3='Phages Added', xlbl='Time (Minutes)', ylbl='Motile Fraction', plotLegend=False);
plt.legend(loc='lower left')
plt.ylim([-0,1.1])
outputFile = outputSaveFileDir+'DDMMotileFractionVsTime';
plt.savefig(outputFile);


# Plot DICF vs tau:
tau, gArray0, gArray1, gArray2
gPlotTimes = [4, 23, 26];
A.plotDataSetsWithErrorBars(tau, gArray0, 't = '+str(gPlotTimes[0])+' mins', x1=tau, y1=gArray1, label1='t = '+str(gPlotTimes[1])+' mins', x2=tau, y2=gArray2, label2='t = '+str(gPlotTimes[2])+' mins', xlbl='Delay Time, tau (seconds)', ylbl='g(q, tau)', plotLegend=False);
plt.legend(loc='upper left')
plt.xscale('log')
#plt.ylim([-0,1.1])
outputFile = outputSaveFileDir+'DDMGVsQ';
plt.savefig(outputFile);



# Plot v, Naq and motile-fraction for Phage 1:
A.plotDataWithTwoYAxes(Time_Pos01, AvVelocityArray_Pos01, 'Average Velocity', y0_L_error=AvVelocityArray_error_Pos01, x1_L=lysisLine_x, y1_L=lysisLine_y, label1_L='Phage Added', x0_R=Time_Pos01, y0_R=MotileFractionArray_Pos01, y0_R_error=NonMotileFractionArray_error_Pos01, label0_R='Motile Fraction', x1_R=Time_Pos01, y1_R=NormalisedAmplitudeArray_Pos01, y1_R_error=NormalisedAmplitudeArray_Pos01, label1_R='Normalised A(q)', xlbl='Time (Minutes)', ylbl_L='Swimming Speed (micrometers per second)', ylbl_R='Motile Fraction or Normalised A(q)', plotLegend=False);
plt.legend(loc='center right')
plt.ylim([-0,1.1])
outputFile = outputSaveFileDir+'DDMPhage1AllData';
plt.savefig(outputFile);


plt.show()
