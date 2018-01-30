#!python
#tools for bugs is written to get started with reading files containing found positions and orientations in frames and trajectories
#available under the MIT license

import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from matplotlib import *
from matplotlib.pyplot import *

def setplotfont(s,s2=None,ax=None ):
    if ax==None: ax=gca()
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(s)
    if s2==None: s2=s
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(s2)

def readTrackingFile(filename):
        """Reads a file containing trajectories
            
            Parameters:
                filename - input file
                
            Returns:
		listForAllTrajectories - list containing all trajectories
		numberOfFrames - number of frames in original dataset
        """

	FILE = open(filename,"r") # open file as "read"
	firstLine=FILE.readline()
	coor=firstLine.split(' ') # split into two parts
        numberOfTracks = int(coor[0][1:])
        numberOfFrames = int(coor[1][1:]) # casts part (from 1 to end) of the string firstLine to an integer number

        print 'Number of tracks: ', numberOfTracks # plots the number of trajectories
        print 'Number of frames: ', numberOfFrames # plots the number of frames
	listForAllTrajectories = []
        PROPERTIES = []
        for i in range(numberOfTracks):
		headerLine=FILE.readline()
                coor=headerLine.split(' ') # split into two parts
                trackNumber=int(coor[0][1:])
                numberOfFramesInThisTrack=int(coor[1][1:])
		
		M=np.zeros((numberOfFramesInThisTrack,13))
                for j in range(numberOfFramesInThisTrack):
                        coordinates=FILE.readline()
                        coor=coordinates.split(' ')
                        thistime = float(coor[0])
			someid = float(coor[1])
                        xvalue = float(coor[2])
                        yvalue = float(coor[3])
                        xor = float(coor[4])
                        yor = float(coor[5])
                        L = float(coor[6])
                        fitpar1 = float(coor[7])
                        fitpar2 = float(coor[8])
                        fitpar3 = float(coor[9])
                        width = float(coor[10])
                        fitpar4 = float(coor[11])
			V=array([xor, yor])
			n=np.linalg.norm(V)
			L=n*2.0
			if (n != 0.0):
				V=V/n
	
                        M[j,0] = thistime # number of frames since start of the video
                        M[j,1] = someid # Some id number
                        M[j,2] = xvalue # x position in pixels
                        M[j,3] = yvalue # y position in pixels
                        M[j,4] = V[0] # normalised orientation x-component
                        M[j,5] = V[1] # normalised orientation y-component
			M[j,6] = L     # length in pixels
			M[j,7] = fitpar1   
			M[j,8] = fitpar2    
			M[j,9] = fitpar3  
			M[j,10] = width   
			M[j,11] = fitpar4
                        M[j,12] = trackNumber  # trajectory id number

                listForAllTrajectories.append(M) # add positions for this frame to the list.
	return listForAllTrajectories, numberOfFrames

def readCoordinateFile(filename, correctionslope=0.0, cut1=0.0, cut2=1000000000.0):
	"""Reads a coordinate file
	    
	    Parameters:
		filename - input file
		correctionslope - if the image is tilted, this can be used to rotate the coordinates
		
	    Returns:
		positionsForAllFrames - list containing coordinates per frame
		propertiesForAllFrames - properties of each frame	
	"""

    	FILE = open(filename,"r") # open file as "read"

	firstLine=FILE.readline()
	numberOfFrames = int(firstLine[1:])  # casts part (from 1 to end) of the string firstLine to an integer number
	print 'Number of frames:', numberOfFrames # plots the number of frames
	
	positionsForAllFrames = []
	propertiesForAllFrames = []
	for i in range(numberOfFrames): # loops through all the frames 	
		headerLine=FILE.readline()
		coor=headerLine.split(' ') # split into two parts
		numberOfTheFrame=int(coor[0][1:])
		numberOfBacteriaInTheFrame=int(coor[1][1:])
		dimensionLine=FILE.readline()
		coor=dimensionLine.split(' ') # split into multiple parts
		xmax=float(coor[0]) 
		ymax=float(coor[1])
		M=np.zeros((numberOfBacteriaInTheFrame,12))
		for j in range(numberOfBacteriaInTheFrame):		
			coordinates=FILE.readline()
			coor=coordinates.split(' ')
			xvalue = float(coor[0])
			yvalue = float(coor[1])
			zvalue = float(coor[2])
			xpol = float(coor[3])
			ypol = float(coor[4])
			zpol = float(coor[5])
			fitpar1 = float(coor[6])
			fitpar2 = float(coor[7])
			fitpar3 = float(coor[8])
			width = float(coor[9])
			length = float(coor[10])
			fitpar4 = float(coor[11])
			M[j,0] = xvalue					# x-coordinate COM
			M[j,1] = yvalue-xvalue*correctionslope		# y-coordinate COM
			M[j,2] = zvalue					# z-coordinate (0 for 3D data set)
			M[j,3] = xpol					# x-coordinate of one of the poles
			M[j,4] = ypol					# y-coordinate of one of the poles
			M[j,5] = zpol					# z-coordinate of one of the poles
			M[j,6] = fitpar1				# fitting parameter
			M[j,7] = fitpar2				# fitting parameter
			M[j,8] = fitpar3				# fittin parameter
			M[j,9] = width					# width of the rod	
			M[j,10] = length				# length of the rod
			M[j,11] = fitpar4				# fitting parameter
		
		cut3 = ymax - cut2
		cut4 = ymax - cut1
		M = M[(M[:,1]>cut3)*(M[:,1]<cut4)]
		M[:,1] = M[:,1]-cut3

		positionsForAllFrames.append(M)
		S = np.shape(M)		
	
		ymax = cut4-cut3
		propertiesForAllFrames.append([xmax, ymax, S[0]]) # max dimension x (pixels), max dimension y (pixels), number of bacteria in frame
	
    	FILE.close()
	return positionsForAllFrames, propertiesForAllFrames

def makeHist(y,x, COL, TYPE, NORM, LW, LAB):
	"""Plots a histogram
            
            Parameters:
                x - input bins
		y - input values
		COL - color
		TYPE = line style
		NORM = normalised or now
		LW = line width
		LAB = label
        """

	if (NORM):
		dx = log(x[1])-log(x[0])
		S= sum(y*dx)
		y=y/S

	x2=np.ones(len(x)*2)
	x2[::2]=x;
	x2[1::2]=x
	y2=np.ones(len(y)*2+2)
	y2[2::2]=y;
	y2[1:-1:2]=y
	y2[0]=0;
	y2[-1]=0

	plot(x2,y2, color=COL, ls=TYPE, lw=LW, label=LAB)

def plotPixelBias(X):
	"""Plots the pixel bias in a dataset
            
            Parameters:
                X - list of coordinates (in units of pixels), lengths, etc. to plot the bias for.
        """

        close('all')
        figure
	
	d=0.02
        R = arange(0.0, 1.0+d, d)
        y2,x2=histogram(mod(X[:,0], 1.0), bins = R, normed=1)
	makehist(y2, x2, 'red', '-', 0, 1,'data in paper (met code voor opschonen)')

	legend()	
