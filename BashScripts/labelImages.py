#from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#from trackReaderClass import readTracks
from toolsForBugs import readCoordinateFile     #to read positions.dat
from toolsForBugs import readTrackingFile       #to read tracks.dat

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

import os

##########################################################################################################
## Define functions
##########################################################################################################

class produceLabelledImages:

    def __init__(self, textSize, initialFrameNum, labelOffset):
        #self.analyseTrajectories = A;   #For calculating the average velocity over one frame
        self.textSize = textSize;
        self.initialFrameNum = initialFrameNum;
        self.labelOffset = labelOffset;   # number of pixels the trajectory label should be offset by in the y direction.
        #self.xPositions = [];
        #self.yPositions = [];
    
    def drawTextOnImage(self, filename, xPos, yPos, textOnImage):
        # Open image
        try:
            base = Image.open(filename).convert('RGBA');
        except:
            print "Unable to load image"
            return

        # make a blank image for the text, initialized to transparent text color
        txt = Image.new('RGBA', base.size, (255,255,255,0))     #colour white

        # get a font
        #fnt = ImageFont.truetype('/usr/share/fonts/truetype/freefont/FreeMono.ttf', self.textSize)
        fnt = ImageFont.truetype('Arial.ttf', self.textSize);   # works on mac
        
        # get a drawing context
        d = ImageDraw.Draw(txt)

        # draw text on image
        d.text((xPos,yPos), textOnImage, font=fnt, fill=(255,255,255,255))

        out = Image.alpha_composite(base, txt)
        out.show()

        return


    # Add trajectory labels to all bacteria in single image
    def labelImage(self, B, filename, BIGLIST, counter):
        frame = counter + self.initialFrameNum;
        #Get position and IDs of trajectories in frame
        IDCounter, xPosArray, yPosArray = B.getTrajInFrame(BIGLIST, frame);
 
        # Open Image
        try:
            base = Image.open(filename).convert('RGBA');
        except:
            print "Unable to load image"
            return

        # make a blank image for the text, initialized to transparent text color
        txt = Image.new('RGBA', base.size, (255,255,255,0))
        
        # get a font
        #fnt = ImageFont.truetype('/usr/share/fonts/truetype/freefont/FreeMono.ttf', self.textSize)
        fnt = ImageFont.truetype('Arial.ttf', self.textSize);   # works on mac

        # get a drawing context
        d = ImageDraw.Draw(txt)
        
        print 'length of IDCounter = '+str(len(IDCounter));

        for i in range(0,len(IDCounter)):
            # draw text on image
            d.text((xPosArray[i],yPosArray[i]), str(int(IDCounter[i])), font=fnt, fill=(255,255,255,255))

        out = Image.alpha_composite(base, txt)
        
        return out


        #Now take trajectories, which presumably have IDs which are given according to the trajectory. Careful, the trajectory ID should be the same throughout the trajectory but don't seem to be in many cases. Print the image given by frame, print all the IDs in IDCounter on this image and save the image in a new directory.

        #So print all IDs on frame
        #save frame to new directory
        #move to next frame.
        # check this before running full loop!

    def getTrajInFrame(self, BIGLIST, frame):
        xPosArray = [];
        yPosArray = [];
        IDCounter = [];
        for i in range(0, len(BIGLIST)):
            for j in range(0, len(BIGLIST[i])):
                if(BIGLIST[i][j, 0] == frame):
                    xPosArray.append(BIGLIST[i][j, 2]);
                    yPosArray.append(BIGLIST[i][j, 3]);
                    IDCounter.append(i)

        IDCounter = np.asarray(IDCounter);
        xPosArray = np.asarray(xPosArray);
        yPosArray = np.asarray(yPosArray);
        return IDCounter, xPosArray, yPosArray

    
    #Only labels the image
    def labelFoundImages(self, B, positionImagesDir, outputImageDir, BIGLIST):
        counter = 0;
        fileList = sorted(os.listdir(positionImagesDir))
        for file in fileList:
            #Ensure only the images beginning with found are outputted
            if (file[0:5] == 'Found'):
                filename = positionImagesDir+file;
                print filename
                annotatedImage = B.labelImage(B, filename, BIGLIST, counter);
                annotatedImage.save(outputImageDir+file)
                counter = counter+1
        
        return

    #Only draw trajectories on image. Could also put trajectory labels?
    def drawTrajOnFoundImages(self, B, positionImagesDir, outputImageDir, BIGLIST):
        fileCounter = 0;
        fileList = sorted(os.listdir(positionImagesDir))
        for file in fileList:
            #Ensure only the images beginning with found are outputted
            if (file[0:5] == 'Found'):
                if (int(file[6:11]) < self.initialFrameNum):
                    continue;

                filename = positionImagesDir+file;
                print filename
                annotatedImage = B.drawTrajectories(B, filename, BIGLIST, fileCounter);
                #annotatedImage = B.labelImage(B, filename, BIGLIST, counter);
                annotatedImage.save(outputImageDir+file)
                fileCounter = fileCounter+1
        
        return
    
    def getPastTrajPositions(self, IDCounter, BIGLIST, frame):
        traj_xPositions = [];
        traj_yPositions = [];

        print 'frame = '+str(frame);
        #print 'IDCounter = '+str(IDCounter);

        for i in range(0, len(IDCounter)):
            #currentIndexInTraj = (int(BIGLIST[IDCounter[i]][0, 0]) - (fileCounter+self.initialFrameNum)) + fileCounter;
            currentIndexInTraj = frame - int(BIGLIST[IDCounter[i]][0, 0]);
            traj_xPositions.append(BIGLIST[IDCounter[i]][0:currentIndexInTraj, 2]);
            traj_yPositions.append(BIGLIST[IDCounter[i]][0:currentIndexInTraj, 3]);

        #print 'traj_xPositions[23]'+str(traj_xPositions[23])
        return traj_xPositions, traj_yPositions


    def drawTrajectories(self, B, filename, BIGLIST, fileCounter):
        frame = fileCounter + self.initialFrameNum;
        #Get position and IDs of trajectories in frame
        IDCounter, xPosArray, yPosArray = B.getTrajInFrame(BIGLIST, frame);
        
        #Get xPositions and yPositions for length of trajectory
        #for i in range(0,len(IDCounter)):     #IDCounter contains all trajectories that exist in the current frame
        traj_xPositions, traj_yPositions = B.getPastTrajPositions(IDCounter, BIGLIST, frame);


        # Open Image
        try:
            base = Image.open(filename).convert('RGBA');
        except:
            print "Unable to load image"
            return

        # make a blank image for the lines and text, initialized to transparent text color
        line = Image.new('RGBA', base.size, (255,255,255,0))
        
        # get a font
        #fnt = ImageFont.truetype('/usr/share/fonts/truetype/freefont/FreeMono.ttf', self.textSize)
        fnt = ImageFont.truetype('Arial.ttf', self.textSize);   # works on mac
        
        # get a drawing context
        d = ImageDraw.Draw(line)
        
        #print 'length of IDCounter = '+str(len(IDCounter));
        #print 'traj_xPositions'+str(traj_xPositions);
        #print 'traj_xPositions[0]'+str(traj_xPositions[0]);
        
        for i in range(0, len(IDCounter)):
            # label bacteria
            #lineColour = (int(255*np.random.uniform()), int(255*np.random.uniform()), 200+int(55*np.random.uniform()), 255)
            lineColour = (255, 255, 255, 255);
            d.text((xPosArray[i],yPosArray[i]+self.labelOffset), str(int(IDCounter[i])), font=fnt, fill=lineColour)
            for j in range(0, len(traj_xPositions[i])-1):
                # draw lines on image
                d.line((traj_xPositions[i][j+1],traj_yPositions[i][j+1],   traj_xPositions[i][j],traj_yPositions[i][j]), fill=lineColour)

        out = Image.alpha_composite(base, line)
        
       # 
       # base = Image.open('a_test.png').convert('RGBA')
       # im = Image.new('RGBA', base.size, (255,255,255,0))
       # draw = ImageDraw.Draw(im)
       # draw.line((100,200, 150,300), fill=(255,255,10,255))
       # out = Image.alpha_composite(base, im)
       # out.show()

        return out

##########################################################################################################
## Main Program begins here
##########################################################################################################

#Declare Variables
NumFramesToAverageOver = 1;
minTrajLen = 30*NumFramesToAverageOver;
fps = 25;
timePerFrame = 1./fps;
pixelsToMicrons = 0.702;    # For x20 Mag
#pixelsToMicrons = 0.354;    # For x40 Mag
textSize = 10;
textColour = 0;


###### MAC ###
#fileDir = '../../../../../../../../Volumes/MyBook/MastersProject/Data/20180227/20180227SurfaceVid1-25fpsx20Mag2000Frames/DDMmovies180227-143931-AsImageSequences/Output-Pos02_Movie0010/';
#fileDir = '../../../../../../../../Volumes/MyBook/MastersProject/Data/20180227/20180227SurfaceVid1-25fpsx20Mag2000Frames/DDMmovies180227-143931-AsImageSequences/Output-Pos03_Movie0011/';
fileDir = '../../../../../../../../Volumes/MyBook/MastersProject/Data/20180227/20180227SurfaceVid2-25fpsx20Mag2000Frames/DDMmovies180227-155825-AsImageSequences/Output-Pos01_Movie0013/';

###### UBUNTU ####
## TEST FILE:
#fileDir = '../../../../../../media/cameron/MyBook/MastersProject/Data/DataForTestingTracking/20180213-x20Surface/';

## Run Files:
#fileDir = '../../../../../../media/cameron/MyBook/MastersProject/Data/20180206/20180206-50fpsx20-SurfacePhage+ControlMeas3/DDMmovies180206-164041-AsImageSequences/Output-Pos01_Movie0001/';
#fileDir = '../../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences/Output-Pos00_Movie0001/';
#fileDir = '../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-154132-AsImageSequences/Output-Pos03_Movie0002/';


trackFilename = fileDir+'filterTracks2DtOutput/tracks_fixed.dat';
positionImagesDir = fileDir+'findRods2DtOutput/';
#outputImageDir = fileDir+'findRods2DtOutput-WithFilteredTrajectories/';
outputImageDir = '../../../../../../../../Volumes/CBOGGONUSB/Data/20180227-Surface/'+'findRods2DtOutput-Vid2-Pos01_Movie0013-WithFilteredTrajectories/';


#Read in tracking data
BIGLIST, numberOfFrames = readTrackingFile(trackFilename)
#COORD_BIGLIST, COORD_PROPERTIES = readCoordinateFile(filename)
#print BIGLIST
#print numberOfFrames

#initialFrameNum = BIGLIST[0][0, 0]
initialFrameNum = 1;
labelOffset = 5;
B = produceLabelledImages(textSize, initialFrameNum, labelOffset);


#BIGLIST_traj = BIGLIST[0];

#IDs = np.zeros(len(BIGLIST_traj));
#xPositions = np.zeros(len(BIGLIST_traj));
#yPositions = np.zeros(len(BIGLIST_traj));
#for j in range(0, len(xPositions)):     # iterate through all tracked frames in this trajectory
#    IDs[j] = BIGLIST_traj[j][1];
#    xPositions[j] = BIGLIST_traj[j][2];
#    yPositions[j] = BIGLIST_traj[j][3];

#### Working with PIL



try:
    os.stat(outputImageDir)
except:
    os.mkdir(outputImageDir);


#B.labelFoundImages(B, positionImagesDir, outputImageDir, BIGLIST);
B.drawTrajOnFoundImages(B, positionImagesDir, outputImageDir, BIGLIST);

#counter = 0;
#fileList = sorted(os.listdir(positionImagesDir))
#for file in fileList:
#    #Ensure only the images beginning with found are outputted
#    if (file[0:5] == 'Found'):
#        filename = positionImagesDir+file;
#        print filename
#        #annotatedImage = B.labelImage(B, filename, BIGLIST, counter);
#        #annotatedImage.save(outputImageDir+file)
#        counter = counter+1
#
#
#annotatedImage = B.labelImage(B, '../Data/171201-DDM/Output-Pos01_Movie0000/findRods2DtOutput/Found_00000.tif', BIGLIST, 0);
#annotatedImage.show();
#out.show()

#B.drawTextOnImage('../Data/171201-DDM/Output-Pos01_Movie0000/findRods2DtOutput/Found_00000.tif', BIGLIST[0][4,2], BIGLIST[0][4,3], str(BIGLIST[0][4,1]))



#lastly need to run:
#ffmpeg -framerate 50 -start_number 00001 -i Found_%4d.tif outputMovie.pf4

