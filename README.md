
*Purpose:* 
Code written for masters project which aimed to use DDM and particle tracking to analyse the swimming behaviour of T4 bacteriophage infected K-12 AB1157 E. Coli.


*Design Architecture:*

*General outline:*
Videos from microscope were given as .avi files and Teun Visser's tracking code (found as ENTER-URL) requires a tif image sequence as inputs so in folder BASHSCRIPTS are the files, written in Bash to convert avi files into image sequences. Also in here is a the python program labelImages.py which creates an image sequence with the tracked trajectories of the bacteria.

The file 'runTracking.sh' automatically runs all the tracking steps of Teun's tracking code and arranges them in the right folders.

Folder PythonScripts contains all the python scripts used for analysing the tracked trajcetories. They can be subdivided into two categories: 
1) Code used for performing tracking on the same videos as I used in the DDM analysis. These videos were focussed in the bulk of the bacterial sample (not on a class surface)because DDM assumes bacteria travel in straight lines and bacteria on surfaces swim in cicular trajectories. As the tracking code projects into 2D, the results of the tracking outputted by this analysis did not agree with the DDM results (See results section in report).

2) Code used for performing tracking on videos taken on the glass surface (so in 2D). 



*Detailed outline:*
(written in chronological order)

In BashScripts folder:
- Run 'convertAviToImageSeq.sh' to convert avi to tif image sequence. Requires ffmpeg to be downloaded and be careful that the frame rate is correct.
- Run 'runTracking.sh' to perform all steps of tracking: findRods2Dt, trackRods2Dt and filterRods2Dt. 
- If you want to few trajectories of bacteria throughout video, run 'labelImages.py' and few image sequence using imageJ. This runs well but does not perform any compression so after all these steps, the files take up a lot of memory. You can use imageJ to compress the image sequence outputted by labelImages.py into a much smaller avi file though.

In PythonScripts folder:
(The names for a lot of these programs are really similar and extremely confusing. Sorry!)

Functions:
- 'analyseDDMTracks.py' contains a list of functions called by master programs. Despite the name, this was also called by surface analysis code.

For tracking on DDM videos:
- Run 'runDDMTrackingAllFiles.py' as master for performing DDM tracking.

For tracking on Surface videos:
- Run 'runSurfaceTrackingAllFiles.py' as master for surface tracking.
- Run 'analyseLysisTracks.py' to view individual tracked trajectories (use output of labelImages to match up trajectory ID to plotted trajectory outputted here.)
- Run 'readStopTimeFile.py' for reading text file created by 'analyseLysisTracks.py' and making Figure 19 in report.


*Other notes:*

To get .sh files to run, first run
	$chmod -x -a file.sh

To use imageJ on ubuntu, go to text notes on ubuntu and copy paste path from there into terminal (it would be better to make an alias and save in .bashrc file but i didn't get round to this). For some reason, it didn't work when i downlaoded it and that was the only way I figured out how to make it work.

There are a few differences between ubuntu and mac which requires small alterations to code in order to work:
1) font types are sometimes not recognised when printing on plots
2) Path of external hard drive is called 'Volumes' on Mac and 'Media' on ubuntu.
3) Ubuntu is much quicker than mac for running tracking so I would run on ubuntu and just ssh in commands if not at computer.





