#!/bin/bash

# Program loops through all image sequence folders within specified file path, creating output folders containing the outputs of findRods2Dt and trackRods2Dt and filterRods.

### Input File Directories:
findRodsDir=../../teunTrackingCode/findRods2Dt/src
trackRodsDir=../../teunTrackingCode/trackRods2Dt/src
filterTracksDir=../../teunTrackingCode/filterTracks2Dt/src

#ImageDir="TestBash-ImConvert_Folder/DDMmovies171215-164612-AsImageSequences"
#ImageDir=../../../../../Volumes/My\ Book/MastersProject/Data/20171215/171215DDMx20Mic150fps-Attempt3/DDMmovies171215-164612/
##ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20171201/20171201Mic1DDM50fps-PhaceConc1e9-Attempt2/DDMmovies171201-154009/DDMmovies171201-154009-AsImageSequences
#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20171201/20171201Mic1DDM50fps-PhaceConc1e9-Attempt2/DDMmovies171201-154521/DDMmovies171201-154521-AsImageSequences
##ImageDir=../../../../../../media/cameron/MyBook/MastersProject/Data/20171201/20171201Mic1DDM50fps-PhaceConc1e9-Attempt2/DDMmovies171201-154345/DDMmovies171201-154345-AsImageSequences
#ImageDir=../../../../../../media/cameron/MyBook/MastersProject/Data/20171201/20171201Mic1DDM50fps-PhaceConc1e9/DDMmovies171201-142252-AsImageSequences
#ImageDir=../../../../../../media/cameron/MyBook/MastersProject/Data/20171024/x40Mic1Phage10fpsVid4+DDM-AsImageSequences

#ImageDir=../../../../../../media/cameron/MyBook/MastersProject/Data/20180206/20180206-50fpsx20-SurfacePhage+ControlMeas3/DDMmovies180206-164041-AsImageSequences
#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180206/20180206-50fpsx20-SurfacePhage+ControlMeas3/DDMmovies180206-164351-AsImageSequences

#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences

### 20180213-Surface:
#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-150927-AsImageSequences
#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-151805-AsImageSequences
ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-152022-AsImageSequences
#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-154132-AsImageSequences

### 20180202-DDM
#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326-AsImageSequences
#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135725-AsImageSequences
#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140000-AsImageSequences
#ImageDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140137-AsImageSequences



echo 'Editing These files:'
echo ' '
ls $ImageDir/Pos*
echo ' '
echo ' '
echo 'which are found in directory: '$ImageDir

# Loop through all image sequence folders
for folder in $ImageDir/Pos*
do
        echo " "        # New line on terminal
        echo "$folder"

        NAME=${folder%-ImSeq*}  # retain the part before the '-ImSeq'
        NAME=${NAME##*/}  # retain the part after the last slash
        echo $NAME
	
	outputFolder="$ImageDir/Output-$NAME"
	mkdir $outputFolder
	findRodsOutputFolder=$outputFolder/findRods2DtOutput
	mkdir $findRodsOutputFolder
	trackRodsOutputFolder=$outputFolder/trackRods2DtOutput
	mkdir $trackRodsOutputFolder
	filterTracksOutputFolder=$outputFolder/filterTracks2DtOutput
	mkdir $filterTracksOutputFolder
	
	

	#### Run findRods2Dt  ######################################################################
	echo ' '
	echo ' '
	echo " #### Running findRods2Dt ######"
	inFile=$findRodsDir/findRods2Dt.in
	echo "$inFile"
        echo "FILEMASK $folder/Pos*.tif" > $inFile
        echo "FILEMASKBG $ImageDir/Background-$NAME/*.tif" >> $inFile
        echo "POSITIONSFILE $findRodsOutputFolder/positions.dat" >> $inFile
        echo "OUTPUTDIR $findRodsOutputFolder/" >> $inFile

#	# x20 Magnification DDM
#        echo "FINDRODS 1" >> $inFile
#        echo "MINVAL 0.10" >> $inFile
#        echo "DIAMETER 3.0" >> $inFile
#        echo "BLURDIAMETER 0.7" >> $inFile
#        echo "SUBBACKGROUND 0.2" >> $inFile
#        echo "BGTHRESHOLD 0.80" >> $inFile
#        echo "DEBUG 1" >> $inFile
#        echo "OVERLAPR 1.0" >> $inFile
#        echo "TOPHATDIAMETER 15" >> $inFile
#        echo "MININTENS 10" >> $inFile
#        echo "WRITEFREQ 1" >> $inFile
#        echo "TRACERR 2" >> $inFile
#        echo "WHITE 0" >> $inFile
#        echo "RGB 0" >> $inFile
#        echo "CHUNK 2000" >> $inFile

	# x20 Magnification at Glass Surface.
        echo "FINDRODS 1" >> $inFile
        echo "MINVAL 0.05" >> $inFile
        echo "DIAMETER 2.8" >> $inFile
        echo "BLURDIAMETER 0.8" >> $inFile
        echo "SUBBACKGROUND 0.2" >> $inFile
        echo "BGTHRESHOLD 0.95" >> $inFile
        echo "DEBUG 1" >> $inFile
        echo "OVERLAPR 0.0" >> $inFile
        echo "TOPHATDIAMETER 15" >> $inFile
        echo "MININTENS 10" >> $inFile
        echo "WRITEFREQ 1" >> $inFile
        echo "TRACERR 2" >> $inFile
        echo "WHITE 0" >> $inFile
        echo "RGB 0" >> $inFile
        echo "CHUNK 2000" >> $inFile
	
#	# x40 Magnification. Note: CHUNK == 5000 !
#        echo "FINDRODS 1" >> $inFile
#        echo "MINVAL 0.20" >> $inFile
#        echo "DIAMETER 20.0" >> $inFile
#        echo "BLURDIAMETER 0.8" >> $inFile
#        echo "SUBBACKGROUND 0.0000001" >> $inFile
#        echo "BGTHRESHOLD 0.70" >> $inFile
#        echo "DEBUG 1" >> $inFile
#        echo "OVERLAPR 1.0" >> $inFile
#        echo "TOPHATDIAMETER 15" >> $inFile
#        echo "MININTENS 5" >> $inFile
#        echo "WRITEFREQ 1" >> $inFile
#        echo "TRACERR 2" >> $inFile
#        echo "WHITE 0" >> $inFile
#        echo "RGB 0" >> $inFile
#        echo "CHUNK 1000" >> $inFile

	runFindRodsFile=${inFile%.*}
	./$runFindRodsFile
				



	####  Run trackRods2Dt   #############################################################
	echo ' '
	echo ' '
	echo " #### Running trackRods2Dt ######"
	inTrackFile=$trackRodsDir/trackRods2Dt.in
	echo "$inTrackFile"

	echo "INPUTFILE $findRodsOutputFolder/positions.dat" > $inTrackFile
	echo "OUTPUTFILE $trackRodsOutputFolder/tracks.dat" >> $inTrackFile
	echo "DEBUG 1" >> ./$inTrackFile
	echo "STARTFRAME 1" >> $inTrackFile
	
	runTrackRodsFile=${inTrackFile%.*}
	./$runTrackRodsFile





	####   Run filterTracks2Dt   ################################################################
	echo ' '
	echo ' '
	echo " #### Running filterTracks2Dt ######"
	inFilterFile=$filterTracksDir/filterTracks2Dt.in
	echo "$inFilterFile"

	echo "TRACKINPUTFILE $trackRodsOutputFolder/tracks.dat" > $inFilterFile
	echo "OUTPUTDIR $filterTracksOutputFolder/" >> $inFilterFile
	echo "DEBUG 1" >> $inFilterFile
	echo "TT 10.0" >> $inFilterFile
	## NB: These values are currently being input manually in filterTracks2Dt.c line 50 and not being read in from here.
	
	runFilterTracksFile=${inFilterFile%.*}
	./$runFilterTracksFile

done


