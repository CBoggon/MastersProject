#!/bin/bash

# Script takes folder of avi movies, creates a new folder and fills with these movies as tif image sequences.
# It then creates a folder of background images, randomly selected from the image folder.


### Declare Variables to be used in scipt
MAXNUMBERBGIMAGES=50	# Max number of background images
count=0			# Count for while loop
#fps=50			# frames per second
fps=50			# frames per second

#folderDir=TestBashScriptFolder/DDMmovies171215-164612
#folderDir=../../../../../../media/cameron/My\ Book/MastersProject/Data/20171215/171215DDMx20Mic150fps-Attempt3/DDMmovies171215-164921
#folderDir=../../../../../../media/cameron/MyBook/MastersProject/Data/20171201/20171201Mic1DDM50fps-PhaceConc1e9-Attempt2/DDMmovies171201-154009
#folderDir=../../../../../../media/cameron/MyBook/MastersProject/Data/20171201/20171201Mic1DDM50fps-PhaceConc1e9-Attempt2/DDMmovies171201-154521
#folderDir=../../../../../../media/cameron/MyBook/MastersProject/Data/20171201/20171201Mic1DDM50fps-PhaceConc1e9-Attempt2/DDMmovies171201-154345
#folderDir=../../../../../../media/cameron/MyBook/MastersProject/Data/20171201/20171201Mic1DDM50fps-PhaceConc1e9/DDMmovies171201-142252
#folderDir=../../../../../../media/cameron/MyBook/MastersProject/Data/20171024/x40Mic1Phage10fpsVid4+DDM
#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326

#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180206/20180206-50fpsx20-SurfacePhage+ControlMeas3/DDMmovies180206-164041
#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180206/20180206-50fpsx20-SurfacePhage+ControlMeas3/DDMmovies180206-164351

## 20180213-Surface:
#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-150927
#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-151805
#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-152022
folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180213/20180213Surface2Samples-50fpsx20Mag/DDMmovies180213-154132


## 20180202-DDM:
#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135326
#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-135725
#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140000
#folderDir=../../../../../media/cameron/MyBook/MastersProject/Data/20180202/20170202DDMx20-50fps/DDMmovies180202-140137


#REMEMBER TO REMOVE BACKSLASH FROM END OF DIRECTORY!!!


echo 'Editing These files:'
ls $folderDir
echo 'which are found in directory: '$folderDir


### Creating new foler folder to output images
newFolderDir=$folderDir-AsImageSequences
mkdir $newFolderDir
###mkdir TestBashScriptFolder/DDMmovies171215-164612/DDMmovies171215-164612-AsImageSequences

### Convert avi images to tif image sequence and save in new directory.
#for number in {1..10}
for file in $folderDir/Pos*.avi
do
	echo " "	# New line on terminal
	echo "$file"
	
	NAME=${file%.*}  # retain the part before the fullstop
	NAME=${NAME##*/}  # retain the part after the last slash
	echo $NAME			#NN
	imageSeqFolder=$newFolderDir/$NAME-ImSeq
	echo $imageSeqFolder	#NN
	mkdir $imageSeqFolder
	
	imType="_%04d.tif"
	imageName="$imageSeqFolder/$NAME$imType"
	echo $imageName		#NN. Name of individual files to be outputted by fmp
	#fmp command goes here. Input folder = $file, Output folder = $imageName
	echo " "	# New line on terminal
	#echo "ffmpeg -i $file -r 10 $imageName"
	ffmpeg -i $file -r $fps $imageName

	###Create Background Images
	backgroundFolder=$newFolderDir/Background-$NAME
	echo "$backgroundFolder"
	mkdir $backgroundFolder
	RANGE=$(ls -l $imageSeqFolder | grep -v ^d | grep -v ^t | wc -l)
	#echo "$RANGE"
	while [ "$count" -le $MAXNUMBERBGIMAGES ]      # Generate 10 ($MAXNUMBERBGIMAGES) random integers.
	do
		#select random image and copy to background folder
		#number=$RANDOM
		#let "number %= $RANGE"
		#echo $number
		number=$(shuf -i 1-$RANGE -n 1)		#get random number between 1 and RANGE
		printf -v number "%04d" $number
		echo $number
		randIm="_$number.tif"
		echo "$imageSeqFolder/$NAME$randIm"	
		cp $imageSeqFolder/$NAME$randIm $backgroundFolder/
		let "count += 1"  # Increment count.
	done
	count=0
done



exit 0


