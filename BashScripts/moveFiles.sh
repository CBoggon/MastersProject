#!/bin/bash


## Note on code: 
## When running on external hard drive, this code seemed to copy over a few wrong files (for 2500 files, roughly 15 of the wrong files had been copied over). I think this was due to the hard drive and not the code but should check this nonetheless!

## WARNING: I didn't use this code, so I wouldn't trust it. 

numFilesToMove=2500

### TO MOVE FILES FROM ONE FOLDER TO ANOTHER FOLDER ONLY ONCE:
#fileDir=Pos00_Movie0000-ImSeq-1/
#newFileDir=Pos00_Movie0000-ImSeq-2/
#find $fileDir -maxdepth 1 -type f |head -$numFilesToMove |xargs cp -t $newFileDir


### TO REPEAT ACTION FOR ALL FOLDERS IN DIRECTORY
for file in Pos*-2
do
	
	fileDir=$file/
	#echo "$fileDir"
	newFileDir=${file%-2*}  # retain the part before the 1
	newFileDir=$newFileDir'-1/'
	#echo "$newFileDir"
	echo " "
	echo " "
	echo "### Moving the first $numFilesToMove files from $fileDir to $newFileDir ###"
	find $fileDir -maxdepth 1 -type f |head -$numFilesToMove |xargs mv -t $newFileDir
	
	#mv $fileDir/* $newFileDir/

done

echo " "
echo " "

