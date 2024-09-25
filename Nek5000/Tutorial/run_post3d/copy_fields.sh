#!/bin/bash

# copy history files:
INPUT_PATH='../run_tutorial'

cp $INPUT_PATH/duct.ma2 .
cp $INPUT_PATH/duct.re2 .

STAT_FOLDER=$(ls $INPUT_PATH | grep 'output_')

OUTPUT_PATH='STAT3D'

rm -rf $OUTPUT_PATH
mkdir $OUTPUT_PATH

i_total_snap_file=0

for i_folder in $STAT_FOLDER
do

	for i_local_stat in {2..6}
	do
        printf -v N_FILE "%05d" $i_local_stat
		TEST_FILE=$INPUT_PATH/$i_folder/s01duct0.f$N_FILE

		if [ -f "$TEST_FILE" ]; then
					
 			((i_total_snap_file++))
	  	        printf -v LAST "%05d" $i_total_snap_file

        		echo "copy s files from $INPUT_PATH/$i_folder"

        		for i_s in {1..11}
        		do
            		printf -v NS "%02d" $i_s
            		COPY_FROM0=$INPUT_PATH/$i_folder/s${NS}duct0.f$N_FILE
            		COPY_TO0=$OUTPUT_PATH/s${NS}duct0.f$LAST
            		echo "$COPY_FROM0 copied in $COPY_TO0"
            		cp $COPY_FROM0 $COPY_TO0
        		done

        		echo "copy t files from $INPUT_PATH/$i_folder"

        		for i_t in {1..4}
        		do
            		printf -v NT "%02d" $i_t
            		COPY_FROM0=$INPUT_PATH/$i_folder/t${NT}duct0.f$N_FILE
            		COPY_TO0=$OUTPUT_PATH/t${NT}duct0.f$LAST
            		echo "$COPY_FROM0 copied in $COPY_TO0"
            		cp $COPY_FROM0 $COPY_TO0
        		done
		else
			echo "$TEST_FILE does not exist."
		fi

	done

done
