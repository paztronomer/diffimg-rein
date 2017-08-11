#!/bin/bash

if [ "$1" != "" ]; then
    readarray -t a < $1
else
    echo "Provide the file with the list of exposures you want to loop through in the command line:"
    echo
    echo "./loopsubmit.sh /path/to/file"
    echo
    exit 1
fi

max=${a[0]}
min=${a[0]}

for i in "${a[@]}"
do

    if [[ "$i" -gt "$max" ]]; then
        max="$i"
    fi

    if [[ "$i" -lt "$min" ]]; then
        min="$i"
    fi
done

echo

startTime=`date`
echo "$startTime"
rm *.list
echo "-----"
echo

count=0

while read e
do
    count=`expr $count + 1`
    echo "Submitting exposure $EXPNUM ($count of ${#a[@]})"
    ./DAGMaker.sh $e > logs/DAGMaker_$e.log 2>&1
    jobsub_submit_dag -G des --role=DESGW file://desgw_pipeline_${e}.dag > logs/jobsub_submit_dag_$e.log 2>&1
done < $1

endTime=`date`
echo "End of bulksubmit.sh script"
echo "start time: $startTime"
echo "end time: $endTime"
echo "number of exposures: $count"
echo "first expnum: $min"
echo "last expnum: $max"
