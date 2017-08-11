
#652686-652796

EXPNUM=652686

ntot=110

n=1
startTime=`date`
echo "$startTime"
rm *.list
while [ $n -le $ntot ]
do  
    echo "Submitting exposure $EXPNUM ($n of $ntot)"
    ./DAGMaker.sh $EXPNUM > logs/DAGMaker_$EXPNUM.log 2>&1
    jobsub_submit_dag -G des --role=DESGW file://desgw_pipeline_${EXPNUM}.dag > logs/jobsub_submit_dag_$EXPNUM.log 2>&1
    n=`expr $n + 1`
    EXPNUM=`expr $EXPNUM + 1`
done
endTime=`date` 
echo "End of bulksubmit.sh script"
echo "start time: $startTime"
echo "end time: $endTime"
echo "number of exposures: $ntot"
echo "last EXPNUM: $EXPNUM"
