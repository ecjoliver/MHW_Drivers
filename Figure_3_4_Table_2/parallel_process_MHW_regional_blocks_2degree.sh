#!/bin/bash

# find regional files
N=`ls -l /g/data1a/e14/asg561/MHW/regional_data/*.nc | wc -l`
echo $N

if [ $N -ne 108 ]; then
  echo only $N of 108 files. Run parallel_process_data_to_regional_blocks.sh
  exit
fi
 
fn=`ls /g/data1a/e14/asg561/MHW/regional_data/*.nc| sort -n | head -1`
y1=`echo $fn | awk -F "." '{print $(NF-2)}'`
y2=`echo $fn | awk -F "." '{print $(NF-1)}'`

# process regional chunks
#for i in $(seq 0 30 350);
for i in $(seq 0 30 350);
  do
  #for j in $(seq 0 20 170);
  for j in $(seq 0 20 170);
  do
     fo=regional_data/avhrr-only-v2.$i\to$(($i+30)).$(($j-90))\to$(($j+20-90)).$y1.$y2.nc
     stem=`echo $fo | awk -F "v2." '{print $(NF)}'`
     fo_out=mhw_data_90pc_2degree/mhw_severity.pc90.$stem
     if [ -f $fo ]; then	
     if [ ! -f $fo_out ]; then
        # RUN AS PARALLEL JOBS
        echo processing $fo_out .....
    	JOB=`qsub  -q express -l walltime=1:00:00,mem=8000MB  -N m90_$i_$(($j-90)) - << EOJ
#PBS -S /bin/bash
cd /g/data1a/e14/asg561/MHW
module use ~access/modules
module load pythonlib/numpy/1.9.3
module load pythonlib/netCDF4/1.0.4
##python regional_MHW_pc90.py $fo
python regional_MHW_pc90_reducedFileSize_2degree.py $fo 

EOJ
`
echo "JobID = ${JOB} for $fo submitted on `date`"
     fi
     fi
  done
done

exit
