#!/bin/bash

# find first and last years
files='/g/data1/ua8/NOAA_OISST/AVHRR/v2-0_modified/*nc'
c=0
for f in $files
do
  filein=$f
  r=`echo $f | awk -F "_" '{print $(NF)}'`
  year[c]=`echo $r | awk -F "." '{print $(1)}'`
  c=$((c+1))
done
IFS=$'\n'
y2=`echo "${year[*]}" | sort -nr | head -n1`
y1=`echo "${year[*]}" | sort -nr | tail -n1`

# concaternat timeseries into regional chunks
echo ____ regional timeseries
for i in $(seq 0 30 350);
  do
  for j in $(seq 0 20 170);
  do
     fo=regional_data/avhrr-only-v2.$i\to$(($i+30)).$(($j-90))\to$(($j+20-90)).$y1.$y2.nc
     if [ ! -f $fo ]; then	
        # RUN AS PARALLEL JOBS
        echo $fo
    	JOB=`qsub  -q normal  -l walltime=0:20:00  -N reg_bloc_$i_$(($j-90)) - << EOJ
#PBS -S /bin/bash

cd /g/data1a/e14/asg561/MHW
module load nco/4.6.4
echo ncrcat -v sst -d lon,$((4*$i)),$((4*$i+4*30-1))  -d lat,$((4*$j)),$((4*$j+4*20-1))  /g/data1/ua8/NOAA_OISST/AVHRR/v2-0_modified/*nc $fo 
ncrcat -v sst  -d lon,$((4*$i)),$((4*$i+4*30-1))  -d lat,$((4*$j)),$((4*$j+4*20-1))  /g/data1/ua8/NOAA_OISST/AVHRR/v2-0_modified/*nc    $fo
EOJ
`
echo "JobID = ${JOB} for $fo submitted on `date`"

     fi
  done
done

exit
