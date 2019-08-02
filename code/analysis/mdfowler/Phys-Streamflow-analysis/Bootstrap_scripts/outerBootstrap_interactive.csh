#!/bin/csh

module load python3
module unload python

#export STARTLON=$((16+$STARTLON))

@ parallelcount = 0 
set xindex = 950
@ xmax = $xindex + 15
#set yindex = 0 

while ($xindex <= 1439)
   set yindex = 0
   while ($yindex <= 719)

       setenv XINDEX_CURRENT $xindex
       setenv YINDEX_CURRENT $yindex

       set fileBase="/work/04268/tg835671/stampede2/PythonBootstrap/bootstrapResults-NoParallel/testBoot1000_lonIndex"
       set fileMid="_latIndex"
       set fileEnd=".pkl"

       set fileName="$fileBase$xindex$fileMid$yindex$fileEnd"
       
       if ( -f $fileName ) then
           echo "$fileName already created. Skipping."
       else
           python3 GEVscript_ParallelSAMPLES_batchable.py & 

           @ parallelcount = $parallelcount + 1
           echo "Parallel count = $parallelcount" 

           if ($parallelcount == 911 ) then
               wait
               @ parallelcount = 0 
           endif
       endif
       @ yindex++
  end
 
  @ xindex++
 
  if ($parallelcount == 911 ) then
      wait
      @ parallelcount = 0
  endif

  echo "Longitude index: $xindex"
end

echo "Done with section of longitudes"

#set i=1
#@ parallelcount = 0
#while ($i <= $nBoot)
#   echo "Computing bootstrap number $i"
#   setenv iboot $i 
#   
#   python3 GEVscript_ParallelLON.py &
#
#   echo "Parallel count = $parallelcount" 
#   @ parallelcount = $parallelcount + 1
#   if ( $parallelcount == nBoot-1 ) then
#      wait
#      @ parallelcount = 0
#   endif

#  @ i++
#end


