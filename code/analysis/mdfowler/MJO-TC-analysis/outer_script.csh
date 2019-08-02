#!/bin/csh
# /gdata/pritchard/gkooperm/melissa_run/cesm_control/atm_h1/spinup.cam.h1.1873-07-05-00000.nc
set rundir="/fast/mdfowler/dailyGPIvars"
set resultsdir=`pwd`/results_OMI
mkdir -p $resultsdir

#set filelist=`ls $rundir/$case/atm_h1/*.cam.h1.*.nc`
#echo $filelist
@ parallelcount = 0
cd $rundir/
foreach localfil ( `ls DailyGPIvars_MJOphase*-OMI_iBoot100.nc` )
#foreach localfil ( `ls DailyGPIvars_climatology*.nc` )
  setenv GPIFILEIN $rundir/$localfil
  setenv GPIFILEOUT $resultsdir/GPI_$localfil
  #setenv LANDFRACFILE /lustre/SCRATCH/pritchard/mdfowler/GPI_daily/LANDFRAC.nc
  if ( ! -e $GPIFILEOUT ) then
    ncl < /DFS-B/DATA/pritchard/mdfowler/GPI_daily/bootstrap_ERAI/calculate_GPI.ncl & 
    echo $GPIFILEOUT
    echo $parallelcount
    @ parallelcount = $parallelcount + 1
    if ( $parallelcount == 8 ) then
      wait
      @ parallelcount = 0
    endif
  endif
end

