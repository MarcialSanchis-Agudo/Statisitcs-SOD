#!/bin/bash
export CASE="duct"
export SOURCE_ROOT="/cfs/klemming/scratch/w/walesiak/Nek5000/Code-dardel"
export FC="ftn" #Dardel uses different compiler ftn
export CC="cc" #cc
export CFLAGS="-std=legacy -mcmodel=large"
export FFLAGS="-std=legacy -fallow-argument-mismatch -O2 -mcmodel=large"
#export CFLAGS="-std=legacy -mcmodel=large"
#export FFLAGS="-mcmodel=large"

export PPLIST=""
export USR="frame.o frame_usr.o io_tools_block.o io_tools.o mntrlog_block.o mntrlog.o mntrtmr_block.o mntrtmr.o rprm_block.o rprm.o map2D.o stat.o stat_IO.o chkpoint.o chkpt_mstp.o trip.o tsrs.o tsrs_IO.o pts_redistribute.o"

for il in "$@"
do
case $il in
	--clean)
		${SOURCE_ROOT}/bin/makenek clean
		shift
		;;
	--compile)
                ${SOURCE_ROOT}/bin/makenek ${CASE}
		shift
		;;
esac
done

