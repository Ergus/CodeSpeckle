#!/bin/bash
#PBS -l walltime=10:00:00
#cd $PBS_O_WORKDIR

ulimit -s unlimited

module load intel mkl openmpi/1.8.3/intel/14.0 cudatoolkit/6.5 

export INSTDIR=/home/jaguilar/Install

export PLASMAPATH=$INSTDIR/Plasma/plasma-installer_2.7.1/install
export PLASMA_LIB=$PLASMAPATH/lib
export LD_LIBRARY_PATH=$PLASMA_LIB:$LD_LIBRARY_PATH

export MAGMAPATH=$INSTDIR/Magma/magma_1.7.0_install
export MAGMA_LIB=$MAGMAPATH/lib
export LD_LIBRARY_PATH=$MAGMA_LIB:$LD_LIBRARY_PATH

#mpirun --map-by ppr:1:node ./speckle_mpi.x --solver ${solver} --file ${solver}_ input
