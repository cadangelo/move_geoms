#/bin/bash

MOAB_PATH=$HOME/opt_n/moab/
MOAB_LIBRARY=$MOAB_PATH"lib"
MOAB_INCLUDE=$MOAB_PATH"include"
DAGMC_PATH=$HOME/opt_n/dagmc/
DAGMC_LIBRARY=$DAGMC_PATH"lib"

echo $MOAB_LIBRARY
echo $DAGMC_LIBRARY

g++ move_geoms.cpp -std=c++11 -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -ldagmc -lMOAB -o move

