#/bin/bash

MOAB_PATH=$HOME/opt/moab/
MOAB_LIBRARY=$MOAB_PATH"lib"
MOAB_INCLUDE=$MOAB_PATH"include"
DAGMC_PATH=$HOME/opt/DAGMC/
DAGMC_LIBRARY=$DAGMC_PATH"lib"
DAGMC_INCLUDE=$DAGMC_PATH"include"

echo $MOAB_LIBRARY
echo $DAGMC_LIBRARY

#g++  -std=c++11 move_geoms.cpp -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -ldagmc -lMOAB -o move
#g++  -std=c++11 mp.cpp -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -lMOAB -I$DAGMC_INCLUDE -L$DAGMC_LIBRARY -ldagmc -o mp
g++  -std=c++11 move_geoms.cpp -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -lMOAB -I$DAGMC_INCLUDE -L$DAGMC_LIBRARY -ldagmc -o move

