#/bin/bash

#MOAB_PATH=$HOME/opt/moab5.0/
MOAB_PATH=$HOME/research/moab/
MOAB_LIBRARY=$MOAB_PATH"lib"
MOAB_INCLUDE=$MOAB_PATH"include"
#DAGMC_PATH=$HOME/opt/DAGMC/
DAGMC_PATH=$HOME/research/DAGMC/
DAGMC_LIBRARY=$DAGMC_PATH"lib"
DAGMC_INCLUDE=$DAGMC_PATH"include"
BOOST_LIBRARY=/usr/lib/x86_64-linux-gnu/

echo $MOAB_LIBRARY
echo $DAGMC_LIBRARY

#g++  -std=c++11 move_geoms.cpp -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -ldagmc -lMOAB -o move
#g++  -std=c++11 mp.cpp -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -lMOAB -I$DAGMC_INCLUDE -L$DAGMC_LIBRARY -ldagmc -o mp
#g++  -std=c++11 translate_geom.cpp -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -lMOAB -I$DAGMC_INCLUDE -L$DAGMC_LIBRARY -ldagmc  -o trans
g++  -std=c++11 translate_geom.cpp -g -I$DAGMC_INCLUDE -L$DAGMC_LIBRARY -ldagmc -I$MOAB_INCLUDE -L$MOAB_LIBRARY -lMOAB -lboost_program_options -o trans

