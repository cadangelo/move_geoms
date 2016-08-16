#/bin/bash

MOAB_PATH=$HOME/research/moab/
MOAB_LIBRARY=$MOAB_PATH"/lib"
MOAB_INCLUDE=$MOAB_PATH"/include"

echo $MOAB_LIBRARY

g++ move_geoms.cpp -g -std=c++11 -I$MOAB_INCLUDE -L$MOAB_LIBRARY -lMOAB -o move
