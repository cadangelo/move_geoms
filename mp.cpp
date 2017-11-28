#include <iostream>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Types.hpp"
#include "MBTagConventions.hpp"

#include <string>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#define CHECK_ERR(A) do { if (moab::MB_SUCCESS != (A)) {		     \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
  << __LINE__ << std::endl; \
  return A; } } while(false)

struct XYZ{
  double x;
  double y;
  double z;
};


//moab::ErrorCode get_verts(moab::Core *mbi, moab::EntityHandle vol, moab::Range &verts)
moab::ErrorCode get_verts(moab::Core *mbi, moab::Range &verts)
{
  moab::ErrorCode rval;
  moab::Range vert_set;
  int num_verts;
  moab::Range::iterator itr;

  vert_set.clear();
  rval =  mbi->get_entities_by_type(0, moab::MBVERTEX, vert_set);
  CHECK_ERR(rval);
  for (itr = vert_set.begin(); itr != vert_set.end(); itr++)
    {
      verts.insert(*itr);
    }

  return moab::MB_SUCCESS;
}

int main(int argc, char* argv[]) 
{
  moab::ErrorCode rval; 
 
  moab::Core *mbi = new moab::Core();

  // load base geometry file that we wish to move
  char* filename = argv[1];
  rval = mbi->load_file(filename);
  CHECK_ERR(rval);

  //base output file name 
  std::string output_file = "moved.h5m";

  moab::Range mv;
  rval = get_verts(mbi, mv);
  CHECK_ERR(rval);

  //Inital point, updated point
  XYZ p_0, p, p_new;
  double xyz[3], xyz_new[3];
  double dx = 30;
  double dy = 0; 
  double dz = 81;

  moab::Range::iterator itt;
  for (itt = mv.begin(); itt != mv.end(); ++itt)
    {
     
   
          //get original position
          rval = mbi->get_coords(&(*itt), 1, xyz);
          CHECK_ERR(rval);
          
          p_0.x = xyz[0];
          p_0.y = xyz[1];
          p_0.z = xyz[2];
   
          //if translation
          p_new.x = p_0.x + dx; // v_0.x*t + (1/2)*b.x*pow(t,2) + (1/6)*c.x*pow(t,3) + (1/12)*d.x*pow(t,4);
          p_new.y = p_0.y + dy; // v_0.y*t + (1/2)*b.y*pow(t,2) + (1/6)*c.y*pow(t,3) + (1/12)*d.y*pow(t,4);
          p_new.z = p_0.z + dz; // v_0.z*t + (1/2)*b.z*pow(t,2) + (1/6)*c.z*pow(t,3) + (1/12)*d.z*pow(t,4);

          //set the coordinates of the updated position  
          xyz_new[0] = p_new.x;
          xyz_new[1] = p_new.y;
          xyz_new[2] = p_new.z;
          rval = mbi->set_coords(&(*itt), 1, xyz_new);
          CHECK_ERR(rval);

    }   

 rval = mbi->write_mesh( ("phan_"+output_file).c_str() );
 CHECK_ERR(rval);

 return 0; 
}
