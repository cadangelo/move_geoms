#include <iostream>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Types.hpp"
#include "MBTagConventions.hpp"
#include <string>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define CHECK_ERR(A) do { if (moab::MB_SUCCESS != (A)) {		     \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
  << __LINE__ << std::endl; \
  return A; } } while(false)

moab::Tag category_tag;
moab::Tag geom_tag;
moab::Tag name_tag;
moab::Tag obj_name_tag;
moab::Tag dim_tag, id_tag;

moab::ErrorCode get_all_handles(moab::Core *mbi)
{
  moab::ErrorCode rval;

  rval = mbi->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
			      name_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
  CHECK_ERR(rval);

  rval = mbi->tag_get_handle( "OBJECT_NAME", 32, moab::MB_TYPE_OPAQUE,
			      obj_name_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
  CHECK_ERR(rval);

  int negone = -1;
  rval = mbi->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
			      geom_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT,&negone);
  CHECK_ERR(rval);

  rval = mbi->tag_get_handle( GLOBAL_ID_TAG_NAME,
			      1, moab::MB_TYPE_INTEGER,
			      id_tag,
			      moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  CHECK_ERR(rval);
  
  rval = mbi->tag_get_handle( CATEGORY_TAG_NAME,
			      CATEGORY_TAG_SIZE,
			      moab::MB_TYPE_OPAQUE,
			      category_tag,
			      moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );

  CHECK_ERR(rval);
  return moab::MB_SUCCESS;
}
  
int main(int argc, char* argv[]) {

  moab::Core *mbi = new moab::Core();
  
  std::string output_file = "moved.h5m";

  get_all_handles(mbi);

  moab::ErrorCode rval;  
  rval = mbi->load_file("big_graveyard.h5m");
  
  moab::Range::iterator it;
  int id_num = 1;// ent_sets.size();
  int dim;
  moab::Range parents;
  moab::Range ent_sets;
  moab::Range vert_set;
  const double x=0.0, y=0.0, z=0.0;
  double xyz[3], xyz_trans[3];

  double total_dist = 5000.0;
  double dist;
  double dir_vec[3] = {1.0, 0.0, 0.0};
  double mag;
  double unit_vec[3];
  double speed = 1.0;
  double time; 
  double dist_per_shot;
  int num_shots = 2;
  int shot_num = 0;
  

  time = total_dist/speed;
  dist_per_shot = total_dist/num_shots;
  mag = sqrt(pow(dir_vec[0], 2.0) + pow(dir_vec[1], 2.0) + pow(dir_vec[2], 2.0));
  unit_vec[0] = dir_vec[0]/mag;
  unit_vec[1] = dir_vec[1]/mag;
  unit_vec[2] = dir_vec[2]/mag;
  std::cout << "unit " << unit_vec[0] << std::endl;
//  double trans_vec[3] = {1.0, 1.0, 1.0};

  char buffer [33];


  rval =  mbi->get_entities_by_type(0, moab::MBTRI, ent_sets);
  CHECK_ERR(rval);
  std::cout << "size " << ent_sets.size() << std::endl;
 
  rval = mbi->get_vertices(ent_sets, vert_set);
  
  dist = dist_per_shot;
  while (dist <= total_dist)
    {
      shot_num++;
      for (it = vert_set.begin(); it != vert_set.end(); ++it)
        {  
          rval = mbi->get_coords(&(*it), 1, xyz);
          xyz_trans[0] = xyz[0] + dist_per_shot*unit_vec[0];
          xyz_trans[1] = xyz[1] + dist_per_shot*unit_vec[1];
          xyz_trans[2] = xyz[2] + dist_per_shot*unit_vec[2];
          rval = mbi->set_coords(&(*it), 1, xyz_trans);
         
        }

      dist = dist + dist_per_shot;
//          rval = mbi->write_mesh((char)shot_num+output_file.c_str());
      rval = mbi->write_mesh( (std::to_string(shot_num)+output_file).c_str());

    }
  
/*  rval =  mbi->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &geom_tag, NULL, 1, ent_sets);
  CHECK_ERR(rval);
  for (it = ent_sets.begin(); it != ent_sets.end(); ++it) {
    rval = mbi->tag_get_data(geom_tag, &(*it), 1, &dim);
    CHECK_ERR(rval);
    if (dim == 2) {
      rval = mbi->tag_set_data(id_tag, &(*it), 1, &(id_num));
      CHECK_ERR(rval);
      parents.clear();
      rval = mbi->get_parent_meshsets(*it, parents);
      CHECK_ERR(rval);
      rval = mbi->tag_set_data(id_tag, &(*parents.begin()), 1, &(id_num));
      CHECK_ERR(rval);
      ++id_num;
    }
  }
  
  
  // save file
  rval = mbi->write_mesh(output_file.c_str());
 */ 

  return 0;
}    
