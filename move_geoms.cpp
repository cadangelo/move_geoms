#include <iostream>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Types.hpp"
#include "MBTagConventions.hpp"

#include "DagMC.hpp"
using moab::DagMC;

#include <string>
#include <math.h>
#include <cmath>
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
moab::Tag move_tag;
moab::Tag obb_tag;
moab::Tag obb_tree_tag;

moab::DagMC *DAG;

moab::ErrorCode get_all_handles(moab::Core *mbi)
{
  moab::ErrorCode rval;

  rval = mbi->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
			      name_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
  CHECK_ERR(rval);

  rval = mbi->tag_get_handle( "OBJECT_NAME", 32, moab::MB_TYPE_OPAQUE,
			      obj_name_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
  CHECK_ERR(rval);

  rval = mbi->tag_get_handle( "MOVE_TAG", 32, moab::MB_TYPE_OPAQUE,
			      move_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
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
/*
  rval = mbi->tag_get_handle( MB_OBB_TREE_TAG_NAME, 1, MB_TYPE_HANDLE, 
                                obb_tree_tag, MB_TAG_DENSE );
  CHECK_ERR(rval);

  rval = mbi->tag_get_handle( MB_OBB_TAG_NAME, 1, MB_TYPE_HANDLE,
                                obb_tag, MB_TAG_DENSE );
  CHECK_ERR(rval);

*/
  return moab::MB_SUCCESS;
}
/*
void get_xyz_dists(double dx, double dy, double dz, double total_dist, int num_shots, double &dist_x, double &dist_y, double &dist_z)
{
  double mag;
  //double unit_vec[3];
  double dist_per_shot;

  mag = sqrt(pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0));
  ux[0] = dx/mag;
  uy[1] = dy[1]/mag;
  uz[2] = dz[2]/mag;

  dist_per_shot = total_dist/num_shots;
  dist_x = dist_per_shot*ux;
  dist_y = dist_per_shot*uy;
  dist_z = dist_per_shot*uz;
}

double* get_dist_uvec_td(double total_dist, double dist_per_shot, double unit_vec[3])
{
  static double dist[3];


  return dist;
}

 */ 
int main(int argc, char* argv[]) {

  moab::Core *mbi = new moab::Core();
  
  std::string output_file = "moved.h5m";

  get_all_handles(mbi);

  moab::ErrorCode rval;  
  rval = mbi->load_file("test_tag.h5m");

  DAG = new moab::DagMC(mbi);

  rval = DAG->load_existing_contents();
  CHECK_ERR(rval);

  rval = DAG->setup_obbs();
  CHECK_ERR(rval);
 
  rval = DAG->setup_indices();
  CHECK_ERR(rval);

  // get moving tag from parse_properties
  std::vector<std::string> group_name;
  std::map<std::string, std::string> group_name_synonyms;

  group_name.push_back("moving");

  rval = DAG->parse_properties(group_name, group_name_synonyms);
  if (moab::MB_SUCCESS != rval) 
    {
      std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
      exit(EXIT_FAILURE);
    }
  
  // get all volumes
  int num_cells = DAG->num_entities( 3 );
  std::cout << "num cells: " << num_cells << std::endl;
 
  // get all triangles in vols that are moving
  moab::EntityHandle moving;
  moab::Range tris;
  moab::Range vert_set;
  moab::Range move_sets;
  move_sets.clear();

  rval = mbi->create_meshset(moab::MESHSET_SET, moving);
  CHECK_ERR(rval);

  for( int i = 1; i <= num_cells; ++i ) 
    {  
      moab::EntityHandle vol = DAG->entity_by_index( 3, i );
      if( DAG->has_prop( vol, "moving"))
        {
          move_sets.insert(vol);

          mbi->add_entities(moving, &vol, 1);
        }
    }
  std::cout << "num move vols " << move_sets.size() << std::endl;

  int num_tris;
  rval =  mbi->get_number_entities_by_type(moving, moab::MBTRI, num_tris);
  CHECK_ERR(rval);
  std::cout<< "num dim 3" << num_tris << std::endl;

  rval =  mbi->get_entities_by_type(moving, moab::MBTRI, tris);
  CHECK_ERR(rval);
  std::cout << "num tris " << tris.size() << std::endl;

  rval = mbi->get_vertices(tris, vert_set);
  CHECK_ERR(rval);
  std::cout << "num moving verts: " << vert_set.size() << std::endl;

  // move verts according to transformation
  moab::Range::iterator it;
  const double x=0.0, y=0.0, z=0.0;
  double xyz[3], xyz_new[3];
 
//  double unit_vec[3];
  double dist[3];
  double v[3] = {100, 0, 0};
  double a[3] = {0, 0, 0};
//  double speed = 1.0;
  double t = 0.0; 
  double ts = 1.0; 
  double end_t = 20.0;
  int shot_num = 0;
//  time = total_dist/speed;
//  dist_per_shot = total_dist/num_shots;
//  double* unit_vec[3] = get_unit_vec(dir_vec);
//  double* dist[3] = get_dist_uvec_td(total_dist, dist_per_shot, unit_vec);
  
//  rval = mbi->tag_delete(obb_tag);
//  rval = mbi->tag_delete(obb_tree_tag);

  while (t < end_t)
    {
      shot_num++;
      t = t + ts;
      for (it = vert_set.begin(); it != vert_set.end(); ++it)
        {  
          rval = mbi->get_coords(&(*it), 1, xyz);
          xyz_new[0] = xyz[0] + v[0]*t + 0.5*a[0]*pow(t,2);
          xyz_new[1] = xyz[1] + v[1]*t + 0.5*a[1]*pow(t,2);
          xyz_new[2] = xyz[2] + v[2]*t + 0.5*a[2]*pow(t,2);
          rval = mbi->set_coords(&(*it), 1, xyz_new);
        
        }
      
      rval = mbi->write_mesh( (std::to_string(shot_num)+output_file).c_str());

    }
  
  delete DAG;
  return 0;
}    
