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

#define dot(u,v)   (u[0]*v[0] + u[1]*v[1] + u[2]*v[2])

moab::Tag category_tag;
moab::Tag geom_tag;
moab::Tag name_tag;
moab::Tag obj_name_tag;
moab::Tag dim_tag, id_tag;
moab::Tag move_tag;
moab::Tag obb_tag;
moab::Tag obb_tree_tag;

moab::DagMC *DAG;

struct XYZ{
  double x;
  double y;
  double z;
};

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

/*
double get_dist_to_axis(double xyz_P[3], double xyz_L0[3], double xyz_L1[3])
//double get_dist_to_axis(xyz_P, xyz_L0, xyz_L1)
{
  double v[3], w[3];
  v[0] = xyz_L1[0] - xyz_L0[0];
  v[1] = xyz_L1[1] - xyz_L0[1];
  v[2] = xyz_L1[2] - xyz_L0[2];

  w[0] = xyz_P[0] - xyz_L0[0];
  w[1] = xyz_P[1] - xyz_L0[1];
  w[2] = xyz_P[2] - xyz_L0[2];

  double c1, c2, b;
  c1 = dot(w,v);
  c2 = dot(v,v);
  b = c1/c2;
  
  double Pb[3];
  Pb[0] = xyz_L0[0] + b*v[0];
  Pb[1] = xyz_L0[1] + b*v[1];
  Pb[2] = xyz_L0[2] + b*v[2];

  //get distance from P to Pb
  double d[3], distance;
  d[0] = Pb[0] - xyz_P[0];   
  d[1] = Pb[1] - xyz_P[1];   
  d[2] = Pb[2] - xyz_P[2];   

  return sqrt(dot(d,d));

}
*/

/* Funtion to rotate a 3D point a distance theta 
   around any line (given by two points)
*/
XYZ rotate_point(XYZ P, double theta, XYZ L1, XYZ L2)
{
   XYZ v, u, q1, q2;
   double m, d;

   //translate so that rotation axis is origin
   q1.x = P.x - L1.x;
   q1.y = P.y - L1.y;
   q1.z = P.z - L1.z;

   // find unit vector of rotation axis
   v.x = L2.x - L1.x;
   v.y = L2.y - L1.y;
   v.z = L2.z - L1.z;
   m = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
   u.x = v.x/m;
   u.y = v.y/m;
   u.z = v.z/m;

   // rotate space about x axis so v lies in xz plane

   // length of u projected onto yz plane
   d = sqrt(u.y*u.y + u.z*u.z);

   //if d = 0, v is already in xz plane
   if (d != 0) 
     {
       q2.x = q1.x;
       q2.y = q1.y * u.z / d - q1.z * u.y / d;
       q2.z = q1.y * u.y / d + q1.z * u.z / d;
     } 
   else 
     {
       q2 = q1;
     }

   // rotate space about y axis so v lies along z axis
   q1.x = q2.x * d - q2.z * u.x;
   q1.y = q2.y;
   q1.z = q2.x * u.x + q2.z * d;
 
   // rotate space by angle theta about z axis
   q2.x = q1.x * cos(theta) - q1.y * sin(theta);
   q2.y = q1.x * sin(theta) + q1.y * cos(theta);
   q2.z = q1.z;

   // inverse of y axis rotation
   q1.x =   q2.x * d + q2.z * u.x;
   q1.y =   q2.y;
   q1.z = - q2.x * u.x + q2.z * d;

   // inverse of x axis rotation
   if (d != 0) 
     {
       q2.x =   q1.x;
       q2.y =   q1.y * u.z / d + q1.z * u.y / d;
       q2.z = - q1.y * u.y / d + q1.z * u.z / d;
     } 
   else
     {
       q2 = q1;
     }

   // inverse of translation to origin
   q1.x = q2.x + L1.x;
   q1.y = q2.y + L1.y;
   q1.z = q2.z + L1.z;

   // return rotated point
   return(q1);
}

moab::Range get_tagged_entities(moab::Core *mbi, int total_cells, std::string tag_name)
{
  moab::ErrorCode rval;

  // get moving tag from parse_properties
  std::vector<std::string> group_name;
  std::map<std::string, std::string> group_name_synonyms;

  group_name.push_back(tag_name);

  rval = DAG->parse_properties(group_name, group_name_synonyms);
  if (moab::MB_SUCCESS != rval) 
    {
      std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
      exit(EXIT_FAILURE);
    }

  // get all triangles in vols that have tag
  moab::EntityHandle tagged_meshset;
  moab::Range vert_set;
  moab::Range surf_set;
  moab::Range tagged_verts;
  int num_tris, num_verts;
  moab::Range::iterator it, itr;

  surf_set.clear();
  vert_set.clear();

  rval = mbi->create_meshset(moab::MESHSET_SET, tagged_meshset);
  //CHECK_ERR(rval);

  for( int i = 1; i <= total_cells; ++i ) 
    {  
      moab::EntityHandle vol = DAG->entity_by_index( 3, i );
      if( DAG->has_prop( vol, tag_name))
        {
          mbi->get_child_meshsets(vol, surf_set);
          for (it = surf_set.begin(); it != surf_set.end(); it++)
           {
             rval =  mbi->get_entities_by_type(*it, moab::MBVERTEX, vert_set);
             for (itr = vert_set.begin(); itr != vert_set.end(); itr++)
               {
                 tagged_verts.insert(*itr);
               }
           }
        }
          
    }

  return tagged_verts;
}

int main(int argc, char* argv[]) 
{
  moab::ErrorCode rval; 
 
  moab::Core *mbi = new moab::Core();
 
  // get all moab tag handles 
  get_all_handles(mbi);

  // load base geometry file that we wish to move
  rval = mbi->load_file("pie.h5m");

  DAG = new moab::DagMC(mbi);

  rval = DAG->load_existing_contents();
  CHECK_ERR(rval);

  //  rval = DAG->setup_obbs();
  //  CHECK_ERR(rval);
 
  rval = DAG->setup_indices();
  //  CHECK_ERR(rval);


  
  // get all volumes
  int num_cells = DAG->num_entities( 3 );
  std::cout << "num cells: " << num_cells << std::endl;
 
  moab::Range mv;
  mv = get_tagged_entities(mbi, num_cells, "moving");

  std::cout << "num moving verts: " << mv.size() << std::endl;
  //std::cout << "num grave verts: " << gv.size() << std::endl;

  //Inital point, updated point
  XYZ p_0, p, p_new;
  double xyz[3], xyz_new[3];

  int transform; // 0 for tranlation, 1 for rotation

  // TRANSLATION
  //transform = 0;

  // initial velocity vector for 3D translation
  XYZ v, v_0;
  v_0.x = 100;
  v_0.y = 0;
  v_0.z = 0;

  // acceleration vector for 3D translation
  // a(t) = b + ct + dt^2
  XYZ b, c, d;
  b.x = 0;
  c.x = 0;
  d.x = 0;
  b.y = 0;
  c.y = 0;
  d.y = 0;
  b.z = 0;
  c.z = 0;
  d.z = 0;


  // ROTATION
  transform = 1;
 
  // alpha(t) = alpha[0]+ alpha[1]*t + alpha[2]t^2
  double alpha[3] = {0, 0, 0}; //angular acceleration [rad/s^2]
  double omega, omega_0 = M_PI/2; //angular velocity [rad/s]
  double theta; //displacement [rad]

  // Two points that define line points rotate about 
  XYZ L0, L1;
  L0.x = 0;
  L0.y = 250;
  L0.z = 0;
  L1.x = 0;
  L1.y = 250;
  L1.z = 1;

  double t = 0.0; //current time [s]
  double ts = 0.25; //length of time step [s]
  double end_t = 24.0; //end time [s]
  int shot_num = 0; //current time step


  std::map<moab::EntityHandle, XYZ> position;
 
  //base output file name 
  std::string output_file = "moved.h5m";

  moab::Range::iterator its;

  while (t <= end_t)
    {
      for (its = mv.begin(); its != mv.end(); ++its)
        {
          if(t == 0)
            {
              //get starting position
              rval = mbi->get_coords(&(*its), 1, xyz);
              CHECK_ERR(rval);
           
              p.x = xyz[0];
              p.y = xyz[1];
              p.z = xyz[2];

              //map original position
              position[*its] = p;
            }

          else
            {
              // get original position
              p_0 = position.find(*its)->second;
 
              //if translation
              if (transform == 0)
                {
                   p_new.x = p_0.x + v_0.x*t + (1/2)*b.x*pow(t,2) + (1/6)*c.x*pow(t,3) + (1/12)*d.x*pow(t,4);
                   p_new.y = p_0.y + v_0.y*t + (1/2)*b.y*pow(t,2) + (1/6)*c.y*pow(t,3) + (1/12)*d.y*pow(t,4);
                   p_new.z = p_0.z + v_0.z*t + (1/2)*b.z*pow(t,2) + (1/6)*c.z*pow(t,3) + (1/12)*d.z*pow(t,4);
                }
           
              //if rotation
              if (transform == 1)
                { 
                  omega = omega_0 + alpha[0]*t + (1/2)*alpha[1]*pow(t,2) + (1/3)*alpha[2]*pow(t,3);
                  theta = omega*t;
                  p_new = rotate_point(p_0, theta, L0, L1);
                }

              xyz_new[0] = p_new.x;
              xyz_new[1] = p_new.y;
              xyz_new[2] = p_new.z;
              rval = mbi->set_coords(&(*its), 1, xyz_new);
              CHECK_ERR(rval);
            }
        }

      rval = mbi->write_mesh( (std::to_string(shot_num)+output_file).c_str());
      shot_num++;
      t = t + ts;
    }
  
  delete DAG;
  return 0;
}    
