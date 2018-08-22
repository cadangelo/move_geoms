#include <iostream>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Types.hpp"
#include "moab/GeomTopoTool.hpp"
#include "MBTagConventions.hpp"
#include "DagMC.hpp"
#include "dagmcmetadata.hpp"
#include <bits/stdc++.h>

#include <string>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

moab::DagMC *DAG;
moab::Core *mbi = new moab::Core;

struct XYZ{
  double x;
  double y;
  double z;
};

moab::ErrorCode setup(moab::Core *mbi, char* filename)
{
  moab::ErrorCode rval;

  DAG = new moab::DagMC(mbi);
  moab::GeomTopoTool *GTT = new moab::GeomTopoTool(mbi);

  // load base geometry file that we wish to move
  rval = mbi->load_file(filename);
  MB_CHK_ERR(rval);

  rval = DAG->load_existing_contents();
  MB_CHK_ERR(rval);

  // find all geometry sets
  rval = GTT->find_geomsets();
  MB_CHK_SET_ERR(rval, "GeomTopoTool could not find the geometry sets");

  // setup indices
  rval = DAG->setup_indices();
  MB_CHK_SET_ERR(rval, "Failed to setup problem indices");


}

void tokenize( const std::string& str, 
               std::vector<std::string>& tokens,
               const char* delimiters)
{
  tokens.clear();

  std::string::size_type next_token_end, next_token_start =
                         str.find_first_not_of( delimiters, 0);

  while ( std::string::npos != next_token_start )
    {
      next_token_end = str.find_first_of( delimiters, next_token_start );
      if ( std::string::npos == next_token_end )
        {
	  tokens.push_back(str.substr(next_token_start));
          next_token_start = std::string::npos;
        }
      else
        {
          tokens.push_back( str.substr( next_token_start, next_token_end -
                                        next_token_start ) );
          next_token_start = str.find_first_not_of( delimiters, next_token_end );
        }
    }
}

moab::ErrorCode get_tagged_vols(moab::Core *mbi,
                                std::map<int, moab::Range> &tagged_vols_map)
{
  moab::ErrorCode rval;

  // Get vols tagged w/ transformation
  std::string tag_name = "tr";
  std::string tag_delims = ":";

  // get tr tag from parse_properties
  std::vector<std::string> group_name;
  std::map<std::string, std::string> group_name_synonyms;
  group_name.push_back(tag_name);
  rval = DAG->parse_properties(group_name, group_name_synonyms, tag_delims.c_str());
  MB_CHK_SET_ERR(rval, "DAGMC failed to parse metadata properties");

  // get total number of volumes
  int all_vols = DAG->num_entities( 3 );

  // for each DAG volume
  for( int i = 1; i <= all_vols; ++i ){ 
   //get DAG index
   moab::EntityHandle vol = DAG->entity_by_index( 3, i );
   // if vol tagged w/ TR number
   if( DAG->has_prop( vol, tag_name)){
     std::vector<std::string> tr_nums;
     rval = DAG->prop_values(vol, tag_name, tr_nums);
     // get the value of each property (i.e. tag number)-- each vol can have more than one TR 
     for (int j = 0 ; j < tr_nums.size() ; j++) {
       bool added = false;
       int val = std::stoi(tr_nums[j]);
       int map_val;
       //if the map is not empty, look for a key that matches the tr # of the volume
       if(tagged_vols_map.size() > 0){
         std::map<int, moab::Range>::iterator itt;
         for(itt = tagged_vols_map.begin(); itt != tagged_vols_map.end(); ++itt){
            map_val = itt->first;
            //if val matches one already in map (map_val), add vol to range 
            if(val == map_val){
              tagged_vols_map[map_val].insert(vol); 
              added = true;
              break;
            }
         }//for each key in map
       }//if tr map not empty
       //if map empty or no matching key found, create key and add vol to range
       if(added == false) {
         map_val = val;
         tagged_vols_map[map_val].insert(vol); 
         added = true;
       }//if empty, create new key
     }// for each prop
   }//has prop
 }// for each DAG vol

  return moab::MB_SUCCESS;
}

moab::ErrorCode get_orig_vert_position(moab::Range mv , std::map<moab::EntityHandle, XYZ> &orig_positions)
{

  moab::ErrorCode rval;
  moab::Range::iterator itt;
  double xyz[3];
  XYZ p;
  //get starting position of each moving vertex
  for (itt = mv.begin(); itt != mv.end(); ++itt)
    {
      rval = mbi->get_coords(&(*itt), 1, xyz);
      MB_CHK_ERR(rval);
      
      p.x = xyz[0];
      p.y = xyz[1];
      p.z = xyz[2];
      
      //keep original position of each vertex
      orig_positions[*itt] = p;
    }
}

moab::ErrorCode get_tagged_verts(std::map<int, moab::Range> tagged_vols_map,
                                 moab::Range &tagged_verts,
                                 std::map<int, moab::Range> &tagged_verts_map)
{
  moab::ErrorCode rval;
  int tr_num;
  moab::Range surf_set, vert_set, edge_set, v_set;
  moab::Range::iterator its, itv, itr, ite, itvt; 
  std::map<int, moab::Range>::iterator ittr;

  for(ittr = tagged_vols_map.begin(); ittr != tagged_vols_map.end(); ++ittr){
    tr_num = ittr->first;
    for(itv = tagged_vols_map[tr_num].begin(); itv != tagged_vols_map[tr_num].end(); ++itv){
      // get the vol's vertices
      surf_set.clear();
      rval = mbi->get_child_meshsets(*itv, surf_set);
      vert_set.clear();
      rval =  mbi->get_entities_by_type(*itv, moab::MBVERTEX, vert_set); 
      MB_CHK_ERR(rval);
      for (its = surf_set.begin(); its != surf_set.end(); its++){
        rval =  mbi->get_entities_by_type(*its, moab::MBVERTEX, vert_set);
        MB_CHK_ERR(rval);
        edge_set.clear();
        rval = mbi->get_child_meshsets(*its, edge_set);
        MB_CHK_ERR(rval);
        for (ite = edge_set.begin(); ite != edge_set.end(); ite++){
          rval =  mbi->get_entities_by_type(*ite, moab::MBVERTEX, vert_set); 
          MB_CHK_ERR(rval);
          rval = mbi->get_child_meshsets(*ite, v_set);
          MB_CHK_ERR(rval);
          for (itvt = v_set.begin(); itvt != v_set.end(); itvt++){
            rval =  mbi->get_entities_by_type(*itvt, moab::MBVERTEX, vert_set);
            MB_CHK_ERR(rval);
          }
        }
      }
      // insert verts into range of all moving verts
      // and into the map of TR #'s to verts
      for (itr = vert_set.begin(); itr != vert_set.end(); itr++){
          tagged_verts.insert(*itr); 
          tagged_verts_map[tr_num].insert(*itr);
      }
    }
  }
  return moab::MB_SUCCESS;
}

// Need to add rotation functionality back in

/* Funtion to rotate a 3D point a distance theta 
   around any line (given by two points)
*/
//XYZ rotate_point(XYZ P, double theta, XYZ L1, XYZ L2)
//{
//   XYZ v, u, q1, q2;
//   double m, d;
//
//   //translate so that rotation axis is origin
//   q1.x = P.x - L1.x;
//   q1.y = P.y - L1.y;
//   q1.z = P.z - L1.z;
//
//   // find unit vector of rotation axis
//   v.x = L2.x - L1.x;
//   v.y = L2.y - L1.y;
//   v.z = L2.z - L1.z;
//   m = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
//   u.x = v.x/m;
//   u.y = v.y/m;
//   u.z = v.z/m;
//
//   // rotate space about x axis so v lies in xz plane
//
//   // length of u projected onto yz plane
//   d = sqrt(u.y*u.y + u.z*u.z);
//
//   //if d = 0, v is already in xz plane
//   if (d != 0) 
//     {
//       q2.x = q1.x;
//       q2.y = q1.y * u.z / d - q1.z * u.y / d;
//       q2.z = q1.y * u.y / d + q1.z * u.z / d;
//     } 
//   else 
//     {
//       q2 = q1;
//     }
//
//   // rotate space about y axis so v lies along z axis
//   q1.x = q2.x * d - q2.z * u.x;
//   q1.y = q2.y;
//   q1.z = q2.x * u.x + q2.z * d;
// 
//   // rotate space by angle theta about z axis
//   q2.x = q1.x * cos(theta) - q1.y * sin(theta);
//   q2.y = q1.x * sin(theta) + q1.y * cos(theta);
//   q2.z = q1.z;
//
//   // inverse of y axis rotation
//   q1.x =   q2.x * d + q2.z * u.x;
//   q1.y =   q2.y;
//   q1.z = - q2.x * u.x + q2.z * d;
//
//   // inverse of x axis rotation
//   if (d != 0) 
//     {
//       q2.x =   q1.x;
//       q2.y =   q1.y * u.z / d + q1.z * u.y / d;
//       q2.z = - q1.y * u.y / d + q1.z * u.z / d;
//     } 
//   else
//     {
//       q2 = q1;
//     }
//
//   // inverse of translation to origin
//   q1.x = q2.x + L1.x;
//   q1.y = q2.y + L1.y;
//   q1.z = q2.z + L1.z;
//
//   // return rotated point
//   return(q1);
//}
void process_input(char* tfilename, 
                   std::map< int, double[5]> &tr_vec_map,
                   int &number_points,
                   double &total_time)
{
  std::ifstream transform_input(tfilename);
  std::string line;
  const char* delimiters = " "; 
  const char* velocity_start_token = "v"; 
  const char* total_time_start_token = "t";
  const char* number_points_start_token = "n";
  int tr_num;

  if (transform_input.is_open())
   {
      while(std::getline(transform_input, line))
        {
          // Skip blank lines in file
	  if (line.length() == 0 ) continue;

          // Tokenize the line
          std::vector<std::string> tokens;
          tokenize(line, tokens, delimiters);
          if (tokens.empty()) continue ; 
  
          // Transformation input       
          if( tokens[0].compare(velocity_start_token ) == 0 && tokens.size() > 1){
              tr_num = atoi(tokens[1].c_str());
              (tr_vec_map[tr_num])[0] = atof(tokens[2].c_str()); // v.x
              (tr_vec_map[tr_num])[1] = atof(tokens[3].c_str()); // v.y
              (tr_vec_map[tr_num])[2] = atof(tokens[4].c_str()); // v.z
              (tr_vec_map[tr_num])[3] = atof(tokens[5].c_str()); // start time
              (tr_vec_map[tr_num])[4] = atof(tokens[6].c_str()); // end time

          }
          if( tokens[0].compare(total_time_start_token ) == 0 && tokens.size() > 1){
            total_time = atof(tokens[1].c_str());
          }
          if( tokens[0].compare(number_points_start_token ) == 0 && tokens.size() > 1){
            number_points = atoi(tokens[1].c_str());
          }
        }
   }
}

void set_transformation_distance(std::map<int, double [5]> tr_vec_map, 
                    int tr_num,
                    double t,
                    XYZ& tr)
{
  tr.x = tr_vec_map[tr_num][0]*t;
  tr.y = tr_vec_map[tr_num][1]*t;
  tr.z = tr_vec_map[tr_num][2]*t;
}

void get_base_filename(std::string path, std::string &base){

  std::string full_filename = path.substr(path.find_last_of("/\\") + 1);
  std::string::size_type const p(full_filename.find_last_of('.'));
  base = full_filename.substr(0, p);
}

int main(int argc, char* argv[]) 
{
  moab::ErrorCode rval; 
 
  char* gfilename = argv[1];
  char* tfilename = argv[2];


  rval = setup(mbi, gfilename);

  //get moving volumes and verts
  std::map<int, moab::Range> tr_vols_map;
  rval = get_tagged_vols(mbi, tr_vols_map);
  moab::Range mv;
  std::map<int, moab::Range> tr_verts_map;
  rval = get_tagged_verts(tr_vols_map, mv, tr_verts_map);MB_CHK_ERR(rval);
 
  // map of vertex eh to original position
  std::map<moab::EntityHandle, XYZ> position;
  rval = get_orig_vert_position(mv, position);

  //process transformation text file info
  std::map<int, double [5]> tr_vec_map;
  int number_points = 1; //default is one time step
  double total_time = 1.0; //default is 1 s
  process_input(tfilename, tr_vec_map, number_points, total_time);
  double time_step_size = total_time/number_points;

  // make sure # trs in geom file == # trs in text file
  if(tr_verts_map.size() != tr_vec_map.size()){
    std::cout<< "Number of transitions found in geometry does not match number found in transformation text file." << std::endl;
    std::cout<< "There are " << tr_verts_map.size() << " transitions in the geometry file." << std::endl;
    std::cout<< "There are " << tr_vec_map.size() << " transitions in the text file." << std::endl;
  }

  //Inital point, updated point
  XYZ p_0;
  double xyz_new[3];

  double t = 0.0; //current time [s]
 
  while( t <= total_time) {
    //for each TR #
    std::map<int, double[5]>::iterator itv;
    for(itv = tr_vec_map.begin(); itv != tr_vec_map.end(); ++itv){
      int tr_num = itv->first;
      double start_time = (tr_vec_map[tr_num])[3];
      double end_time = (tr_vec_map[tr_num])[4];
      //if t is between start and end time, update pos of vol based on this TR
      if( t >= start_time && t <= end_time){
       //create map of TR #'s to motion vectors
       XYZ trans_vec;
       set_transformation_distance(tr_vec_map, tr_num, t, trans_vec);
       //for each vert
       moab::Range::iterator itvt;
       for(itvt =  tr_verts_map[tr_num].begin(); itvt !=  tr_verts_map[tr_num].end(); ++itvt){
         //get original position
         p_0 = position[*itvt];
       
         //set the coordinates of the updated position  
         xyz_new[0] = p_0.x + trans_vec.x;
         xyz_new[1] = p_0.y + trans_vec.y;
         xyz_new[2] = p_0.z + trans_vec.z;

         rval = mbi->set_coords(&(*itvt), 1, xyz_new);
         MB_CHK_ERR(rval);
       }//move each vertex   
      }//if btwn st/et
    }//for each TR #

    std::string base;
    get_base_filename(std::string(gfilename), base);
    std::string filenum = std::to_string(int(t));
    rval = mbi->write_mesh( (filenum+"_"+base+".h5m").c_str());
 
    // update time
    t = t+time_step_size; 
  }//while
      

  delete DAG;
  return 0;
}
