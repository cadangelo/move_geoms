#include <iostream>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Types.hpp"
#include "moab/GeomTopoTool.hpp"
#include "MBTagConventions.hpp"
#include "DagMC.hpp"
#include "dagmcmetadata.hpp"

#include <string>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

struct XYZ{
  double x;
  double y;
  double z;
};

moab::Tag category_tag;
moab::Tag geom_tag;
//moab::Tag name_tag;
//moab::Tag obj_name_tag;
moab::Tag dim_tag, id_tag;
moab::Tag move_tag;
//moab::Tag sense_tag;
//moab::Tag obb_tag;
//moab::Tag obb_tree_tag;

moab::DagMC *DAG;

dagmcMetaData* DMD;

moab::Core *mbi = new moab::Core;


moab::ErrorCode get_all_handles(moab::Core *mbi)
{
  moab::ErrorCode rval;

//  rval = mbi->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
//			      name_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//  MB_CHK_ERR(rval);
//
//  rval = mbi->tag_get_handle( "OBJECT_NAME", 32, moab::MB_TYPE_OPAQUE,
//			      obj_name_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//  MB_CHK_ERR(rval);

  rval = mbi->tag_get_handle( "MOVE_TAG", 32, moab::MB_TYPE_OPAQUE,
			      move_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
  MB_CHK_ERR(rval);

  int negone = -1;
  rval = mbi->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
			      geom_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT,&negone);
  MB_CHK_ERR(rval);

  rval = mbi->tag_get_handle( GLOBAL_ID_TAG_NAME,
			      1, moab::MB_TYPE_INTEGER,
			      id_tag,
			      moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  MB_CHK_ERR(rval);
  
  rval = mbi->tag_get_handle( CATEGORY_TAG_NAME,
			      CATEGORY_TAG_SIZE,
			      moab::MB_TYPE_OPAQUE,
			      category_tag,
			      moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );

  MB_CHK_ERR(rval);

//  rval = mbi->tag_get_handle("GEOM_SENSE_2", 2, moab::MB_TYPE_HANDLE,
//                             sense_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
//  MB_CHK_ERR(rval);

  return moab::MB_SUCCESS;
}

moab::ErrorCode setup(moab::Core *mbi, char* filename)
{
  moab::ErrorCode rval;

  // get all moab tag handles 
  rval = get_all_handles(mbi);
  MB_CHK_ERR(rval);

  DAG = new moab::DagMC(mbi);
  moab::GeomTopoTool *GTT = new moab::GeomTopoTool(mbi);

  // load base geometry file that we wish to move
  rval = mbi->load_file(filename);
  MB_CHK_ERR(rval);

  rval = DAG->load_existing_contents();
  MB_CHK_ERR(rval);

//  rval = DAG->init_OBBTree();
//  MB_CHK_ERR(rval);

  // find all geometry sets
  rval = GTT->find_geomsets();
  MB_CHK_SET_ERR(rval, "GeomTopoTool could not find the geometry sets");

  // implicit compliment
  // EntityHandle implicit_complement;
  //  rval = GTT->get_implicit_complement(implicit_complement, true);
  rval = DAG->setup_impl_compl();
  MB_CHK_SET_ERR(rval, "Failed to setup the implicit compliment");


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
                                std::string tag_name,
                                std::string tag_delims,
                                std::map<int, moab::Range> &tagged_vols_map)
{
  moab::ErrorCode rval;

  // get tr tag from parse_properties
  std::vector<std::string> group_name;
  std::map<std::string, std::string> group_name_synonyms;
  group_name.push_back(tag_name);
  rval = DAG->parse_properties(group_name, group_name_synonyms, tag_delims.c_str());
  MB_CHK_SET_ERR(rval, "DAGMC failed to parse metadata properties");

  // get desired tagged entities
  moab::EntityHandle tagged_meshset;
  moab::Range surf_set, vert_set;
  int num_verts;
  moab::Range::iterator it, itr;
  rval = mbi->create_meshset(moab::MESHSET_SET, tagged_meshset);
  MB_CHK_ERR(rval);

  // get total number of volumes
  int all_vols = DAG->num_entities( 3 );

  // for each DAG volume
  for( int i = 1; i <= all_vols; ++i ){ 
   //get DAG index
   moab::EntityHandle vol = DAG->entity_by_index( 3, i );

   if( DAG->has_prop( vol, tag_name)){
     std::vector<std::string> properties;
     rval = DAG->prop_values(vol, tag_name, properties);

     // get the value of each property 
     for (int j = 0 ; j < properties.size() ; j++) {
       int val = std::stoi(properties[j]);
       bool added = false;
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
            }
            // if added not true, keep looping, if true exit loop
            if(added == true)
              break;
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
  moab::Range surf_set, vert_set;
  moab::Range::iterator its, itv, itr; 
  std::map<int, moab::Range>::iterator ittr;

  for(ittr = tagged_vols_map.begin(); ittr != tagged_vols_map.end(); ++ittr){
    tr_num = ittr->first;
    std::cout << "TR # " << tr_num << std::endl;
    for(itv = tagged_vols_map[tr_num].begin(); itv != tagged_vols_map[tr_num].end(); ++itv){
      std::cout << "VOL # " << *itv << std::endl;
      // get the vol's vertices
      surf_set.clear();
      rval = mbi->get_child_meshsets(*itv, surf_set);
      MB_CHK_ERR(rval);
      for (its = surf_set.begin(); its != surf_set.end(); its++){
        vert_set.clear();
        rval =  mbi->get_entities_by_type(*its, moab::MBVERTEX, vert_set);
        MB_CHK_ERR(rval);
        // insert verts into range of all moving verts
        // and into the map of TR #'s to verts
        for (itr = vert_set.begin(); itr != vert_set.end(); itr++){
            tagged_verts.insert(*itr);
            tagged_verts_map[tr_num].insert(*itr);
        }
      }
    }
   // std::cout << "num verts in this tr num " << tagged_verts_map[tr_num].size() << std::endl;
      std::cout << "first gtverts " << *tagged_verts_map[tr_num].begin() << std::endl;
  }
    std::cout << "num verts in this tr num " << tagged_verts_map[tr_num].size() << std::endl;
  std::cout << "num all moving verts" << tagged_verts.size() << std::endl;
  return moab::MB_SUCCESS;
}

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
void process_input(char* tfilename, 
                   std::map< int, double[12]> &tr_vec_map,
                   std::map<int, std::vector<int>> &step_tr_map,
                   XYZ& v_0, XYZ& b, XYZ& c, XYZ& d, 
                   XYZ& L0, XYZ& L1, double& omega_0, double alpha[3],
                   double& ts, double& end_t, int& rotation, int& translation)
{
  std::ifstream transform_input(tfilename);
  std::string line;
  const char* delimiters = " "; 
  const char* velocity_start_token = "v"; 
  const char* a0_start_token = "b";
  const char* a1_start_token = "c";
  const char* a2_start_token = "d";
  const char* L0_start_token = "i";
  const char* L1_start_token = "j";
  const char* ang_velocity_start_token = "w";
  const char* ang_acc_start_token = "a";
  const char* time_step_start_token = "s";
  const char* end_time_start_token = "e";
  const char* rotation_start_token = "r";
  const char* translation_start_token = "x";
  const char* mcnp_start_token = "t";
  int tr_num;
  int step_num;
  std::vector<int> tr_nums;
//  std::map<int, std::vector<int>> step_tr_map;

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
  
          // MCNP TR card input       
          if( tokens[0].compare(mcnp_start_token ) == 0 && tokens.size() > 1){
              step_num = atoi(tokens[1].c_str());
              tr_num = atoi(tokens[2].c_str());
              step_tr_map[step_num].push_back(tr_num);
              std::cout << "s t m " << step_num << " " << step_tr_map[step_num].size() << std::endl; 
              (tr_vec_map[tr_num])[0] = atof(tokens[3].c_str());
              (tr_vec_map[tr_num])[1] = atof(tokens[4].c_str());
              (tr_vec_map[tr_num])[2] = atof(tokens[5].c_str());

              //if reading from TR card style input, set end_t and ts to 1
              end_t = 1.0;
              ts = 1.0; 

          }
          // Initial velocity-- need to update the rest of these to be like MCNP TR input above
          // first token should be TR number 
          // create map of TR nums to motion vectors
          // will also need to add these to set params fxn
          // can only have one end_t and ts per transformation text file, so no map for these, obvs
          if( tokens[0].compare(velocity_start_token ) == 0 && tokens.size() > 1)
            {
             // v_0.x = atof(tokens[1].c_str());
             // v_0.y = atof(tokens[2].c_str());
             // v_0.z = atof(tokens[3].c_str());
              tr_num = atoi(tokens[1].c_str());
              (tr_vec_map[tr_num])[0] = atof(tokens[2].c_str());
              (tr_vec_map[tr_num])[1] = atof(tokens[3].c_str());
              (tr_vec_map[tr_num])[2] = atof(tokens[4].c_str());
            }
          // Acceleration 
          // a(t) = b + ct + dt^2
          if( tokens[0].compare(a0_start_token ) == 0 && tokens.size() > 1)
            {
             // b.x = atof(tokens[1].c_str());
             // b.y = atof(tokens[2].c_str());
             // b.z = atof(tokens[3].c_str());
              tr_num = atoi(tokens[1].c_str());
              (tr_vec_map[tr_num])[3] = atof(tokens[2].c_str());
              (tr_vec_map[tr_num])[4] = atof(tokens[3].c_str());
              (tr_vec_map[tr_num])[5] = atof(tokens[4].c_str());
            }
          if( tokens[0].compare(a1_start_token ) == 0 && tokens.size() > 1)
            {
             // c.x = atof(tokens[1].c_str());
             // c.y = atof(tokens[2].c_str());
             // c.z = atof(tokens[3].c_str());
              tr_num = atoi(tokens[1].c_str());
              (tr_vec_map[tr_num])[6] = atof(tokens[2].c_str());
              (tr_vec_map[tr_num])[7] = atof(tokens[3].c_str());
              (tr_vec_map[tr_num])[8] = atof(tokens[4].c_str());
            }
          if( tokens[0].compare(a2_start_token ) == 0 && tokens.size() > 1)
            {
             // d.x = atof(tokens[1].c_str());
             // d.y = atof(tokens[2].c_str());
             // d.z = atof(tokens[3].c_str());
              tr_num = atoi(tokens[1].c_str());
              (tr_vec_map[tr_num])[9] = atof(tokens[2].c_str());
              (tr_vec_map[tr_num])[10] = atof(tokens[3].c_str());
              (tr_vec_map[tr_num])[11] = atof(tokens[4].c_str());
            }
          // Two points that define line points rotate about 
          if( tokens[0].compare(L0_start_token ) == 0 && tokens.size() > 1)
            {
              L0.x = atof(tokens[1].c_str());
              L0.y = atof(tokens[2].c_str());
              L0.z = atof(tokens[3].c_str());
            }
          if( tokens[0].compare(L1_start_token ) == 0 && tokens.size() > 1)
            {
              L1.x = atof(tokens[1].c_str());
              L1.y = atof(tokens[2].c_str());
              L1.z = atof(tokens[3].c_str());
            }
          //  Angular Acceleration         
          // alpha(t) = alpha[0]+ alpha[1]*t + alpha[2]t^2
          if( tokens[0].compare(ang_acc_start_token ) == 0 && tokens.size() > 1)
            {
              alpha[0] = atof(tokens[1].c_str());
              alpha[1] = atof(tokens[2].c_str());
              alpha[2] = atof(tokens[3].c_str());
            }
          if( tokens[0].compare(ang_velocity_start_token ) == 0 && tokens.size() > 1)
            {
              omega_0 = atof(tokens[1].c_str());
            }
          if( tokens[0].compare(time_step_start_token ) == 0 && tokens.size() > 1)
            {
              ts = atof(tokens[1].c_str());
              std::cout << "time step " << ts << std::endl;
            }
          if( tokens[0].compare(end_time_start_token ) == 0 && tokens.size() > 1)
            {
              end_t = atof(tokens[1].c_str());
              std::cout << "end time " << end_t << std::endl;
            }
          if( tokens[0].compare(rotation_start_token ) == 0 && tokens.size() > 1)
            {
              rotation = stoi(tokens[1]);
            }
          if( tokens[0].compare(translation_start_token ) == 0 && tokens.size() > 1)
            {
              translation = stoi(tokens[1]);
            }
        }
   }
}

void set_parameters(std::map<int, double [12]> tr_vec_map, 
                    int tr_num,
                    XYZ& v_0, XYZ& b, XYZ& c, XYZ& d, 
                    XYZ& L0, XYZ& L1, double& omega_0, double alpha[3],
                    double& ts, double& end_t, int& rotation, int& translation)
{

  v_0.x = tr_vec_map[tr_num][0];
  v_0.y = tr_vec_map[tr_num][1];
  v_0.z = tr_vec_map[tr_num][2];
  
  b.x = tr_vec_map[tr_num][3];
  b.y = tr_vec_map[tr_num][4];
  b.z = tr_vec_map[tr_num][5];

  c.x = tr_vec_map[tr_num][6];
  c.y = tr_vec_map[tr_num][7];
  c.z = tr_vec_map[tr_num][8];

  d.x = tr_vec_map[tr_num][9];
  d.y = tr_vec_map[tr_num][10];
  d.z = tr_vec_map[tr_num][11];


  if( v_0.x + v_0.y + v_0.z != 0.0){
    translation = 1;
  }

}


int main(int argc, char* argv[]) 
{
  moab::ErrorCode rval; 
 
  char* gfilename = argv[1];
  char* tfilename = argv[2];

  rval = setup(mbi, gfilename);


  //get moving volumes and verts
  moab::Range vols;
  std::map<int, moab::Range> tr_vols_map;
  std::string tag_name = "tr";
  std::string tag_delims = ":";
  rval = get_tagged_vols(mbi, tag_name, tag_delims, tr_vols_map);
 

  //process input function
  XYZ v, v_0;
  XYZ b, c, d;
  XYZ L0, L1;
  double omega_0;
  double alpha[3] = {};
  double theta;
  double ts; //length of time step [s]
  double end_t; //end time [s]
  int translation;
  int rotation; 
  std::map<int, double [12]> tr_vec_map;
  std::map<int, double [12]>::iterator itrv;
  std::map<int, std::vector<int>>::iterator its;
  std::map<int, std::vector<int>> step_tr_map;
  std::vector<int> tr_nums;
  process_input(tfilename, tr_vec_map, step_tr_map, v_0, b, c, d, L0, L1, omega_0, alpha, ts, end_t, rotation, translation);
  
  moab::Range mv;
  std::map<int, moab::Range> tr_verts_map;
  rval = get_tagged_verts(tr_vols_map, mv, tr_verts_map);MB_CHK_ERR(rval);
      std::cout << "num verts in tr3 " << tr_verts_map[3].size() << std::endl;
 
  
// map of vertex eh to original position
//  std::map<moab::EntityHandle, XYZ> orig_positions;
  std::map<moab::EntityHandle, XYZ> position;
//  std::map<int, moab::Range> tr_verts_map;
//  rval = get_moving_verts(tr_vols_map, mv, tr_verts_map);
  rval = get_orig_vert_position(mv, position);
  std::cout << "position " << position[*mv.begin()].x << std::endl;

  //Inital point, updated point
  XYZ p_0, p, p_new;
  double xyz[3], xyz_new[3];

  double t = 0.0; //current time [s]
  int shot_num = 0; //current time step

  int tr_num;
  int step_num;
 
  //base output file name 
  std::string output_file = "moved.h5m";

  moab::Range::iterator itt, itv, itx, itz, itvt;
  std::map<int, moab::Range>::iterator ittr;
  std::cout << "t step " << ts << std::endl;
  std::cout << "end time " << end_t << std::endl;

  std::cout << "verts map size1 " << tr_verts_map.size() << std::endl;
  // make sure # trs in geom file == # trs in text file
  if(tr_verts_map.size() != tr_vec_map.size()){
    std::cout<< "Number of transitions found in geometry does not match number found in transformation text file." << std::endl;
    std::cout<< "There are " << tr_vols_map.size() << " transitions in the geometry file." << std::endl;
    std::cout<< "There are " << tr_vec_map.size() << " transitions in the text file." << std::endl;
  }

      std::cout << "verts map size2 " << tr_verts_map.size() << std::endl;
  //set end_t to 1 for relocation translations (single time step)
  // for velocity vectors that need broken into time steps,
  // set start_t and end_t to full length of time
//  while (t <= end_t){
//    std::cout << "time step t= " << t << std::endl;
    //for each time step 
    for(its = step_tr_map.begin(); its != step_tr_map.end(); ++its){
      tr_nums = step_tr_map[its->first];
      step_num =  its->first;
      std::cout << "step # " << step_num << std::endl;
      std::cout << "verts map size3 " << tr_verts_map.size() << std::endl;
    //for each TR #
    for(int i = 0; i < tr_nums.size(); i++){
       tr_num = tr_nums[i];
       std::cout << "tr num " << tr_num << std::endl;
       
       //create map of TR #'s to motion vectors
       // then can move each vert w/ correct TR
  //     set_parameters(tr_vec_map, tr_num, v_0, b, c, d, L0, L1, omega_0, alpha, ts, end_t, rotation, translation);
  v_0.x = tr_vec_map[tr_num][0];
  v_0.y = tr_vec_map[tr_num][1];
  v_0.z = tr_vec_map[tr_num][2];
      std::cout << "v 0 x " << v_0.x << std::endl;
      //for each vert
      std::cout << "verts map size " << tr_verts_map.size() << std::endl;
      std::cout << "# verts in this tr " << tr_verts_map[tr_num].size() << std::endl;
      std::cout << "verts x " << position[*tr_verts_map[tr_num].end()].x << std::endl;
      std::cout << "orig x pos " << p.x << std::endl;
      for(itvt =  tr_verts_map[tr_num].begin(); itvt !=  tr_verts_map[tr_num].end(); ++itvt){
        std::cout << "vert # " << *itvt << std::endl;
        //get original position
        p_0 = position[*itvt];
        std::cout << "p 0 x " << p_0.x << std::endl;
      
        //if translation
     //   if (translation == 1){
             //p_new.x = p_0.x + v_0.x*t + (1/2)*b.x*pow(t,2) + (1/6)*c.x*pow(t,3) + (1/12)*d.x*pow(t,4);
             p_new.x = p_0.x + v_0.x;
             //p_new.y = p_0.y + v_0.y*t + (1/2)*b.y*pow(t,2) + (1/6)*c.y*pow(t,3) + (1/12)*d.y*pow(t,4);
             p_new.y = p_0.y + v_0.y;
             //p_new.z = p_0.z + v_0.z*t + (1/2)*b.z*pow(t,2) + (1/6)*c.z*pow(t,3) + (1/12)*d.z*pow(t,4);
             p_new.z = p_0.z + v_0.z;

   //     }
        
   //     //if rotation
   //     if (rotation == 1){ 
   //         double omega = omega_0 + alpha[0]*t + (1/2)*alpha[1]*pow(t,2) + (1/3)*alpha[2]*pow(t,3);
   //         theta = omega*t;
   //         p_new = rotate_point(p_0, theta, L0, L1);
   //     }
      
        //set the coordinates of the updated position  
        xyz_new[0] = p_new.x;
        xyz_new[1] = p_new.y;
        xyz_new[2] = p_new.z;
        rval = mbi->set_coords(&(*itvt), 1, xyz_new);
        MB_CHK_ERR(rval);

      }//move each vertex   
 
      std::cout << "p orig " << p_0.x << " " << p_0.y <<  std::endl;
      std::cout << "p new " << xyz_new[0] << " " << xyz_new[1] << std::endl;
    }//for each TR #
      
    std::string filenum;
    //filenum = std::to_string(tr_num)+std::to_string(shot_num);
    filenum = std::to_string(step_num);
    //rval = mbi->write_mesh( (std::to_string(shot_num)+output_file).c_str());
    //rval = mbi->write_mesh( (filenum+output_file).c_str());
    std::cout << "file num " << filenum << std::endl;
//    t = t + ts;
    } // for each time step
      
//  }//while

  delete DAG;
  return 0;
}
