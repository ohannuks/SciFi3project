/*
 * Copyright 2015 <copyright holder> <email>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "Poisson.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

using namespace SAMRAI;
using namespace std;

Poisson::Poisson(
  const string& _object_name,
  const tbox::Dimension& _dimensions,
  boost::shared_ptr<solv::LocationIndexRobinBcCoefs>& _boundary_condition_coefficients):
  d_object_name(_object_name),
  d_dimensions(_dimensions),
  d_boundary_condition_coefficients(_boundary_condition_coefficients)
{

   hier::VariableDatabase* variable_database = hier::VariableDatabase::getDatabase();

   /*
    * Define class member context
    */
   d_variable_context = variable_database->getContext(object_name + ":Context");

   /*
    * Register variables with hier::VariableDatabase
    * and get the descriptor indices for those variables.
    * (Descriptor indices mean the ID we use to fetch things from variable_database)
    */
   boost::shared_ptr<pdat::CellVariable<double> > computed_solution(
      new pdat::CellVariable<double>(
         d_dimensions,
         d_object_name + ":computed solution",
         1));
   
   const int ghost_cell_stencil_width = 1;
   
   d_computed_solution_id =
      variable_database->registerVariableAndContext(
         computed_solution,
         d_variable_context,
         hier::IntVector(d_dimensions, ghost_cell_stencil_width) /* ghost cell width is 1 for stencil widths */);

   /*
    * Also register right hand side variable into the variable data base and save the id (id is used to fetch variables from the variable database)
    */
   boost::shared_ptr<pdat::CellVariable<double> > right_hand_side_variable(
      new pdat::CellVariable<double>(
         d_dimensions,
         d_object_name
         + ":linear system right hand side"));

   /*
    * 
    */
   d_right_hand_side_id =
      variable_database->registerVariableAndContext(
         right_hand_side_variable,
         d_variable_context,
         hier::IntVector(d_dimensions, 0) /* ghost cell width is 0 */);
}

Poisson::~Poisson()
{

}

void set_right_hand_side_data(hier::Box & patch_box,
                              boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
                              boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry,
                              boost::shared_ptr<pdat::CellData<double> > right_hand_side_data) {
  // Get the cell lengths in the patch (in here grid geometry is patch geometry because we dont have mesh refinement)
  const int dimensions = 3;
  const double cell_lengths[dimensions] = grid_geometry->getDx();
  const double dx = cell_lengths[0];
  const double dy = cell_lengths[1];
  const double dz = cell_lengths[2];
  
  // Get minimum coordinates (these will be used for iteration)
  const double cell_min_coordinates[dimensions] = patch_geometry->getXLower();
  const double min_x = cell_min_coordinates[0];
  const double min_y = cell_min_coordinates[1];
  const double min_z = cell_min_coordinates[2];
  
  // Get the indices of the cells in x- y- and z-direction (will be used for iteration)
  const hier::Index cell_indices_lower = patch_box.lower();
  const hier::Index cell_indices_upper = patch_box.upper();

  //patch_box.getLocalId();
  //patch_box.
  //for( int i = cell_indices_lower[0]; i 
  //patch_geometry->
}


// Initializer function for patch data (this is called whenever a new patch is created)
void Poisson::initializeLevelData(const boost::shared_ptr< hier::PatchHierarchy >& hierarchy, 
                                  const int level_number, 
                                  const double init_data_time, 
                                  const bool can_be_refined, 
                                  const bool initial_time, 
                                  const boost::shared_ptr< hier::PatchLevel > & old_level, 
                                  const bool allocate_data)
{
  // Call null_use on unused data ( this is so we don't get compiler warnings )
  NULL_USE( init_data_time );
  NULL_USE( can_be_refined );
  NULL_USE( initial_time );
  NULL_USE( old_level );
  
  // Copy the patch hierarchy pointer (this I guess needs to be done)
  boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy = hierarchy;
  
  // Get the grid geometry
  boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
     BOOST_CAST<geom::CartesianGridGeometry, hier::BaseGridGeometry>(
        patch_hierarchy->getGridGeometry()));
  
  TBOX_ASSERT(grid_geometry); // Error checking
  
  // Get the patch level data
  boost::shared_ptr<hier::PatchLevel> level(
   hierarchy->getPatchLevel(level_number));
  
  // If required, allocate all patch data on the level.
  if (allocate_data) {
    level->allocatePatchData(d_computed_solution_id);
    level->allocatePatchData(d_right_hand_side_id);
  }
  
  /*
   * Initialize data in all patches in the level.
   */
  // Iterate through all patches in the given level
  // 
  for (hier::PatchLevel::iterator pi(level->begin());
       pi != level->end(); ++pi) {
    
    // Get the next patch in iteration
    const boost::shared_ptr<hier::Patch>& patch = *pi;
  
    // Check for errors:
    if (!patch) {
      TBOX_ERROR(d_object_name
         << ": Cannot find patch.  Null patch pointer.");
    }
    
    // Get the box (dimensions, etc, of the box)
    hier::Box patch_box = patch->getBox();
  
    // Get the patch geomety
    boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch->getPatchGeometry()));

    // Get the right hand side data:
    boost::shared_ptr<pdat::CellData<double> > right_hand_side_data(
      BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
         patch->getPatchData(d_right_hand_side_id)));
    
    // Check for errors:
    TBOX_ASSERT(patch_geometry);
    TBOX_ASSERT(right_hand_side_data);

    /*
     * Set source function and exact solution.
     */
    //TODO: Put grid geometry -> patch geometry for AMR
    set_right_hand_side_data( 
      patch_box,
      grid_geometry,
      patch_geometry,
      right_hand_side_data
    );

  }    // End patch loop.
}





