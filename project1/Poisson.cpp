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
   d_variable_context = variable_database->getContext(_object_name + ":Context");

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

void Poisson::resetHierarchyConfiguration(
     const boost::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
     int coarsest_level,
     int finest_level) {
  return;
}

void set_right_hand_side_data(hier::Box & patch_box,
                              boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
                              boost::shared_ptr<geom::CartesianPatchGeometry> patch_geometry,
                              boost::shared_ptr<pdat::CellData<double> > right_hand_side_data) {
  // Get the cell lengths in the patch (in here grid geometry is patch geometry because we dont have mesh refinement)
  const int dimensions = 3;
  const double * cell_lengths = grid_geometry->getDx();
  const double dx = cell_lengths[0];
  const double dy = cell_lengths[1];
  const double dz = cell_lengths[2];
  
  // Get minimum coordinates (these will be used for iteration)
  const double * patch_min_coordinates = patch_geometry->getXLower();
  const double min_x = patch_min_coordinates[0];
  const double min_y = patch_min_coordinates[1];
  const double min_z = patch_min_coordinates[2];
  
  // Get the indices of the cells in x- y- and z-direction (will be used for iteration)
  const hier::Index cell_indices_lower = patch_box.lower();
  const hier::Index cell_indices_upper = patch_box.upper();
  // i = x-dimension indice, j = y-dimension indice, k = z-dimension indice
  const int i_min = cell_indices_lower[0];
  const int i_max = cell_indices_upper[0];
  const int j_min = cell_indices_lower[1];
  const int j_max = cell_indices_upper[1];
  const int k_min = cell_indices_lower[2];
  const int k_max = cell_indices_upper[2];
  
  double * rho = right_hand_side_data->getPointer();
  
  int i, j, k;
  for( k = k_min; k <= k_max; ++k ) {
    for ( j = j_min; j <= j_max; ++j ) {
      for ( i = i_min; i < i_max; ++i ) {
        const double x = min_x + 0.5*dx + i*dx;
        const double y = min_y + 0.5*dy + i*dy;
        const double z = min_z + 0.5*dz + i*dz;
        rho[k + (k_max - k_min)*j + (k_max - k_min)*(j_max - j_min)*i] = 0;
      }
    }
  }
  
  

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

/*
 *************************************************************************
 * Set up VisIt to plot internal data from this class.
 * Tell the plotter about the refinement ratios.  Register variables
 * appropriate for plotting.
 *************************************************************************
 */
int Poisson::registerVariablesWithPlotter(
   appu::VisItDataWriter& visit_writer) const {

   /*
    * This must be done once.
    */
   if (!d_patch_hierarchy) {
      TBOX_ERROR(
         d_object_name << ": No hierarchy in\n"
                       << " HyprePoisson::registerVariablesWithPlotter\n"
                       << "The hierarchy must be built before calling\n"
                       << "this function.\n");
   }
   
   /*
    * Register variables with plotter.
    */
   visit_writer.registerPlotQuantity("Computed solution",
      "SCALAR",
      d_computed_solution_id);
   visit_writer.registerDerivedPlotQuantity("Error",
      "SCALAR",
      (appu::VisDerivedDataStrategy *)this);
   
   visit_writer.registerPlotQuantity("Poisson source",
      "SCALAR",
      d_right_hand_side_id);
   visit_writer.registerDerivedPlotQuantity("Patch level number",
      "SCALAR",
      (appu::VisDerivedDataStrategy *)this);

   return 0;
}



/*
 *************************************************************************
 * Solve the Poisson problem.
 *************************************************************************
 */
bool Poisson::solvePoisson()
{

   if (!d_patch_hierarchy) {
      TBOX_ERROR("Cannot solve using an uninitialized object.\n");
   }

   const int level_number = 0;

   /*
    * Fill in the initial guess and Dirichlet boundary condition data.
    * For this example, we want u=0 on all boundaries.
    * The easiest way to do this is to just write 0 everywhere,
    * simultaneous setting the boundary values and initial guess.
    */
   boost::shared_ptr<hier::PatchLevel> level(d_patch_hierarchy->getPatchLevel(
                                                level_number));
   // Iterate through patches
   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
     // Get the patch
     const boost::shared_ptr<hier::Patch>& patch = *ip;
     // Get the patch data
     boost::shared_ptr<pdat::CellData<double> > data(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_computed_solution_id)));
     // Error checking
     TBOX_ASSERT(data);
     // Set zeros as the initial state
     data->fill(0.0);
   }
   
   return true;
//   // d_poisson_hypre->setBoundaries( "Dirichlet" );
//   d_poisson_hypre->setPhysicalBcCoefObject(d_bc_coefs.get());
//
//   /*
//    * Set up HYPRE solver object.
//    * The problem specification is set using the
//    * CellPoissonSpecifications object then passed to the solver
//    * for setting the coefficients.
//    */
//   d_poisson_hypre->initializeSolverState(d_hierarchy,
//      level_number);
//   solv::PoissonSpecifications sps("Hypre Poisson solver");
//   sps.setCZero();
//   sps.setDConstant(1.0);
//   d_poisson_hypre->setMatrixCoefficients(sps);
//
//   /*
//    * Solve the system.
//    */
//   tbox::plog << "solving..." << std::endl;
//   int solver_ret;
//   solver_ret = d_poisson_hypre->solveSystem(d_comp_soln_id,
//         d_rhs_id);
//   /*
//    * Present data on the solve.
//    */
//   tbox::plog << "\t" << (solver_ret ? "" : "NOT ") << "converged " << "\n"
//              << "      iterations: " << d_poisson_hypre->getNumberOfIterations()
//              << "\n"
//              << "      residual: " << d_poisson_hypre->getRelativeResidualNorm()
//              << "\n"
//              << std::flush;
//
//   /*
//    * Deallocate state.
//    */
//   d_poisson_hypre->deallocateSolverState();
//
//   /*
//    * Return whether solver converged.
//    */
//   return solver_ret ? true : false;
}