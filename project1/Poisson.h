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

#ifndef POISSON_H
#define POISSON_H

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/solv/CellPoissonHypreSolver.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#include "boost/shared_ptr.hpp"

namespace SAMRAI {


class Poisson:
   public mesh::StandardTagAndInitStrategy
{
public:
  Poisson(
          const std::string& object_name,
          const tbox::Dimension& dim,
          boost::shared_ptr<solv::LocationIndexRobinBcCoefs>& bc_coefs);
  
//  Poisson(const Poisson& other);
  ~Poisson();
  
  //@{ @name mesh::StandardTagAndInitStrategy virtuals

  /*!
   * @brief Allocate and initialize data for a new level 
   * in the patch hierarchy.
   *
   * This is where you implement the code for initialize data on 
   * the grid.  All the information needed to initialize the grid 
   * are in the arguments.
   *
   * @see mesh::StandardTagAndInitStrategy::initializeLevelData()
   */
  virtual void
  initializeLevelData(
     const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
     const int level_number,
     const double init_data_time,
     const bool can_be_refined,
     const bool initial_time,
     const boost::shared_ptr<hier::PatchLevel>& old_level =
        boost::shared_ptr<hier::PatchLevel>(),
     const bool allocate_data = true);

  /*!
   * @brief Reset any internal hierarchy-dependent information.
   */
  virtual void
  resetHierarchyConfiguration(
     const boost::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
     int coarsest_level,
     int finest_level);

  //@}
  
  /*!
   * @brief Set up external plotter to plot internal
   * data from this class.
   *
   * After calling this function, the external
   * data writer may be used to write the
   * visit file for this object.
   *
   * The internal hierarchy is used and must be
   * established before calling this function.
   * (This is commonly done by building a hierarchy
   * with the mesh::StandardTagAndInitStrategy virtual
   * functions implemented by this class.)
   *
   * @param visit_writer VisIt data writer
   */
  int
  registerVariablesWithPlotter(
     appu::VisItDataWriter& visit_writer) const;
     
  /*!
   * @brief Solve using custom Poisson solver
   *
   * Set up the linear algebra problem and use a
   * solv::CellPoissonHypreSolver object to solve it.
   * -# Set initial guess
   * -# Set boundary conditions
   * -# Specify Poisson equation parameters
   * -# Call solver
   *
   * @return whether solver converged
   */
  bool
  solvePoisson();

  
private:
  std::string d_object_name;

  const tbox::Dimension d_dimensions;

  boost::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
  
  //@{
  /*!
   * @name Major algorithm objects.
   */

  /*!
   * @brief Boundary condition coefficient implementation.
   */
  boost::shared_ptr<solv::LocationIndexRobinBcCoefs> d_boundary_condition_coefficients;

  //@}
  
  /*!
 * @name Private state variables for solution.
 */

 /*!
  * @brief Context owned by this object.
  */
 boost::shared_ptr<hier::VariableContext> d_variable_context;

 /*!
  * @brief Descriptor indices of internal data.
  *
  * These are initialized in the constructor and never change.
  */
 int d_computed_solution_id, d_right_hand_side_id;
};

}

#endif // POISSON_H
