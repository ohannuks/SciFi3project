//SAMRAI INCLUDES<
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
//SAMRAI INCLUDES>

//BOOST:
#include "boost/shared_ptr.hpp"

#include <iostream>
#include <string>


// Local includes:
#include "poisson.h"

using namespace std;
using namespace SAMRAI;

int main(int argc, char **argv) {
  /*
   * Initialize MPI, SAMRAI.
   */

  tbox::SAMRAI_MPI::init(&argc, &argv);
  tbox::SAMRAIManager::initialize();
  tbox::SAMRAIManager::startup();
  //std::cout << "Hello, world!" << std::endl;
  { // Block to make sure there's no problem with deallocation
    string input_filename;
    if (argc != 2) {
      TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n"
                            << "  options:\n"
                            << "  none at this time" << endl);
    } else {
      input_filename = argv[1];
    }
    
    // Create input database, and parse all data in input file:
    boost::shared_ptr<tbox::InputDatabase> input_db(
       new tbox::InputDatabase("input_db"));
    tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);
    
    //TODO: Find out what timer manager is
    // Set up TIMER MANAGER
    
    if (input_db->isDatabase("TimerManager")) {
       tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
    }
    
    // Parse "MAIN " section from the input database file
    boost::shared_ptr<tbox::Database> main_db(input_db->getDatabase("Main"));
    
    // Get dimensions of the problem from the main database
    const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));
    // Read the base name for the project (save file name)
    string base_name = "unnamed";
    base_name = main_db->getStringWithDefault("base_name", base_name);
        // Log all nodes in a parallel run:
        const string log_file_name = base_name + ".log";
        bool log_all_nodes = false;
        log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
              log_all_nodes);
        if (log_all_nodes) {
           tbox::PIO::logAllNodes(log_file_name);
        } else {
           tbox::PIO::logOnlyNodeZero(log_file_name);
        }
    
    // Define grid geometry (We use Cartesian coordinates (not e.g. spherical)
    boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
    new geom::CartesianGridGeometry(
                                    dim,
                                    base_name + "CartesianGridGeometry",
                                    input_db->getDatabase("CartesianGridGeometry"))
                                   );
        // LOG FILE STUFF NOT INTERESTING
        tbox::plog << "Cartesian Geometry:" << endl;
        grid_geometry->printClassData(tbox::plog);
    
    // Define patch hierarchy for the grid (read from input file)
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
                                                            new hier::PatchHierarchy(
                                                               base_name + "::PatchHierarchy",
                                                               grid_geometry,
                                                               input_db->getDatabase("PatchHierarchy"))
                                                           );
  }
  
  tbox::SAMRAIManager::shutdown();
  tbox::SAMRAIManager::finalize();
  tbox::SAMRAI_MPI::finalize();
  return 0;
}
