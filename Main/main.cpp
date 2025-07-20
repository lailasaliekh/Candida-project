// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <cassert>
#include <random>
#include <chrono>
#include <memory>
#include <numeric>
#include <omp.h>
#include <filesystem>

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "RandUtil.hpp"
#include "IBacterium.hpp"
#include "SphericalBacteria.hpp"
#include "RodShapedBacteria.hpp"
#include "Candida.hpp"
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "IO.hpp"
#include "PolyBiofilm.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <unistd.h> 
#include <unordered_set>



//---well mixed--
std::vector<IBacterium*> initialiseBiofilm(double linking1, double linking2, 
  int numTypeA, int numTypeB, 
  double centerX, double centerY) {  // One shared center
    // Set the initial conditions
    std::vector<IBacterium*> initial_conditions;
    std::srand(std::time(0));

    // Define wall constraints scaled by the diameter of PA Bacteria
    double xMin = -73 / 1.17;
    double xMax = 73 / 1.17;
    double yMax = 147/ 1.17;
    double yMin=3;
    double normDistance = 4.5; // Minimum spacing
    double maxRadius = 73.0/1.17; // Max spread distance from center
    int maxAttempts = 300; // Limit placement retries

    std::vector<std::pair<double, double>> positions; // Track placed positions

    // auto generateMixedBacteria = [&](int totalCount) {
    //     for (int i = 0; i < totalCount; ++i) {
    //         bool validPosition = false;
    //         double x, y;

    //         int attempts = 0;
    //         while (!validPosition && attempts < maxAttempts) {
    //             // Generate polar coordinates (radius and angle) for mixing
    //             double r = (std::rand() / (double)RAND_MAX) * maxRadius;  // Random radius
    //             double theta = (std::rand() / (double)RAND_MAX) * 2 * constants::pi; // Random angle

    //             // Convert polar to Cartesian
    //             x = centerX + r * std::cos(theta);
    //             y = centerY + r * std::sin(theta);

    //             // Ensure x and y remain inside the walls
    //             if (x < xMin) x = xMin + normDistance;
    //             if (x > xMax) x = xMax - normDistance;
    //             if (y > yMax) y = yMax - normDistance;

    //             // Check for minimum separation from existing bacteria
    //             validPosition = true;
    //             for (const auto& pos : positions) {
    //                 double dx = x - pos.first;
    //                 double dy = y - pos.second;
    //                 if (std::sqrt(dx * dx + dy * dy) < normDistance) {
    //                     validPosition = false;
    //                     break;
    //                 }
    //             }
    //             attempts++;
    //         }

    //         if (attempts == maxAttempts) {
    //             std::cerr << "Warning: Could not place all bacteria without overlap.\n";
    //             break;
    //         }

    //         double angle = (std::rand() / (double)RAND_MAX) * 2 * constants::pi; // Random orientation

    //         // **Interleave Type A and Type B**
    //         bool isTypeA = (i % 2 == 0 && numTypeA > 0) || (numTypeB == 0); 
    //         bool isTypeB = (i % 2 == 1 && numTypeB > 0) || (numTypeA == 0);

    //         if (isTypeA) {
    //             auto* rod = new RodShapedBacterium{
    //                 x, y, 0,  // x, y, z (z=0)
    //                 angle, constants::pi * 0.5, // Random angle
    //                 RodShapedBacterium::mAvgGrwthRate,
    //                 4, // Type A property
    //                 1, 
    //                 Candida::mRadius
    //             };
    //             initial_conditions.push_back(rod);
    //             numTypeA--; 
    //         } 
    //         else if (isTypeB) {
    //             auto* rod = new RodShapedBacterium{
    //                 x, y, 0,  // x, y, z (z=0)
    //                 angle, constants::pi * 0.5, // Random angle
    //                 RodShapedBacterium::mAvgGrwthRate,
    //                 3, // Type B property
    //                 0,
    //                 0.5
    //             };
    //             initial_conditions.push_back(rod);
    //             numTypeB--; 
    //         }

    //         positions.push_back({x, y}); // Store the placed position
    //     }
    // };

    //--------------------------------------------------------
    //-----uniformly distribute cells/segments in the box-----
    //--------------------------------------------------------

    // Create a global random number generator with a unique seed
    std::random_device rd;
    std::mt19937 gen(rd() ^ std::mt19937::result_type(std::time(0)) ^ std::mt19937::result_type(getpid()));
    
    auto generateMixedBacteria = [&](int totalCount) {
        // Define uniform distributions
        std::uniform_real_distribution<double> distX(xMin, xMax);
        std::uniform_real_distribution<double> distY(yMin, yMax);
        std::uniform_real_distribution<double> distAngle(0, 2 * constants::pi);
    
        for (int i = 0; i < totalCount; ++i) {
            bool validPosition = false;
            double x, y, angle;
            std::array<Vec3, 2> poles;
    
            int attempts = 0;
            while (!validPosition && attempts < maxAttempts) {
                // Generate random position inside the square using mt19937
                x = distX(gen);
                y = distY(gen);
                angle = distAngle(gen);
    
                // Create a temporary bacterium to get end poles
                RodShapedBacterium tempCell{
                    x, y, 0, // x, y, z
                    angle, constants::pi * 0.5, // Random orientation
                    RodShapedBacterium::mAvgGrwthRate,
                    5, // Temporary type
                    1, 
                    Candida::mRadius
                };
    
                tempCell.getMyEndVecs(poles[0], poles[1]); // Get end poles
    
                double x1 = poles[0].x, y1 = poles[0].y;
                double x2 = poles[1].x, y2 = poles[1].y;
    
                // Ensure both end poles are inside boundaries
                if (x1 < xMin || x2 > xMax || y1 < yMin || y2 > yMax) {
                    validPosition = false;
                } else {
                    // Check for minimum separation from existing bacteria
                    validPosition = true;
                    for (const auto& pos : positions) {
                        double dx = x - pos.first;
                        double dy = y - pos.second;
                        if (std::sqrt(dx * dx + dy * dy) < normDistance) {
                            validPosition = false;
                            break;
                        }
                    }
                }
                attempts++;
            }
    
            if (attempts == maxAttempts) {
                std::cerr << "Warning: Could not place all bacteria without overlap.\n";
                break;
            }
    
            // **Interleave Type A and Type B**
            bool isTypeA = (i % 2 == 0 && numTypeA > 0) || (numTypeB == 0);
            bool isTypeB = (i % 2 == 1 && numTypeB > 0) || (numTypeA == 0);
    
            if (isTypeA) {
                auto* rod = new RodShapedBacterium{
                    3/1.17, 75/1.17, 0,  // x, y, z (z=0)
                    constants::pi * 0.5, constants::pi * 0.5, // Random angle
                    RodShapedBacterium::mAvgGrwthRate,
                    4, // Type A property
                    1, // if hyphal ca, 0 if yeast-locked ca
                    Candida::mRadius
                };
                initial_conditions.push_back(rod);
                numTypeA--;
            } 
            else if (isTypeB) {
                auto* rod = new RodShapedBacterium{
                  -3/1.17, 75/1.17, 0,  // x, y, z (z=0)
                    constants::pi * 0.5, constants::pi * 0.5, // Random angle
                    RodShapedBacterium::mAvgGrwthRate,
                    3, // Type B property
                    0, // non-chaining PA bacteria
                    0.5
                };
                initial_conditions.push_back(rod);
                numTypeB--;
            }
    
            positions.push_back({x, y}); // Store the placed position
        }
    };
    
  
  

    // **Generate a mixed biofilm in a single area**
    generateMixedBacteria(numTypeA + numTypeB);

    return initial_conditions;
}

std::vector<IBacterium*> initialiseBiofilmFromFile(double linking1, double linking2, const std::string& filename)
{
    std::vector<IBacterium*> initial_conditions;
    std::unordered_map<int, IBacterium*> id_to_bacterium;

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Unable to open file " << filename << '\n';
        return initial_conditions;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string cell_type;
        int cell_id;
        double length, radius, pos_x, pos_y, pos_z, ori_x, ori_y, ori_z;
        std::string lower_link_str, upper_link_str;

        if (!(iss >> cell_type >> cell_id >> length >> radius >> pos_x >> pos_y >> pos_z >> ori_x >> ori_y >> ori_z >> lower_link_str >> upper_link_str))
            continue;

        IBacterium* bacterium = nullptr;
        RodShapedBacterium* rod = nullptr;

        // Determine linking type from radius
        int linkingType = 0;
        // linking based on the type chained or not
        if (lower_link_str == "None" && upper_link_str == "None") {
            linkingType = 0;  // No links --> type 0
        } else {
            linkingType = 1;  // Has links --> type 1
        }
        double angle = std::atan2(ori_y, ori_x);

        rod = new RodShapedBacterium(pos_x, pos_y, pos_z, angle, constants::pi * 0.5,
                                     RodShapedBacterium::mAvgGrwthRate, length,
                                     linkingType, radius);
        rod->mId = cell_id;
        bacterium = rod;

#ifdef CHAINING
        rod->tmp_lower_link_id = (lower_link_str == "None") ? -1 : std::stoi(lower_link_str);
        rod->tmp_upper_link_id = (upper_link_str == "None") ? -1 : std::stoi(upper_link_str);
#endif

        if (bacterium) {
            initial_conditions.push_back(bacterium);
            id_to_bacterium[cell_id] = bacterium;
        }
    }
    file.close();

    // -------------------- Validation --------------------------
    std::unordered_set<int> all_ids;
    std::vector<std::pair<int, int>> link_pairs;
    bool validation_passed = true;

    for (IBacterium* b : initial_conditions) {
        int id = b->getID();
        if (all_ids.count(id)) {
            std::cerr << "Duplicate ID found: " << id << '\n';
            validation_passed = false;
        } else {
            all_ids.insert(id);
        }

#ifdef CHAINING
        auto* rod = dynamic_cast<RodShapedBacterium*>(b);
        if (!rod) continue;

        if (rod->tmp_upper_link_id == id || rod->tmp_lower_link_id == id) {
            std::cerr << "Cell " << id << " is linked to itself.\n";
            validation_passed = false;
        }

        if (rod->tmp_upper_link_id != -1) link_pairs.emplace_back(id, rod->tmp_upper_link_id);
        if (rod->tmp_lower_link_id != -1) link_pairs.emplace_back(id, rod->tmp_lower_link_id);
#endif
    }

    for (auto& [from_id, to_id] : link_pairs) {
        if (!all_ids.count(to_id)) {
            std::cerr << " Invalid link: Bacterium " << from_id << " links to non-existent ID " << to_id << '\n';
            validation_passed = false;
        }
    }

    if (!validation_passed) {
        std::cerr << "Input file validation failed. Aborting simulation.\n";
        std::exit(EXIT_FAILURE);
    } else {
        std::cout << "Input file validation passed.\n";
    }

#ifdef CHAINING
    for (IBacterium* b : initial_conditions) {
        auto* rod = dynamic_cast<RodShapedBacterium*>(b);
        if (!rod) continue;

        if (rod->tmp_upper_link_id != -1) {
            auto it = id_to_bacterium.find(rod->tmp_upper_link_id);
            rod->setUpperLink((it != id_to_bacterium.end()) ? it->second : nullptr);
        }

        if (rod->tmp_lower_link_id != -1) {
            auto it = id_to_bacterium.find(rod->tmp_lower_link_id);
            rod->setLowerLink((it != id_to_bacterium.end()) ? it->second : nullptr);
        }
    }
#endif

    return initial_conditions;
  }
int main(int argc, char const *argv[])
{

  int num_A; //default number of candida
  int num_B; //default number of pa
  double linking1=0;
  double linking2=1.0; //defualt
#ifdef RANDOM_SEED
  std::cout << "setting random seed" << '\n';
  gen_rand.setSeed(
    std::chrono::high_resolution_clock::now().time_since_epoch().count()
  );
#endif

#if defined(CHAINING)
  std::string run_dir; // Run directory
  if ( argc==7 )
  {
    run_dir = argv[1];                              // Run directory
    double kappa           { std::stod( argv[2]) }; // Spring tension
    double bend_rig        { std::stod( argv[3]) }; // Bending rigidity
    double linking1   { std::stod( argv[4]) }; // Probability daughters link
    double linking2   { std::stod( argv[5]) }; 
    double force_thresh    { std::stod( argv[6]) }; // Threshold force before breaking
  }
  else if ( argc==2 )
  {
    run_dir = argv[1];                              // Run directory
  }
 // else
 // {
   // std::cout << "Expected 5 command line arguents! Received " << argc-1 << '\n';
    //std::cout << "Example usage:\n"
        //      << "./main.out run_dir kappa bend_rig linking_prob force_thresh" << '\n';
   // exit(EXIT_FAILURE);
  //}
  if ( !std::filesystem::exists(sim_out_dir) )
  {
    std::filesystem::create_directories(sim_out_dir);
  }

  sim_out_dir += "/" + run_dir + "/";
  int numTypeA = 1;      // Number of ca
  int numTypeB = 1;      // Number of pa
  double centerX = 0.0;   // Shared center X for mixed distribution
  double centerY = 75.0/1.17;   // Shared center Y for mixed distribution
  

// well-mixed initial conditions
  std::vector<IBacterium*> bacteria_population = initialiseBiofilm(linking1, linking2, 
                                                                    numTypeA, numTypeB, 
                                                                    centerX, centerY);

  //FROM FILE for longer runs, reading from file
  std::string dataFilePath = "/storage/datastore-personal/s2507701/Leonado_paper/NewTestCAndida/GeneratedOutput/SimOutput/data_production/VERTICAL_ORI/CAm1_PA1/repeat1/biofilm_00107.dat"; // 
  
  // std::vector<IBacterium*> bacteria_population = initialiseBiofilmFromFile(linking1, linking2,dataFilePath);

  PolyBiofilm pb { bacteria_population };
  pb.runSim();
#else
  std::cout << "Please define one of the following MACROS: CHAINING or nothing" << '\n';
  exit(42);
#endif // End control input parameters

  return 0;
}

