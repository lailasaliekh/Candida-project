// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
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

// Shortest distance between two line segments (p1,q1) and (p2,q2).
static double segmentSegmentDistance(const Vec3& p1, const Vec3& q1,
                                     const Vec3& p2, const Vec3& q2)
{
    const Vec3 d1 { q1 - p1 };
    const Vec3 d2 { q2 - p2 };
    const Vec3 r  { p1 - p2 };
    const double a { dot(d1,d1) }, e { dot(d2,d2) }, f { dot(d2,r) };
    const double eps { 1e-12 };
    double s { 0.0 }, t { 0.0 };

    if ( a<=eps && e<=eps ) return (p1-p2).norm();
    if ( a<=eps ) { t = std::clamp( f/e, 0.0, 1.0 ); }
    else
    {
        const double c { dot(d1,r) };
        if ( e<=eps ) { s = std::clamp( -c/a, 0.0, 1.0 ); }
        else
        {
            const double b { dot(d1,d2) };
            const double denom { a*e - b*b };
            s = ( denom>eps ) ? std::clamp( (b*f - c*e)/denom, 0.0, 1.0 ) : 0.0;
            t = ( b*s + f )/e;
            if      ( t<0.0 ) { t = 0.0; s = std::clamp( -c/a, 0.0, 1.0 ); }
            else if ( t>1.0 ) { t = 1.0; s = std::clamp( (b-c)/a, 0.0, 1.0 ); }
        }
    }
    return ( (p1 + d1*s) - (p2 + d2*t) ).norm();
}

std::vector<IBacterium*> initialiseBiofilm(double linking1, double linking2,
  int numTypeA, int numTypeB,
  double centerX, double centerY,
  double requestedDropletRadius = -1.0)
{
    std::vector<IBacterium*> initial_conditions;

    // Wall constraints scaled by the diameter of PA bacteria
    const double xMin { -73.0 / 1.17 };
    const double xMax {  73.0 / 1.17 };
    const double yMin {   3.0 };
    const double yMax {  147.0 / 1.17 };

    // Surface-to-surface gap demanded between two cells at t=0
    const double minGap { 0.3 };

    // Largest droplet centred on (centerX,centerY) that still fits in the box.
    const double maxFit {
        std::min( std::min( centerX-xMin, xMax-centerX ),
                  std::min( centerY-yMin, yMax-centerY ) )
    };
    double dropletRadius { requestedDropletRadius>0.0 ? requestedDropletRadius
                                                      : maxFit - 4.0 };
    if ( dropletRadius > maxFit )
    {
        std::cerr << "Warning: droplet of radius " << dropletRadius
                  << " does not fit at (" << centerX << ',' << centerY
                  << "); clipping to " << maxFit << ".\n";
        dropletRadius = maxFit;
    }

    const int maxAttempts { 20000 };

    struct PlacedCell { Vec3 p, q; double radius; };
    std::vector<PlacedCell> placed;
    placed.reserve( numTypeA + numTypeB );

    std::random_device rd;
    std::mt19937 gen( rd() ^ std::mt19937::result_type(std::time(0))
                          ^ std::mt19937::result_type(getpid()) );

    struct InitialCellData {
        bool isCandida;
        double length;
        double linkingProb;
        double radius;
    };

    std::vector<InitialCellData> plannedCells;
    plannedCells.reserve( numTypeA + numTypeB );
    for ( int i=0; i<numTypeA; ++i )
        plannedCells.push_back( InitialCellData{ true, 4.0, 1.0, Candida::mRadius } );
    for ( int i=0; i<numTypeB; ++i )
        plannedCells.push_back( InitialCellData{ false, 3.0, 0.0, 0.5 } );
    std::shuffle( plannedCells.begin(), plannedCells.end(), gen );

    std::uniform_real_distribution<double> distAngle( 0.0, 2*constants::pi );
    std::uniform_real_distribution<double> distUnit ( 0.0, 1.0 );

    int numFailed { 0 };

    for ( const auto& cellData : plannedCells )
    {
        bool validPosition { false };
        double x{0}, y{0}, angle{0};
        Vec3 p, q;

        for ( int attempt=0; attempt<maxAttempts && !validPosition; ++attempt )
        {
            // Uniform over the droplet area
            const double radialDistance { dropletRadius*std::sqrt( distUnit(gen) ) };
            const double polarAngle     { distAngle(gen) };
            x = centerX + radialDistance*std::cos( polarAngle );
            y = centerY + radialDistance*std::sin( polarAngle );
            angle = distAngle(gen);

            // End poles of the spherocylinder (no temporary cell needed, so the
            // static id counter is not burnt on rejected trials)
            const Vec3 halfAxis { 0.5*cellData.length*std::cos(angle),
                                  0.5*cellData.length*std::sin(angle), 0.0 };
            p = Vec3{x,y,0.0} - halfAxis;
            q = Vec3{x,y,0.0} + halfAxis;

            // Whole cell (both caps) must sit inside the box: compare against
            // the min/max of BOTH poles, not pole 0 vs pole 1.
            const double cellXMin { std::min(p.x,q.x) - cellData.radius };
            const double cellXMax { std::max(p.x,q.x) + cellData.radius };
            const double cellYMin { std::min(p.y,q.y) - cellData.radius };
            const double cellYMax { std::max(p.y,q.y) + cellData.radius };
            if ( cellXMin<xMin || cellXMax>xMax || cellYMin<yMin || cellYMax>yMax )
                continue;

            // Centre of mass inside the droplet
            const double distanceFromCentre {
                std::hypot( x-centerX, y-centerY )
            };
            if ( distanceFromCentre>dropletRadius ) continue;

            // True rod-rod separation, so long cells cannot interpenetrate
            validPosition = true;
            for ( const auto& other : placed )
            {
                const double surfaceGap {
                    segmentSegmentDistance(p,q,other.p,other.q)
                    - cellData.radius - other.radius
                };
                if ( surfaceGap<minGap ) { validPosition = false; break; }
            }
        }

        if ( !validPosition ) { ++numFailed; continue; } // keep trying the rest

        auto* rod = new RodShapedBacterium{
            x, y, 0,
            angle,
            constants::pi*0.5,
            RodShapedBacterium::mAvgGrwthRate,
            cellData.length,
            cellData.linkingProb,
            cellData.radius
        };
        initial_conditions.push_back( rod );
        placed.push_back( PlacedCell{ p, q, cellData.radius } );
    }

    if ( numFailed>0 )
    {
        std::cerr << "Warning: could only place "
                  << initial_conditions.size() << " of "
                  << ( numTypeA + numTypeB )
                  << " cells in a droplet of radius " << dropletRadius
                  << " (" << numFailed << " failed). "
                  << "Increase the droplet radius or reduce the cell number.\n";
    }
    std::cout << "Initialised " << initial_conditions.size()
              << " cells in a droplet of radius " << dropletRadius
              << " centred at (" << centerX << ", " << centerY << ")\n";

    return initial_conditions;
}

//--------------------------------------------------------
//-----reading from a file-----
//--------------------------------------------------------
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
  double growthRateMultiplier1 { RodShapedBacterium::mHyphalGrowthRateMultiplier };
  int numTypeA = 1;
  int numTypeB = 1;
  if ( argc==10 )
  {
    run_dir = argv[1];                              // Run directory
    double kappa           { std::stod( argv[2]) }; // Spring tension
    double bend_rig        { std::stod( argv[3]) }; // Bending rigidity
    double linking1   { std::stod( argv[4]) }; // Probability daughters link
    double linking2   { std::stod( argv[5]) }; 
    double force_thresh    { std::stod( argv[6]) }; // Threshold force before breaking
    growthRateMultiplier1 = std::stod( argv[7] );
    numTypeA = std::stoi( argv[8] );
    numTypeB = std::stoi( argv[9] );
  }
  else if ( argc==8 )
  {
    run_dir = argv[1];                              // Run directory
    double kappa           { std::stod( argv[2]) }; // Spring tension
    double bend_rig        { std::stod( argv[3]) }; // Bending rigidity
    double linking1   { std::stod( argv[4]) }; // Probability daughters link
    double linking2   { std::stod( argv[5]) }; 
    double force_thresh    { std::stod( argv[6]) }; // Threshold force before breaking
    growthRateMultiplier1 = std::stod( argv[7] );
  }
  else if ( argc==7 )
  {
    run_dir = argv[1];                              // Run directory
    double kappa           { std::stod( argv[2]) }; // Spring tension
    double bend_rig        { std::stod( argv[3]) }; // Bending rigidity
    double linking1   { std::stod( argv[4]) }; // Probability daughters link
    double linking2   { std::stod( argv[5]) }; 
    double force_thresh    { std::stod( argv[6]) }; // Threshold force before breaking
  }
  else if ( argc==3 )
  {
    run_dir = argv[1];                              // Run directory
    growthRateMultiplier1 = std::stod( argv[2] );
  }
  else if ( argc==5 )
  {
    run_dir = argv[1];                              // Run directory
    growthRateMultiplier1 = std::stod( argv[2] );
    numTypeA = std::stoi( argv[3] );
    numTypeB = std::stoi( argv[4] );
  }
  else if ( argc==2 )
  {
    run_dir = argv[1];                              // Run directory
  }
  else
  {
    std::cout << "Unexpected number of command line arguments: " << argc-1 << '\n';
    std::cout << "Example usage:\n"
             << "./main.out run_dir [growthRateMultiplier1]\n"
             << "./main.out run_dir [growthRateMultiplier1] [numTypeA] [numTypeB]\n"
             << "./main.out run_dir kappa bend_rig linking1 linking2 force_thresh [growthRateMultiplier1] [numTypeA] [numTypeB]\n";
    exit(EXIT_FAILURE);
  }

  RodShapedBacterium::mHyphalGrowthRateMultiplier = growthRateMultiplier1;
  std::cout << "growthRateMultiplier1\t: "
            << RodShapedBacterium::mHyphalGrowthRateMultiplier << '\n';
  std::cout << "numTypeA\t: " << numTypeA << '\n';
  std::cout << "numTypeB\t: " << numTypeB << '\n';

  if ( !std::filesystem::exists(sim_out_dir) )
  {
    std::filesystem::create_directories(sim_out_dir);
  }

  sim_out_dir += "/" + run_dir + "/";
  ///////------------------------------
  ////----------------------------------
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
