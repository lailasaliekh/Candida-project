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
#include "Geometry.hpp"
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

/*
  Minimum distance between the axes of two rods, given their centres, unit
  orientations and half lengths. closestApproachLineSegments degenerates for
  (anti)parallel rods, so the end caps are tested against the opposite axis too
  and the smallest separation is kept. Working on raw geometry (rather than on
  IBacterium*) means rejected trial positions do not burn ids from the static
  RodShapedBacterium counter.
*/
static double axisSeparation(
  const Vec3 &ri, const Vec3 &ni, double di,
  const Vec3 &rj, const Vec3 &nj, double dj
)
{
    double ti, tj, closest_sq;
    closestApproachLineSegments( ri, ni, di, rj, nj, dj, ti, tj, closest_sq );

    double dist_sq, tt;
    Vec3 cap;
    for ( const double end : { -1.0, 1.0 } )
    {
        cap = ri + ( end * di ) * ni;
        findSqrDistToLine( rj, nj, dj, cap, dist_sq, tt );
        closest_sq = std::min( closest_sq, dist_sq );

        cap = rj + ( end * dj ) * nj;
        findSqrDistToLine( ri, ni, di, cap, dist_sq, tt );
        closest_sq = std::min( closest_sq, dist_sq );
    }
    return std::sqrt( closest_sq );
}

/*
  Fill a circular droplet with exactly numTypeA Candida and numTypeB PA cells.

  The loop is driven by the target population, never by a trial budget: every
  planned cell is placed before this returns, so the composition of the initial
  condition is exactly what was asked for. Positions are drawn uniformly over
  the droplet and accepted only if the whole rod clears every cell already down.

  To make that terminate, a hexagonal set of reserved sites is built up front at
  a pitch wide enough that two cells on distinct sites cannot touch whatever
  their orientations. A random position is accepted only while enough sites stay
  free for every cell still queued, and a cell that cannot find a random spot is
  dropped onto a free site instead. The fallback therefore can never run out.
  If the trap genuinely cannot hold the population we say so and stop, rather
  than silently returning a partly filled droplet.
*/
std::vector<IBacterium*> initialiseBiofilm(double linking1, double linking2, 
  int numTypeA, int numTypeB, 
  double centerX, double centerY) {  // One shared center
    // Set the initial conditions
    std::vector<IBacterium*> initial_conditions;

    // Define wall constraints scaled by the diameter of PA Bacteria
    const double xMin { -73.0 / 1.17 };
    const double xMax {  73.0 / 1.17 };
    const double yMax {  147.0 / 1.17 };
    const double yMin { 3.0 };

    // Surface-to-surface clearance demanded between two cells at t=0. Anything
    // less and the pair starts overlapping, which fires a large Hertzian force
    // on the first step.
    const double minGap { 0.3 };

    // Keep the droplet this far clear of the trap walls
    const double dropletMargin { 4.0 };

    struct InitialCellData {
        bool isCandida;
        double length;
        double linkingProb;
        double radius;
    };

    std::vector<InitialCellData> plannedCells;
    plannedCells.reserve(numTypeA + numTypeB);

    for (int ii = 0; ii < numTypeA; ++ii) {
        plannedCells.push_back(InitialCellData{true, 4.0, 1.0, Candida::mRadius});
    }
    for (int ii = 0; ii < numTypeB; ++ii) {
        plannedCells.push_back(InitialCellData{false, 3.0, 0.0, 0.5});
    }

    const int targetCount { numTypeA + numTypeB };
    if ( targetCount <= 0 ) return initial_conditions;

    // Half the tip-to-tip extent of the largest cell: a centre at least this far
    // inside a wall keeps the whole cell in the box whatever its orientation.
    double maxExtent { 0.0 };
    for ( const auto& cellData : plannedCells ) {
        maxExtent = std::max( maxExtent, 0.5*cellData.length + cellData.radius );
    }

    //--------------------------------------------------------
    //-----uniformly distribute cells/segments in a droplet----
    //--------------------------------------------------------

    // Distance from the droplet centre to the nearest wall. Note this is
    // measured from where the droplet actually sits, not from the half width of
    // the box, so the droplet stays inside the trap if the centre is moved.
    const double wallClearance {
        std::min( std::min( centerX - xMin, xMax - centerX ),
                  std::min( centerY - yMin, yMax - centerY ) )
    };
    const double maxDropletRadius { wallClearance - maxExtent };
    if ( maxDropletRadius <= 0.0 ) {
        std::cerr << "Error: droplet centre (" << centerX << ',' << centerY
                  << ") is too close to the trap walls to hold a single cell.\n";
        std::exit(EXIT_FAILURE);
    }
    double dropletRadius {
        std::min( wallClearance - dropletMargin, maxDropletRadius )
    };

    // Two cells whose centres are at least sitePitch apart are clear of each
    // other by at least minGap for every pair of orientations.
    const double sitePitch { 2.0*maxExtent + minGap };
    // Compared against squared distances; the slack keeps exact lattice
    // neighbours, which sit at exactly sitePitch, from being discarded.
    const double reserveDistSq { sitePitch * sitePitch * (1.0 - 1e-9) };

    auto buildSites = [&](double radius) {
        std::vector<Vec3> sites;
        const double rowHeight { sitePitch * std::sqrt(3.0) * 0.5 };
        const int numRows { static_cast<int>( std::ceil(radius/rowHeight) ) };
        const int numCols { static_cast<int>( std::ceil(radius/sitePitch) ) + 1 };
        for ( int row = -numRows; row <= numRows; ++row ) {
            const double yy { centerY + row*rowHeight };
            const double rowOffset { (row % 2 == 0) ? 0.0 : 0.5*sitePitch };
            for ( int col = -numCols; col <= numCols; ++col ) {
                const double xx { centerX + col*sitePitch + rowOffset };
                if ( std::hypot(xx - centerX, yy - centerY) <= radius ) {
                    sites.push_back( Vec3{xx, yy, 0.0} );
                }
            }
        }
        return sites;
    };

    std::vector<Vec3> freeSites { buildSites(dropletRadius) };

    // Grow the droplet, never past the walls, until there is a reserved site for
    // every cell. Only then is the full population guaranteed to be deliverable.
    while ( static_cast<int>(freeSites.size()) < targetCount &&
            dropletRadius < maxDropletRadius )
    {
        dropletRadius = std::min( dropletRadius + sitePitch, maxDropletRadius );
        freeSites = buildSites(dropletRadius);
    }

    if ( static_cast<int>(freeSites.size()) < targetCount ) {
        std::cerr << "Error: cannot guarantee " << targetCount
                  << " cells in the trap. The largest droplet centred on ("
                  << centerX << ", " << centerY << ") has radius " << dropletRadius
                  << ", which holds " << freeSites.size() << " cells at the "
                  << sitePitch << " pitch needed to keep any two cells "
                  << minGap << " apart whatever their orientation. Ask for fewer "
                  << "cells, move the droplet centre, or lower minGap.\n";
        std::exit(EXIT_FAILURE);
    }

    // Create a global random number generator with a unique seed
    std::random_device rd;
    std::mt19937 gen(rd() ^ std::mt19937::result_type(std::time(0)) ^ std::mt19937::result_type(getpid()));

    std::shuffle(plannedCells.begin(), plannedCells.end(), gen);
    std::shuffle(freeSites.begin(), freeSites.end(), gen);

    std::uniform_real_distribution<double> distAngle(0, 2 * constants::pi);
    std::uniform_real_distribution<double> distUnit(0.0, 1.0);

    struct PlacedCell {
        Vec3 pos;
        Vec3 orientation;
        double halfLength;
        double radius;
    };
    std::vector<PlacedCell> placed;
    placed.reserve(targetCount);

    const int attemptsPerCell { 5000 };
    int numFallback { 0 };
    int numPlacedA { 0 }, numPlacedB { 0 };

    for ( int ii = 0; ii < targetCount; ++ii )
    {
        const InitialCellData& cellData { plannedCells[ii] };
        const double halfLength { 0.5 * cellData.length };
        const int cellsQueued { targetCount - ii - 1 }; // still waiting after this one

        double angle { 0.0 };
        Vec3 pos, orientation;
        bool placedRandomly { false };

        for ( int attempt = 0; attempt < attemptsPerCell && !placedRandomly; ++attempt )
        {
            // Sample uniformly over the droplet area. Because dropletRadius is
            // at most wallClearance-maxExtent, any centre drawn here keeps the
            // whole cell inside the box for every orientation, so no separate
            // end-pole boundary test is needed.
            const double radialDistance { dropletRadius * std::sqrt(distUnit(gen)) };
            const double polarAngle { distAngle(gen) };
            angle = distAngle(gen);

            pos = Vec3{ centerX + radialDistance * std::cos(polarAngle),
                        centerY + radialDistance * std::sin(polarAngle), 0.0 };
            orientation = Vec3{ std::cos(angle), std::sin(angle), 0.0 };

            // True rod-rod clearance, so long cells cannot interpenetrate the
            // way a centre-to-centre test lets them.
            bool clear { true };
            for ( const auto& other : placed )
            {
                const double surfaceGap {
                    axisSeparation( pos, orientation, halfLength,
                                    other.pos, other.orientation, other.halfLength )
                    - cellData.radius - other.radius
                };
                if ( surfaceGap < minGap ) { clear = false; break; }
            }
            if ( !clear ) continue;

            // Accept only while every cell still queued keeps a reserved site.
            int sitesLeft { 0 };
            for ( const auto& site : freeSites )
            {
                const double ddx { site.x - pos.x }, ddy { site.y - pos.y };
                if ( ddx*ddx + ddy*ddy >= reserveDistSq && ++sitesLeft >= cellsQueued ) break;
            }
            if ( sitesLeft < cellsQueued ) continue;

            placedRandomly = true;
        }

        if ( !placedRandomly )
        {
            // Fall back to a reserved site: free means it is at least sitePitch
            // from every cell already placed, so it is clear by construction.
            pos = freeSites.front();
            angle = distAngle(gen);
            orientation = Vec3{ std::cos(angle), std::sin(angle), 0.0 };
            ++numFallback;
        }

        // Sites this cell now covers can no longer host one of the queued cells
        freeSites.erase(
            std::remove_if( freeSites.begin(), freeSites.end(),
                [&](const Vec3 &site) {
                    const double ddx { site.x - pos.x }, ddy { site.y - pos.y };
                    return ddx*ddx + ddy*ddy < reserveDistSq;
                } ),
            freeSites.end()
        );

        placed.push_back( PlacedCell{ pos, orientation, halfLength, cellData.radius } );

        auto* rod = new RodShapedBacterium{
            pos.x, pos.y, 0,  // position x, y, z (in 2D z=0)
            angle,            // random angle
            constants::pi * 0.5,
            RodShapedBacterium::mAvgGrwthRate,
            cellData.length,
            cellData.linkingProb,
            cellData.radius
        };
        initial_conditions.push_back(rod);
        if ( cellData.isCandida ) ++numPlacedA; else ++numPlacedB;
    }

    // The population is the target by construction; check it anyway so a future
    // change to the loop cannot quietly go back to under-filling the droplet.
    assert( static_cast<int>(initial_conditions.size()) == targetCount );
    assert( numPlacedA == numTypeA && numPlacedB == numTypeB );

    std::cout << "droplet centre\t: (" << centerX << ", " << centerY << ")\n"
              << "droplet radius\t: " << dropletRadius
              << " (max that fits " << maxDropletRadius << ")\n"
              << "placed\t\t: " << initial_conditions.size() << '/' << targetCount
              << " cells (" << numPlacedA << " Candida, " << numPlacedB << " PA), "
              << numFallback << " on reserved sites\n";

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
