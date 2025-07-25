// Standard libraries
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <cassert>
#include <random>
#include <chrono>
#include <string>
#include <memory>
#include <numeric>
#include <filesystem>

// Custom classes
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"
#include "IBacterium.hpp"
#include "RodShapedBacteria.hpp"
#include "RandUtil.hpp"
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "IO.hpp"
#include "PolyBiofilm.hpp"


void PolyBiofilm::createLogFile()
{
  std::filesystem::create_directories(sim_out_dir);
  std::string filename=sim_out_dir+'/'+"biofilm.log";
  std::ofstream out_file;
  out_file.open( filename, std::ios::out );
  if (!out_file)
  {
    std::cerr << "Failed to open "
              << filename
              << " for writing!"
              << std::endl;
    exit (10);
  }
  std::vector<std::string> headers{};
  std::vector<double> values{};

  /*=== save basics ===*/
  headers.push_back("OutputFrequency");
  values.push_back(mOutFreq);

  //headers.push_back("RodAspectRatio");
  //values.push_back(mAvgDivLen/(2*RodShapedBacterium::mRadius));

  headers.push_back("RodModE");
  values.push_back(RodShapedBacterium::mRodModE);

  // In microns / hr
  headers.push_back("RodGrowthRate");
  values.push_back(RodShapedBacterium::mAvgGrwthRate);

#if defined(CHAINING)
  headers.push_back("Kappa");
  values.push_back(RodShapedBacterium::mKappa);

  headers.push_back("BendRig");
  values.push_back(RodShapedBacterium::mBendRig);


  headers.push_back("ForceThresh");
  values.push_back(RodShapedBacterium::mForceThresh);
#endif

  for (const auto & header : headers)
  {
    out_file << header << '\t';
  }
  out_file << "SEED" << '\n';

  for (double value : values)
  {
    out_file << value << '\t';
  }
  out_file << gen_rand.getSeed() << '\n';
  out_file.close();
  std::cout << "Saved log file to " << filename << '\n';
}


long mTotalDivisions = 0; // Initialize the global division counter
void PolyBiofilm::updateOneTimeStep(bool &update_neighbours, uint &verlet_counter)
{
  int divisionCount = 0;
  // ---------------------- Grow and divide cells ------------------------------
  for ( int ii=mCells.size()-1; ii>=0; --ii )
  {
    {
      mCells[ii]->grow(mDt);
      if ( mCells[ii]->signalDivide() )
      {
        mCells[ii]->divide(mCells);
        // Increment division counter
        divisionCount++;
        // Any division will require a re-binning of the cells
        update_neighbours=true;
      }
    }
    // Update the total number of divisions
    mTotalDivisions += divisionCount;

    // Optional debug output:
    // std::cerr << "Divisions in this step: " << divisionCount << "\n";
  }

  // ---------------------- Generate Verlet List ------------------------------
  // cell lists need to be updated
  if ( update_neighbours )
  {
    mGrid.updateVerletList(mCells);
    verlet_counter=0;
    update_neighbours=false;
  }

  // Update the force and torque for each cell to the current time
  polyInteractParticles(mCells);

  for ( auto &cell : mCells )
  {
    cell->move(mDt);
    double moved {
      dot2( cell->getLoggedPos()-cell->getPos() )
    };
    if ( moved>dot2( mGrid.mVerletMoveR ) )
    {
      update_neighbours=true;
    }
    cell->reset();
  }

}

void PolyBiofilm::runSim()
{
  createLogFile();

  std::cout << "------------------------------------------------------" << '\n';
  std::cout << "------------------Begin Simulation--------------------" << '\n';
  std::cout << "------------------------------------------------------" << '\n';

  uint verlet_counter { 0 };      // Count successful steps since last rebinning
  long output_counter { 0 };      // Index the outputs
  bool update_neighbours { true };// Always bin neighbours on first step
  const ulong num_outputs { static_cast<ulong>(ceil(mTimeSteps/mOutFreq)) };
  for ( ulong tt=0; tt<=mTimeSteps; ++tt )
  {
    // check for outputting
    if ( tt%mOutFreq==0 )
    {
      updateStatus(output_counter,num_outputs);

      //Output to file
      printCellsToFile(
        createFileName(
          output_counter,
          "biofilm_",
          sim_out_dir
        ),
        mCells,
        false,
        false
      );

      ++output_counter;

      if ( tt>=mTimeSteps )
      {
        break;
      }
    }

    // refresh list if every N timesteps
    if ( verlet_counter>=mGrid.mVerletUptime )
    {
      update_neighbours=true;
    }
    else ++verlet_counter;

    updateOneTimeStep(update_neighbours,verlet_counter);

    if ( getLength()>=mTargetSize  )
    {
      printCellsToFile(
        createFileName(
          output_counter,
          "final_",
          sim_out_dir
        ),
        mCells,
        false,
        false);
        break;
    }
  }
  std::cout << '\n';
      // Print the total number of divisions
  std::cout << "Total divisions during simulation: " << mTotalDivisions << '\n';
}
