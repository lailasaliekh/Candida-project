/*******************************************************************************
    BiofilmDES  - a program that simulates a growing colony of microbial cells

    Contributing author:
    Rory Claydon, University of Edinburgh, rory.claydon@ed.ac.uk

    Copyright (2020) The University of Edinburgh.

    The software is based on algorithms described in:

    Mechanically driven growth of quasi-two dimensional microbial colonies,
    F.D.C. Farrell, O. Hallatschek, D. Marenduzzo, B. Waclaw,
    Phys. Rev. Lett. 111, 168101 (2013).

    Three-dimensional distinct element simulation of spherocylinder crystallization.
    Pournin, L., Weber, M., Tsukahara, M. et al.
    Granul. Matter 7, 119–126 (2005).

    A fast algorithm to evaluate the shortest distance between rods,
    C. Vega, S. Lago,
    Comput. Chem., 18(1), 55-59 (1994)

    I would like to thank Bartlomiej Waclaw from Edinburgh University for some
    very useful discussions on algorithm stability, timestep choice and some
    potential optimisations to try out in future.

    This file is part of BiofilmDES.

    BiofilmDES is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BiofilmDES is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file
    License.txt or at <http://www.gnu.org/licenses/>.

    Compilation and run from current directory:
      make && ./biofilm.out 0 1.1 0.95

    Further details in the documentation

*******************************************************************************/

// Standard libraries
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <array>
#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include <random>
#include <memory>
#include <omp.h>
#include <unordered_set>
// User defined
#include "VerletGrid.hpp"
#include "forces.hpp"
#include "RodShapedBacteria.hpp"
#include "Candida.hpp"
#include "constants.hpp"         // definition of constants namespace
#include "MathUtility.hpp"       // definition of Vec3 class

//*----------------------------------------------------------------------------------------------------
//this is to insure we are deleting cells and if they are in a chain the links being deleted properly/////
void removeCellsOutsideTrap(std::vector<IBacterium*> &cells, double trapYLimit) {
  auto isOutsideTrap = [trapYLimit](IBacterium* cell) -> bool {
      std::array<Vec3, 2> poles;
      cell->getMyEndVecs(poles[0], poles[1]); // Get the end points (poles) of the bacterium
      
      // Check if both poles are outside the y < trapYLimit
      return poles[0].y < trapYLimit && poles[1].y < trapYLimit;
  };

  // Identify cells to be removed
  std::vector<IBacterium*> cellsToRemove;
  for (auto& cell : cells) {
      if (isOutsideTrap(cell)) {
          cellsToRemove.push_back(cell);
      }
  }

  // Update links and clean up before deleting cells
  for (auto& cell : cellsToRemove) {
      // Get links of the current cell
      IBacterium* lowerLink = cell->getLowerLink();
      IBacterium* upperLink = cell->getUpperLink();

      // Handle the links of neighboring cells
      if (lowerLink) {
          lowerLink->setUpperLink(nullptr); // Disconnect lower neighbor's upper link
      }
      if (upperLink) {
          upperLink->setLowerLink(nullptr); // Disconnect upper neighbor's lower link
      }

      // Print details of the removed cell for debugging
      std::array<Vec3, 2> poles;
      cell->getMyEndVecs(poles[0], poles[1]);
      std::cout << "Deleting cell with ID: " << cell->getID()
                << ", Position: " << cell->getPos()
                << ", Poles: (" << poles[0] << ") and (" << poles[1] << ")" 
                << std::endl;
  }

  // Remove cells completely outside the trap
  cells.erase(std::remove_if(cells.begin(), cells.end(), isOutsideTrap), cells.end());
}
//!!!!!------------- remove cells when overlapping to prevent sims from crashing--------
void removeOverlappingCells(std::vector<IBacterium*> &cells, double sepThreshold) {
    std::unordered_set<IBacterium*> cellsToRemove;

    // Use neighbor lists instead of full pairwise check
    for (IBacterium* cell : cells) {
        for (IBacterium* neighbor : cell->getNeighbourList()) {  // Only check nearby cells
            double sep;
            Vec3 cA, cB;

            getMinDist(cell, neighbor, sep, cA, cB);  // Compute separation using the virtual contact points
            Vec3 cv = cA - cB;  // Compute separation vector
            sep = cv.norm();  // Compute separation distance

            if (sep <= sepThreshold) {
                // Choose only one cell to remove (based on length or ID)
                IBacterium* toDelete = (cell->getLength() < neighbor->getLength()) ? cell : neighbor;
                cellsToRemove.insert(toDelete);
            }
        }
    }

    // Ensure chain integrity before removing cells
    for (auto& cell : cellsToRemove) {
        IBacterium* lowerLink = cell->getLowerLink();
        IBacterium* upperLink = cell->getUpperLink();

        if (lowerLink) lowerLink->setUpperLink(nullptr);
        if (upperLink) upperLink->setLowerLink(nullptr);
    }

    // Remove selected cells (O(n) deletion instead of O(n²))
    cells.erase(std::remove_if(cells.begin(), cells.end(),
                               [&cellsToRemove](IBacterium* cell) {
                                   return cellsToRemove.count(cell) > 0;
                               }),
                cells.end());
}


// void removeOverlappingCells(std::vector<IBacterium*> &cells, double sepThreshold) {
//   std::vector<IBacterium*> cellsToRemove;

//   for (size_t i = 0; i < cells.size(); ++i) {
//       for (size_t j = i + 1; j < cells.size(); ++j) {
//           double sep;  // Separation distance
//           Vec3 cA, cB;  // Virtual contact points

//           getMinDist(cells[i], cells[j], sep, cA, cB);  // Compute separation

//           Vec3 cv = cA - cB;  // Compute separation vector
//           sep = cv.norm();  // Compute separation distance

//           if (sep <= sepThreshold) {
//               // Determine which cell to remove based on length, ID, or another property
//               IBacterium* toDelete;
//               IBacterium* toKeep;

//               if (cells[i]->getLength() < cells[j]->getLength()) {
//                   toDelete = cells[i];
//                   toKeep = cells[j];
//               } else if (cells[i]->getLength() > cells[j]->getLength()) {
//                   toDelete = cells[j];
//                   toKeep = cells[i];
//               } else {
//                   // If lengths are equal, remove the cell with the higher ID
//                   toDelete = (cells[i]->getID() > cells[j]->getID()) ? cells[i] : cells[j];
//                   toKeep = (cells[i]->getID() > cells[j]->getID()) ? cells[j] : cells[i];
//               }

//               std::cerr << "Deleting overlapping cell ID: " << toDelete->getID()
//                         << " (Keeping ID: " << toKeep->getID() << ")\n";

//               cellsToRemove.push_back(toDelete);
//           }
//       }
//   }

//   // Ensure chain integrity before removing cells
//   for (auto& cell : cellsToRemove) {
//       IBacterium* lowerLink = cell->getLowerLink();
//       IBacterium* upperLink = cell->getUpperLink();

//       if (lowerLink) lowerLink->setUpperLink(nullptr);
//       if (upperLink) upperLink->setLowerLink(nullptr);

//       std::cerr << "Deleted cell ID: " << cell->getID() << "\n";
//   }

//   // Remove the selected cells from the main container
//   cells.erase(std::remove_if(cells.begin(), cells.end(),
//                              [&cellsToRemove](IBacterium* cell) {
//                                  return std::find(cellsToRemove.begin(), cellsToRemove.end(), cell) != cellsToRemove.end();
//                              }),
//               cells.end());
// }
//-------------------------------
// void removeOverlappingCells(std::vector<IBacterium*> &cells, double sepThreshold) {
//   std::vector<IBacterium*> cellsToRemove;

//   // Identify overlapping cells
//   for (size_t i = 0; i < cells.size(); ++i) {
//       for (size_t j = i + 1; j < cells.size(); ++j) {
//           double sep;  // Separation distance
//           Vec3 cA, cB;  // Virtual contact points

//           getMinDist(cells[i], cells[j], sep, cA, cB);  // Correct function call

//           Vec3 cv = cA - cB;  // Compute separation vector
//           sep = cv.norm();  // Compute separation distance (ensures consistency)

//           if (sep <= sepThreshold) {
//               std::cerr << "Deleting overlapping cells (sep <= " << sepThreshold << "):\n";
//               std::cerr << "Cell A ID: " << cells[i]->getID() << " Contact Position: " << cA << "\n";
//               std::cerr << "Cell B ID: " << cells[j]->getID() << " Contact Position: " << cB << "\n";
//               std::cerr << "Separation vector: " << cv << "\n";
//               std::cerr << "Separation distance: " << sep << "\n";

//               cellsToRemove.push_back(cells[i]);
//               cellsToRemove.push_back(cells[j]);
//           }
//       }
//   }

//   // Ensure chain integrity before removing cells
//   for (auto& cell : cellsToRemove) {
//       IBacterium* lowerLink = cell->getLowerLink();
//       IBacterium* upperLink = cell->getUpperLink();

//       if (lowerLink) lowerLink->setUpperLink(nullptr);
//       if (upperLink) upperLink->setLowerLink(nullptr);

//       std::cerr << "Deleted cell ID: " << cell->getID() << "\n";
//   }

//   // Remove the overlapping cells from the main container
//   cells.erase(std::remove_if(cells.begin(), cells.end(),
//                              [&cellsToRemove](IBacterium* cell) {
//                                  return std::find(cellsToRemove.begin(), cellsToRemove.end(), cell) != cellsToRemove.end();
//                              }),
//               cells.end());
// }


//-----------------------------------------------!!>>>>><<<<


//*to use it for cells and walls closest contact point */
//-----------------------------------------------//
Vec3 closestPointOnRodToVerticalWall(const IBacterium* cell, double wallX) {
    Vec3 pos = cell->getPos();
    Vec3 orientation = cell->getOrientation();
    Vec3 rodStart = pos - (0.5 * cell->getLength()+cell->getRadius()) * orientation;
    Vec3 rodEnd = pos + (0.5 * cell->getLength()+cell->getRadius()) * orientation;

    // Project rod onto the x-axis to find closest point to wallX
    double t = (wallX - rodStart.x) / (rodEnd.x - rodStart.x);
    t = std::clamp(t, 0.0, 1.0); // Clamp to rod segment
    return rodStart + t * (rodEnd - rodStart);
}

Vec3 closestPointOnRodToTrapUpperWall(const IBacterium* cell, double wallY) {
    Vec3 pos = cell->getPos();
    Vec3 orientation = cell->getOrientation();
    Vec3 rodStart = pos - (0.5 * cell->getLength()+cell->getRadius()) * orientation;
    Vec3 rodEnd = pos + (0.5 * cell->getLength()+cell->getRadius()) * orientation;

    // Project rod onto the x-axis to find closest point to wallX
    double t = (wallY - rodStart.y) / (rodEnd.y - rodStart.y);
    t = std::clamp(t, 0.0, 1.0); // Clamp to rod segment
    return rodStart + t * (rodEnd - rodStart);
}
//-----------------------------------------------------//

void interactTrapWall(IBacterium *cell, double wallY, bool isLowerWall, double wallHertzianConstant, 
                      double trapStartX, double trapEndX, double trapHeight) {
  std::array<Vec3, 2> poles;
  cell->getMyEndVecs(poles[0], poles[1]);  // Get the end points (poles) of the rod-shaped bacterium

  const double Rstar = cell->getRadius();  // Cell's radius
  const double Estar = getEffectiveQ(cell->getModE(), constants::nondim_agar_mod_E);  // Effective modulus (from cell-surface interaction)

  for (int ii = 0; ii < 2; ++ii) {  // Loop over both poles (ends) of the rod-shaped bacterium
    double currentWallY = wallY;

    // Adjust the upper wall height if within the trap region
    if (!isLowerWall && poles[ii].x >= trapStartX && poles[ii].x <= trapEndX) {
      currentWallY += trapHeight;
    }


    //adding new logic for walls forces/torques
    // upper horizontal trap wall interaction
    //------------------------------------
    Vec3 closestPointLeft = closestPointOnRodToTrapUpperWall(cell, currentWallY);
    double upper_sep = fabs(closestPointLeft.y - currentWallY);
    double upper_overlap = std::max(cell->getRadius() - upper_sep, 0.0);
    if (upper_overlap > 0) {
        double force_mag = Estar * upper_overlap * sqrt(Rstar);
        Vec3 force(0, -force_mag, 0); // Push right (+x direction)
        cell->addForce(force);

        // Apply force at the surface of the cell (closest point adjusted by radius)
        Vec3 force_pos = closestPointLeft;
        Vec3 torque = cross(force_pos - cell->getPos(), force);
        cell->addTorque(torque);
        //---------------------------------------------
    }
    // Handle vertical walls (left at trapStartX, right at trapEndX)
    if (!isLowerWall) {
      // // Left vertical wall at trapStartX
      // double left_sep = fabs(poles[ii].x - trapStartX);  // Distance from left vertical wall
      // double left_overlap = std::max(cell->getRadius() - left_sep, 0.0);
      // if (left_overlap > 0 && poles[ii].y >= wallY) {  // Ensure it acts only above the upper wall
      //   double force_mag = Estar * left_overlap * sqrt(Rstar * left_overlap);
      //   //Vec3 force((poles[ii].x < trapStartX ? force_mag : -force_mag), 0, 0);  // Force direction in x-axis
      //   Vec3 force(force_mag, 0, 0); // Always push right
      //   cell->addForce(force);

      //   Vec3 force_pos = Vec3(trapStartX, poles[ii].y, 0);  // Position of left wall
      //   Vec3 torque = cross(force_pos - cell->getPos(), force);
      //   cell->addTorque(torque);
      // }

      // // Right vertical wall at trapEndX
      // double right_sep = fabs(poles[ii].x - trapEndX);  // Distance from right vertical wall
      // double right_overlap = std::max(cell->getRadius() - right_sep, 0.0);
      // if (right_overlap > 0 && poles[ii].y >= wallY) {  // Ensure it acts only above the upper wall
      //   double force_mag = Estar * right_overlap * sqrt(Rstar * right_overlap);
      //   Vec3 force(-force_mag, 0, 0); // Always push left
      //   //Vec3 force((poles[ii].x > trapEndX ? -force_mag : force_mag), 0, 0);  // Force direction in x-axis
      //   cell->addForce(force);

      //   Vec3 force_pos = Vec3(trapEndX, poles[ii].y, 0);  // Position of right wall
      //   Vec3 torque = cross(force_pos - cell->getPos(), force);
      //   cell->addTorque(torque);
      // }


      //--------------- new logic-----------------
      // Left vertical wall interaction
      Vec3 closestPointLeft = closestPointOnRodToVerticalWall(cell, trapStartX);
      double left_sep = fabs(closestPointLeft.x - trapStartX);
      double left_overlap = std::max(cell->getRadius() - left_sep, 0.0);
      if (left_overlap > 0) {
          double force_mag = Estar * left_overlap * sqrt(Rstar);
          Vec3 force(force_mag, 0, 0); // Push right (+x direction)
          cell->addForce(force);

          // Apply force at the surface of the cell (closest point adjusted by radius)
          Vec3 force_pos = closestPointLeft;
          Vec3 torque = cross(force_pos - cell->getPos(), force);
          cell->addTorque(torque);
      }

      // Right vertical wall interaction similarly
      Vec3 closestPointRight = closestPointOnRodToVerticalWall(cell, trapEndX);
      double right_sep = fabs(closestPointRight.x - trapEndX);
      double right_overlap = std::max(cell->getRadius() - right_sep, 0.0);
      if (right_overlap > 0) {
          double force_mag = Estar * right_overlap * sqrt(Rstar);
          Vec3 force(-force_mag, 0, 0); // Push left (-x direction)
          cell->addForce(force);

          Vec3 force_pos = closestPointRight;
          Vec3 torque = cross(force_pos - cell->getPos(), force);
          cell->addTorque(torque);
      }
    }
  }
}
// this is for the channel confinement only//
//**-----------------------------------------------------**/
void interactWall(IBacterium *cell, double wallY, bool isLowerWall, double wallHertzianConstant) {
  std::array<Vec3, 2> poles;
  cell->getMyEndVecs(poles[0], poles[1]);  // Get the end points (poles) of the rod-shaped bacterium

  const double Rstar = cell->getRadius();  // Cell's radius
  const double Estar = getEffectiveQ(cell->getModE(), constants::nondim_agar_mod_E);  // Effective modulus (from cell-surface interaction)

  for (int ii = 0; ii < 2; ++ii) {  // Loop over both poles (ends) of the rod-shaped bacterium
    double sep = fabs(poles[ii].y - wallY);  // Distance from the wall in the y-direction

    double overlap = std::max(cell->getRadius() - sep, 0.0);  // Calculate overlap
    if (overlap <= 0) continue;  // No overlap, no force

    double force_mag = Estar* overlap * sqrt(Rstar * overlap);  // Weaker Hertzian repulsive force

    Vec3 force(0, (isLowerWall ? force_mag : -force_mag), 0);  // Positive for lower wall, negative for upper wall
    cell->addForce(force);

    Vec3 force_pos = Vec3(poles[ii].x, wallY, 0);  // Position where force is applied (at the wall)
    Vec3 torque = cross(force_pos - cell->getPos(), force);  // Calculate torque based on applied force
    cell->addTorque(torque);
  }
}
//*----------------------------------------------------*/
//stacking walls so that the cells don't break through them/
void interactMultiWallWithTrap(IBacterium *cell, double yMin, double yMax, int layers, double trapStartX, double trapEndX, double trapHeight) {
  double wallHertzianConstant = 100.0;  // Weaker constant for each wall layer
  double totalForceMultiplier = 10.0 / layers;  // Split the force among all layers

  // // Interact with the lower wall (stacked walls at yMin with no gaps)
  // for (int i = 0; i < layers; ++i) {
  //   interactWall(cell, yMin, true, wallHertzianConstant * totalForceMultiplier);
  // }

  // Interact with the upper wall with a trap
  for (int i = 0; i < layers; ++i) {
    interactTrapWall(cell, yMax, false, wallHertzianConstant * totalForceMultiplier, trapStartX, trapEndX, trapHeight);
  }
}
#include <unordered_set>
#include <iostream>
#include <sstream>
#include <fstream>  // Include for file operations

void polyInteractParticles(std::vector<IBacterium*> &pars) {
    double yMin = -100000.0;  // Minimum y-boundary
    double yMax = 0.0;  // Maximum y-boundary
    int wallLayers = 10;  // Number of wall layers on each side
    double trapStartX = -75.0 / 1.17;  // Start of the trap region
    double trapEndX = 75.0 / 1.17;     // End of the trap region
    double trapHeight = 150.0 / 1.17;  // Height of the trap from the original upper wall


    // Remove trapped cells before interactions
    removeCellsOutsideTrap(pars, yMax);

    // Use an unordered set to track cells marked for deletion
    std::unordered_set<IBacterium*> cellsToRemove;

    #pragma omp parallel for shared(pars, yMin, yMax, wallLayers, trapStartX, trapEndX, trapHeight, cellsToRemove) \
            schedule(static) default(none)
    for (uint ii = 0; ii < pars.size(); ++ii) {
        IBacterium* cell = pars[ii];

        // Cell-to-cell interactions
        for (auto &neighbour_cell : cell->getNeighbourList()) {
            const bool different_cells {
                (cell->getID() != neighbour_cell->getID()) ||
                (cell->getMyType() != neighbour_cell->getMyType())
            };
            if (different_cells) {
                pairCollision(cell, neighbour_cell);
            }
        }

        // Interact with overlapping stacked walls, including the trap region for the upper wall
        interactMultiWallWithTrap(cell, yMin, yMax, wallLayers, trapStartX, trapEndX, trapHeight);

    #ifdef CHAINING
        computeChainingInteractions(cell);
    #endif

    #ifndef MOVE_3D
        interactSurface(cell);
        interactInterface(cell);
    #endif

        // Update cell velocity and angular velocity
        cell->setVel();
        cell->setAngVel();
    }

}

inline double computeHertzianForceMag(
  const IBacterium *A, const IBacterium *B, double sep, bool useLennardJones
)
{
  double sigma { A->getRadius() + B->getRadius() };
  double overlap { std::max(sigma-sep,0.0) }; // Cell overlap

  if ( sep>=pow(2,1/6)*sigma  ) return 0.0; // this is for the interaction cut off 

  // Lennard-Jones and Hertzian parameters
  constexpr double epsilon = 1e-6;  // strength of the interaction
  double sigma_r=sigma / sep;
  if (useLennardJones)
  {
    
    //wca force magnitude
    double lj_force = 24.0 * epsilon * 
                      (2 * pow(sigma_r, 13)/sigma- pow(sigma_r, 7)/sigma);
    return lj_force;
  }
  else
  {
    // Hertzian force magnitude
    // double overlap = std::max(sigma - sep, 0.0);  // Cell overlap
    // if (overlap <= 0)
    //     return 0.0;
    double hertzian_force = overlap * sqrt(overlap);
    return hertzian_force;
  }
}

inline void computeHertzianForceTorque(
  const IBacterium *A,
  const IBacterium *B,
  double &sep,
  Vec3 &cA,
  Vec3 &cB,
  Vec3 &force_pos,
  Vec3 &force,
  Vec3 &torque, 
  bool useLennardJones
)
{
  std::cout << "Update force with parallel force" << '\n';
  exit(12);
  getMinDist(A,B,sep,cA,cB);  // Get the minimum distance between the bacteria
  // const double min_sep=1e-3;
  // sep=std::max(sep,min_sep);
  assert(sep>0);
  // assert(sep>0.1);

  const Vec3 normal_BA { (cA-cB)/sep }; // Normal direction from B to A

  const double steric_force { computeHertzianForceMag(A,B,sep,useLennardJones) };
  force = steric_force*normal_BA;

  // Torque 0 is due to force of B on A
  // lever arm is vector from rod center to virt_sphere center on that rod
  // Take force at the surface of the cell
  force_pos = cA - A->getRadius()*normal_BA;
  torque = cross(force_pos-A->getPos(),force);
}

//inline int deletedCellCount = 0; // Tracks the number of deleted cells

inline void computeHertzianForceTorque(
  IBacterium *A,
  const IBacterium *B,
  Vec3 &cA,
  Vec3 &cB,
  bool useLennardJones
)
{
  Vec3 cv = cA-cB;
  double sep=cv.norm();
  // const double min_sep=1e-3;
  // sep=std::max(sep,min_sep);
  if (sep <= 1e-6) {
  // Print debug information
  std::cerr << "sep <= 0.1\n";
  std::cerr << "Details for A:\n";
  std::cerr << "  ID: " << A->getID() << "\n";
  std::cerr << "  Position: " << A->getPos() << "\n";
  std::cerr << "  Orientation: " << A->getOrientation() << "\n";
  std::cerr << "  Length: " << A->getLength() << "\n";
  std::cerr << "Details for B:\n";
  std::cerr << "  ID: " << B->getID() << "\n";
  std::cerr << "  Position: " << B->getPos() << "\n";
  std::cerr << "  Orientation: " << B->getOrientation() << "\n";
  std::cerr << "  Length: " << B->getLength() << "\n";
  std::cerr << "Separation vector (cA - cB): " << cv << "\n";
  std::cerr << "Separation distance: " << sep << "\n";
  }
  assert(sep > 0);
  Vec3 normal_BA { cv/sep }; // Normal direction from B to A

  double steric_force { computeHertzianForceMag(A,B,sep,useLennardJones) };
  Vec3 force = steric_force*normal_BA;

  // Torque 0 is due to force of B on A
  // lever arm is vector from rod center to virt_sphere center on that rod
  // Take force at the surface of the cell
  Vec3 force_pos = cA - A->getRadius()*normal_BA;
  Vec3 torque = cross(force_pos-A->getPos(),force);

  A->addForce(force);
  A->addTorque(torque);
}

inline void interactSurface(IBacterium *cell)
{
  // agar lowest and highest possible values
  constexpr double agar_low {
    constants::nondim_min_z - 0.5*constants::nondim_agar_roughness
  };
  constexpr double agar_high { agar_low+constants::nondim_agar_roughness };

  std::array<Vec3,2> poles;
  cell->getMyEndVecs(poles[0],poles[1]);

  // Reduced quantities for use in the Hertzian force calculation
  const double Rstar { cell->getRadius() };
  const double Estar {
    getEffectiveQ(cell->getModE(),constants::nondim_agar_mod_E)
  };

  const int num_poles { 1 + ( dot2(poles[1]-poles[0])!=0 ) };
  for ( int ii=0; ii<num_poles; ++ii )
  {
    // Height of the agar at the poles indexed by 1,2
    const double local_agar_height {
      gen_rand.getUniformRand(agar_low,agar_high)
    };

    const double sep { poles[ii].z-local_agar_height };

    // Cell pole overlap with agar
    const double overlap { std::max(cell->getRadius()-sep,0.0) };
    if ( overlap<=0 ) return;

    // Note: assumed that cells are essentialy incompressible
    // Rstar is just the cell radius here as the wall is an infinite half plane
    const double force_mag { Estar * overlap * sqrt(Rstar*overlap) };
    const Vec3 force { 0,0,force_mag };

    cell->addForce(force);
    const Vec3 force_pos { poles[ii].x, poles[ii].y, local_agar_height };
    const Vec3 torque = cross(force_pos-cell->getPos(),force);
    cell->addTorque(torque);
  }
}

inline void interactInterface(
  IBacterium *cell
)
{
  const Vec3 gravity { 0,0,-constants::effective_g };
  cell->addForce(gravity);
}

/* === Better parallel handling === */
void pairCollision(IBacterium *A, const IBacterium *B)
{

  double s1,t1,s2,t2; // Line parameter values at which they are closest
  Vec3 c1,c2;         // centres of the closet approach on each line resp.

  // COMs
  Vec3 pos_A { A->getPos() };
  Vec3 pos_B { B->getPos() };

  // find the direction Segment of S1 and S2
  Vec3 v1 = 0.5*A->getLength()*A->getOrientation();
  Vec3 v2 = 0.5*B->getLength()*B->getOrientation();

  // Check line was parallel
  bool par = closestApproachLineSegmentsParallel(pos_A,v1,pos_B,v2,s1,t1,s2,t2);

  c1 = pos_A+v1*s1; // centre of the virtual sphere on cell A
  c2 = pos_B+v2*t1; // centre of the virtual sphere on cell B
  computeHertzianForceTorque(A,B,c1,c2,true);

  // Parallel segements so two points were returned
  if ( par )
  {
    c1 = pos_A+v1*s2; // centre of the virtual sphere on cell A
    c2 = pos_B+v2*t2; // centre of the virtual sphere on cell B
    computeHertzianForceTorque(A,B,c1,c2,true);
  }
}

inline double computeHertzianEnergyMag(
  IBacterium *A,
  const IBacterium *B,
  Vec3 cA,
  Vec3 cB, bool useLennardJones
)
{
  Vec3 cv = cA-cB;
  double sep=cv.norm();
  assert(sep>0);

  // assert(sep>0.1);
  double sigma { A->getRadius() + B->getRadius() };
  if ( sep>=pow(2,1/6)*sigma  ) return 0.0;
  // double overlap { std::max(sigma-sep,0.0) }; // Cell overlap
  // double energy { (2.0/5.0)*overlap*overlap*sqrt(overlap) };
  // return energy;
  constexpr double epsilon = 1e-6;  // Pre-computed epsilon for matching potentials
  double sigma_r=sigma / sep;
  double overlap = std::max(sigma - sep, 0.0);  // Cell overlap
  if (useLennardJones)
  {
      // Compute WCA WITH THE ATTRACTIVE PART
      double lj_energy = 4 * epsilon * (pow(sigma_r, 12) - pow(sigma_r, 6)) + epsilon;
      return lj_energy;
  }
  else
  {
      // Compute Hertzian energy
      
      double hertzian_energy = (2.0 / 5.0) * overlap * overlap * sqrt(overlap);
      return hertzian_energy;
  }
}

/* === Better parallel handling === */
double getPairHertzianEnergy(IBacterium *A, const IBacterium *B)
{

  double energy { 0.0 }; // Energy due to an overlap
  double s1,t1,s2,t2;    // Line parameter values at which they are closest
  Vec3 c1,c2;            // centres of the closet approach on each line resp.

  // COMs
  Vec3 pos_A { A->getPos() };
  Vec3 pos_B { B->getPos() };

  // find the direction Segment of S1 and S2
  Vec3 v1 = 0.5*A->getLength()*A->getOrientation();
  Vec3 v2 = 0.5*B->getLength()*B->getOrientation();

  // Check line was parallel
  bool par = closestApproachLineSegmentsParallel(pos_A,v1,pos_B,v2,s1,t1,s2,t2);

  c1 = pos_A+v1*s1; // centre of the virtual sphere on cell A
  c2 = pos_B+v2*t1; // centre of the virtual sphere on cell B
  energy+=computeHertzianEnergyMag(A,B,c1,c2,true);

  // Parallel segements so two points were returned
  if ( par )
  {
    c1 = pos_A+v1*s2; // centre of the virtual sphere on cell A
    c2 = pos_B+v2*t2; // centre of the virtual sphere on cell B
    energy+=computeHertzianEnergyMag(A,B,c1,c2,true);
  }
  return energy;
}

/* === Chaining Bacteria Definitions === */
#ifdef CHAINING
inline Vec3 calcBiNormal(const Vec3 &t0, const Vec3 &t1)
{
  return 2 * cross(t0,t1) / ( t0.norm()*t1.norm() + dot(t0,t1) );
}

inline Vec3 calckb1t0(const Vec3 &t0, const Vec3 &t1)
{
  const Vec3 kb_1 { calcBiNormal(t0,t1) };
  return {
    ( 2 * cross(t1,kb_1) - dot(kb_1,kb_1)*( t1 + t0*(t1.norm()/t0.norm()) ) )
      / ( t0.norm()*t1.norm() + dot(t0,t1) )
  };
}

inline Vec3 calckb1t1(const Vec3 &t0, const Vec3 &t1)
{
  const Vec3 kb_1 { calcBiNormal(t0,t1) };
  return {
    ( -2 * cross(t0,kb_1) - dot(kb_1,kb_1)*( t0 + t1*(t0.norm()/t1.norm()) ) )
      / ( t0.norm()*t1.norm() + dot(t0,t1) )
  };
}

inline void computeSingleChainInteraction(
  const IBacterium* lower_cell,
  const IBacterium* upper_cell,
  Vec3 &lower_force,
  Vec3 &upper_force,
  Vec3 &lower_torque,
  Vec3 &upper_torque
)
{
  assert(upper_cell->getLowerLink()->getID()==lower_cell->getID());
  assert(lower_cell->getUpperLink()->getID()==upper_cell->getID());

  // constexpr double radius { 0.5*constants::nondim_rodSpheroDiam };
  const Vec3 lower_cell_n { lower_cell->getOrientation() };
  const Vec3 upper_cell_n { upper_cell->getOrientation() };

  // The top of the lower cell's head to which the bottom of the spring will attach
  const Vec3 lower_head {
    lower_cell->getPos()
    + lower_cell_n*(
        // radius
        + 0.5*lower_cell->getLength()
      )
  };

  // The tail of the uppers cell's head to which the top of the spring will attach
  const Vec3 upper_tail {
    upper_cell->getPos()
    - upper_cell_n*(
        // radius
        + 0.5*upper_cell->getLength()
      )
  };

  const Vec3 low_to_high { upper_tail-lower_head };
  const double mod_l_to_h { low_to_high.norm() };

  const double inv_mod_l_to_h { 1.0 / mod_l_to_h };
  const Vec3 ux12 { low_to_high * inv_mod_l_to_h };


  /*-------------------*/
  /*   Stiff linking   */
  /*-------------------*/
  const Vec3 lower_rod_base {
    lower_head - lower_cell_n*constants::nondim_candidaRadius*2
  };
  const Vec3 upper_rod_base {
    upper_tail + upper_cell_n*constants::nondim_candidaRadius*2
  };

  const Vec3 t0 { lower_head - lower_rod_base };

  // Check there is no performance penalty for this
  const Vec3 t1 { low_to_high };

  const Vec3 t2 { upper_rod_base - upper_tail };

  // define inverse lengths
  const double ds1_inv { 2/( t0.norm() + t1.norm() ) };
  const double ds2_inv { 2/( t1.norm() + t2.norm() ) };

  // Alias for bending rigidity
  const double K { RodShapedBacterium::mBendRig };
  const Vec3 stiff_force0 {
    0.5*K* ( calckb1t0(t0,t1) * ds1_inv )
  };
  const Vec3 stiff_force1 {
    -0.5*K* (
      ( calckb1t0(t0,t1) - calckb1t1(t0,t1) ) * ds1_inv - calckb1t0(t1,t2) * ds2_inv
    )
  };
  const Vec3 stiff_force2 {
    -0.5*K* (
      calckb1t1(t0,t1) * ds1_inv + ( calckb1t0(t1,t2) - calckb1t1(t1,t2)) * ds2_inv
    )
  };
  const Vec3 stiff_force3 {
   -0.5*K* ( calckb1t1(t1,t2) * ds2_inv )
 };

  // Both the forces here are conservative, hence forces on 2 is minus total on 1
  lower_force = stiff_force0 + stiff_force1;
  // upper_force = -lower_force; // Calculate this below after adding spring force

  lower_torque  = cross( lower_rod_base - lower_cell->getPos() , stiff_force0 );
  lower_torque += cross( lower_head     - lower_cell->getPos() , stiff_force1 );
  upper_torque  = cross( upper_tail     - upper_cell->getPos() , stiff_force2 );
  upper_torque += cross( upper_rod_base - upper_cell->getPos() , stiff_force3 );

 //  /*-------------------*/
 //  /*   Spring linking  */
 //  /*-------------------*/
 //
  const Vec3 spring_force { // we are letting the chained cells to overlap slightly so that the smaller PA cells don't cross the chains
    RodShapedBacterium::mKappa * ( mod_l_to_h - constants::nondim_candidaRadius*1.3 ) * ux12 //constants::nondim_candidaRadius*2
  };

  // Both the forces here are conservative, hence forces on 2 is minus total on 1
  lower_force += spring_force;
  upper_force = -lower_force;

  lower_torque += cross( lower_head - lower_cell->getPos(),  spring_force );
  upper_torque += cross( upper_tail - upper_cell->getPos(), -spring_force );
}

void getSpringEnergy(
  IBacterium* lower_cell, IBacterium* upper_cell,
  double &bend_energy, double &spring_energy
)
{
  // constexpr double radius { 0.5*constants::nondim_rodSpheroDiam };
  const Vec3 lower_cell_n { lower_cell->getOrientation() };
  const Vec3 upper_cell_n { upper_cell->getOrientation() };

  // The top of the lower cell's head to which the bottom of the spring will attach
  const Vec3 lower_head {
    lower_cell->getPos()
    + lower_cell_n*(
        // radius
        + 0.5*lower_cell->getLength()
      )
  };

  // The tail of the uppers cell's head to which the top of the spring will attach
  const Vec3 upper_tail {
    upper_cell->getPos()
    - upper_cell_n*(
        // radius
        + 0.5*upper_cell->getLength()
      )
  };

  const Vec3 low_to_high { upper_tail-lower_head };
  const double mod_l_to_h { low_to_high.norm() };

  const double inv_mod_l_to_h { 1.0 / mod_l_to_h };
  const Vec3 ux12 { low_to_high * inv_mod_l_to_h };


  /*-------------------*/
  /*   Stiff linking   */
  /*-------------------*/
  const Vec3 lower_rod_base {
    lower_head - lower_cell_n*constants::nondim_candidaRadius*2
  };
  const Vec3 upper_rod_base {
    upper_tail + upper_cell_n*constants::nondim_candidaRadius*2
  };

  const Vec3 t0 { lower_head - lower_rod_base };

  // Check there is no performance penalty for this
  const Vec3 t1 { low_to_high };

  const Vec3 t2 { upper_rod_base - upper_tail };

  // define inverse lengths
  const double ds1_inv { 2/( t0.norm() + t1.norm() ) };
  const double ds2_inv { 2/( t1.norm() + t2.norm() ) };

  // Alias for bending rigidity
  const double K { RodShapedBacterium::mBendRig };

  bend_energy=0.0;
  bend_energy+=0.5*K*dot2(calcBiNormal(t0,t1))*ds1_inv;
  bend_energy+=0.5*K*dot2(calcBiNormal(t1,t2))*ds2_inv;

  double en_pr_fac { 0.5*RodShapedBacterium::mKappa };
  spring_energy=en_pr_fac*dot2( mod_l_to_h -constants::nondim_candidaRadius*1.3 );//previously was ***-> constants::nondim_candidaRadius*2 );
}

inline
void computeChainingInteractions( IBacterium* cell )
{
  // Temporary variables
  Vec3 lower_force;
  Vec3 lower_torque;

  Vec3 upper_force;
  Vec3 upper_torque;

  if ( cell->getUpperLink() )
  {
    computeSingleChainInteraction(
      cell,
      cell->getUpperLink(),
      lower_force,
      upper_force,
      lower_torque,
      upper_torque
    );
    cell->addForce(lower_force);
    cell->addTorque(lower_torque);
  }

  if ( cell->getLowerLink() )
  {
    computeSingleChainInteraction(
      cell->getLowerLink(),
      cell,
      lower_force,
      upper_force,
      lower_torque,
      upper_torque
    );
    cell->addForce(upper_force);
    cell->addTorque(upper_torque);
  }
}
#endif
