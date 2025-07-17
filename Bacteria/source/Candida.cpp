#include "RodShapedBacteria.hpp"
#include "Candida.hpp"
#include "constants.hpp"
// Initialize static members
double Candida::mRadius {constants::nondim_candidaRadius};  // larger radius
double Candida::mAvgDivLen  {constants::nondim_candidaAvgDivLen};  //  larger division length
double Candida::mAvgGrwthRate {constants::nondim_candidaGrwthRtePreFac};  //  slower growth rate
// uint Candida:: counter { 0 }; 
 
Candida::Candida(
  double _x,
  double _y,
  double _z,
  double theta,
  double alpha,
  double grwthPreFac,
  double init_length,
  double linking_prob
) : RodShapedBacterium(_x, _y, _z, theta, alpha, grwthPreFac, init_length, linking_prob)
//  mPos{_x,_y,_z}, mAngles {theta,alpha,0}, mLength {init_length},
//  mCandidaGrwthRtePreFac{grwthPreFac},mLinkingProb{linking_prob}
{
  mId = counter++;           // increment unique counter
  
  // Prevent avalanche of divisions
  assert( mLength <= mAvgDivLen );

}

Candida::Candida(
    const Vec3 &rcm,
    double theta,
    double alpha,
    double grwthPreFac,
    double init_length,
    double linking_prob):
     RodShapedBacterium{rcm.x,rcm.y,rcm.z, theta, alpha, grwthPreFac, init_length, linking_prob}
{}
std::ostream &printCellDataToFile(std::ostream &out, const Candida &cell)
{
    out << "Candida\t";
    out << cell.mId << "\t";
    out << cell.mLength << "\t";
    out << cell.mRadius << "\t"; // Specific radius for Candida
    out << cell.mPos << "\t";
    out << cell.getOrientation() << "\n";

    return out;
}
