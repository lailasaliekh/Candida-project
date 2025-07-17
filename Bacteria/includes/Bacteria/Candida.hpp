#ifndef CANDIDA_HPP
#define CANDIDA_HPP

#include "RodShapedBacteria.hpp"

class Candida : public RodShapedBacterium
{
public:
    static double mRadius;      //!< Candida-specific radius
    static double mAvgDivLen;      //!< Candida-specific division length
    static double mAvgGrwthRate;   //!< Candida-specific growth rate
    double mGrwthRtePreFac;
    double mLinkingProb; 
    /*---------------------- Set cell dynamic properties -----------------------*/
    // // Cell posn and dynamics
    // Vec3 mPos;                          //!< centre of mass vector
    // Vec3 mLogPos;                       //!< centre of mass vec at time of binning
    // Vec3 mVel{0,0,0};                   //!< centre of mass velocity
    // Vec3 mForce{0,0,0};                 //!< net force experienced by this particle
    // Vec3 mAngles{0,0.5*constants::pi,0};//!< phi (xy), theta (z), psi (body axis)
    // Vec3 mAngVel{0,0,0};                //!< angular velocity in the body frame
    // Vec3 mTorque{0,0,0};                //!< net torque experienced by this particle (fixed frame)
    // double mLength;                     //!< Length of the cylindircal part
 
    // static uint counter;        
    // uint mId;       
    Candida(
        double rcm_x = 0,
        double rcm_y = 0,
        double rcm_z = 0,
        double theta = 0,
        double alpha = constants::pi * 0.5,
        double grwthPreFac = constants::nondim_candidaGrwthRtePreFac,
        double init_length = constants::nondim_candidaInitLength,
        double linking_prob = 1.0
    );

    Candida(
        const Vec3 &rcm,
        double theta = 0,
        double alpha = constants::pi * 0.5,
        double grwthPreFac = constants::nondim_candidaGrwthRtePreFac,
        double init_length = constants::nondim_candidaInitLength,
        double linking_prob = 1.0
    );
    virtual std::string getMyType() const override
    {
        return "Candida";
    }
    // void splitChainIfNecessary(int splitThreshold, std::vector<IBacterium*>& cell_list);
//     virtual double getRadius() const override
//     {
//         return mRadius;
//     }
//     virtual double getEffectiveR() const override
//     {
//         return mRadius + 0.5 * mLength; // Use mCandidaRadius instead of mRadius
//     }
//     void grow(double dt) override
//     {
//         mLength += dt * mCandidaGrwthRtePreFac; // Use Candida's specific growth rate
//     }
//     virtual bool signalDivide() override
//     {
//         return mLength > mCandidaDivLen;
//     }
    virtual void divide(std::vector<IBacterium*>& cell_list) override;
//     virtual double getLength() const override
//     {
//         return mLength;
//     }
//   virtual Vec3& getForce() override
//   {
//     return mForce;
//   }
//   virtual void addForce(Vec3 force) override
//   {
//     mForce+=force;
//   }
//   virtual Vec3& getTorque() override
//   {
//     return mTorque;
//   }
//   virtual void addTorque(Vec3 torque) override
//   {
//     mTorque+=torque;
//   }
//   virtual Vec3& getVel() override
//   {
//     return mVel;
//   }
//   virtual void setVel() override
//   {
//     mVel=mForce/mLength;
//   }

//   // Return the lab frame angvel
//   virtual Vec3& getAngVel() override
//   {
//     return mAngVel;
//   }

//   // Set the lab frame angvel from the torque on this rod
//   virtual void setAngVel() override
//   {
// #ifndef ANISOTROPIC
//     // const double effective_length_3{ mRadius*mRadius*( mRadius + 0.75*mLength ) };
//     // mAngVel=(4.0/3.0)*mTorque/effective_length_3;
//      mAngVel = mTorque * 12 / (mLength*mLength*mLength);
// #else
//     std::cout << "exit from set angvel rod shaped" << '\n';
//     exit(1);
// #endif
//   }
//   virtual Vec3 getPos() const override
//   {
//     return mPos;
//   }

//   virtual void setPos(double x, double y, double z=0.0) override
//   {
//     mPos.x=x; mPos.y=y; mPos.z=z;
//   }

//   virtual void setAngles(double theta,double alpha=0.5*constants::pi) override
//   {
//     mAngles.x=theta;
//     mAngles.y=alpha;
//   }

//   virtual Vec3 getLoggedPos() const override
//   {
//     return mLogPos;
//   }

//   virtual void setLoggedPos() override
//   {
//     mLogPos=mPos;
//   }
//   virtual void move(double dt) override
//   {
    
//     // Get the vector ( dphi / dt, dtheta /dt , dpsi / dt )
//     Vec3 bodyFrameAngVel { projectAngVelRod(mAngVel,mAngles) };
//     #ifndef MOVE_3D
//     bodyFrameAngVel.y=0;
//     mVel.z=0;
//     #endif
//     // Should phi be mod 2 pi?
//     mAngles+=bodyFrameAngVel*dt;
//     mPos+=mVel*dt;

//   }
//   void getMyEndVecs(Vec3& p, Vec3& q) const override
//     {
//     getEndVecs(*this,p,q);
//     }
//     virtual Vec3 getOrientation() const override
//     {
//           Vec3 n_hat { getOrientationFromAngles(mAngles) };
// #ifndef MOVE_3D
//   n_hat.z=0;
// #endif
//   return n_hat;
//     }
    virtual void printToFile(std::ostream &out) override //to correctly print data of Candida
    {
    #ifdef CHAINING
        out << getMyType()       << "\t";          // Candida
        out << mId               << "\t";          // Unique ID
        out << mLength           << "\t";          // Current length
        out << mRadius    << "\t";          // Specific radius for Candida
        out << mPos              << "\t";          // Position vector
        out << getOrientation()  << "\t";          // Orientation vector
        if (getLowerLink()) out << getLowerLink()->getID() << "\t"; // Link to lower cell
        else out << "None" << "\t";
        if (getUpperLink()) out << getUpperLink()->getID() << "\n"; // Link to upper cell
        else out << "None" << "\n";
    #else
        out << getMyType()       << "\t";
        out << mId               << "\t";
        out << mLength           << "\t";
        out << mRadius    << "\t";
        out << mPos              << "\t";
        out << getOrientation()  << "\n";
    #endif
    }

};
std::ostream &printCellDataToFile(std::ostream &out, const Candida &cell);

#endif // CANDIDA_HPP
