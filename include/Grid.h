//----------------------------------------------------------------------------------------------------------------------
/// @file Grid.h
/// @brief contains all the Material Point Method functions that need the grid
/// @author Esther de Jong
/// @version 2.0
/// @date 23-08-2015
/// Revision History :
/// Added and changed some functions
//----------------------------------------------------------------------------------------------------------------------

#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cstring>
#include <math.h>

#include <omp.h>

#include <ngl/Vec3.h>
#include <ngl/Random.h>

#include "Particle.h"
#include "InterpolationWeight.h"
#include "ComputePolarDecomposition.h"
#include "ExportParticleData.h"

typedef struct GridNode{
  float m_nodeMass;
  ngl::Vec3 m_nodeVelocity;
  ngl::Vec3 m_newNodeVelocity;
  ngl::Vec3 m_nodeForce;
  ngl::Vec3 m_nodePosition;

  //implicit variables
  ngl::Vec3 m_r;
  ngl::Vec3 m_s;
  ngl::Vec3 m_p;
  ngl::Vec3 m_v;
  ngl::Vec3 m_df;
  ngl::Vec3 m_dFE;
  ngl::Vec3 m_gamma;
  ngl::Vec3 m_Ar;
  ngl::Vec3 m_Ap;
  float m_beta;
  float m_alpha;
  float m_residual;
  bool m_implicitSolved;
} GridNode;

class Grid
{
public:
  Grid(int _numParticles, ngl::Vec3 _pMin, ngl::Vec3 _pMax, float _gridSpacing, float _timer);

  //Compute velocity and mass for grid points
  void ComputeGridMass();
  //compute the grid node velocities from the particle velocities
  void ComputeGridVelocity();
  //first timestep compute particle volumes an densities
  void ComputeParticleVolumes();
  //compute grid velocities call the ParticleForce
  void ComputeGridForces();
  void ComputeFeReSeHat();
  //update the grid velocities after the grid forces are calculated
  void UpdateGridVelocities();
  //grid collisions
  void GridCollision();

  //solve the linear system explicitly
  void ExplicitSolver();
  //solve the linear system for semi-implicit integration
  void ImplicitSolver();
  //compute Ar for the implicit integration
  void ComputeAr();
  //calculate the plastic and elastic deformation gradients
  void CalculateDeformationGradient();
  //update the particle velocties using the weighted grid nodes
  void UpdateParticleVelocities();

private:
  bool m_outsideRange;
  bool m_implicitDone;

  int m_nodesAmount;
  int m_newNodesAmount;
  int m_weightDistance;
  int m_frameNumber;
  int m_frameUpdate;
  int m_partamount;

  float m_timeStep;
  float m_gridSpacing;
  float m_gridSpacingCubed;
  float m_poissons, m_youngsmodulus;

  //parameters are localized for openmp
  //float m_intermediateVelocity;  
  //float m_xMin, m_xMax;
  //float m_yMin, m_yMax;
  //float m_zMin, m_zMax;
  //float m_weightX, m_weightY, m_weightZ, m_weight;
  //float m_dWeightX, m_dWeightY, m_dWeightZ;
  //float m_partGridDist;  
  //float m_tempValue;
  //float m_nodeForce;
  //float m_jacobian, m_jElastic, m_jPlastic;
  //float m_pVolume;
  //float m_pDensity;
  //float m_lameFirstParameter, m_lameSecondParameter;
  //float m_mu, m_lambda;
  //float m_vn, m_vtLength;
  //float m_elasticDefG, m_plasticDefG;

  ngl::Vec3 m_centrePos;
  ngl::Vec3 m_gMin;
  ngl::Vec3 m_gMax;
  ngl::Vec3 m_newGMin, m_newGMax;
  ngl::Vec3 m_gridSize;
  ngl::Vec3 m_newGridSize;
  ngl::Vec3 m_gravity;
  ngl::Vec3 m_pVelocity;

  //parameters are localized for openmp
  //ngl::Vec3 m_gridPosition;
  //ngl::Vec3 m_weightGradient;
  //ngl::Vec3 m_tempVec;
  //ngl::Vec3 m_picVelocity;
  //ngl::Vec3 m_flipVelocity;
  //ngl::Vec3 m_collisionNormal;
  //ngl::Vec3 m_collisionVelocity;
  //ngl::Vec3 m_velocityTang;
  //ngl::Vec3 m_df;

  std::vector <Particle> particles;
  std::vector <GridNode> m_nodes;
  ExportParticleData *m_expPart;
  ComputePolarDecomposition *m_polarDecomposition;
  InterpolationWeight *m_iPWeight;

  //parameters are localized for openmp
  //ngl::Mat3 m_pVelocityGradient;
  //ngl::Mat3 m_gDeformationGradient;
  //ngl::Mat3 m_gElasticDefG;
  //ngl::Mat3 m_dElasticDefG;
  //ngl::Mat3 m_gPlasticDefG;
  //ngl::Mat3 m_rotationPD;
  //ngl::Mat3 m_sPD;
  //ngl::Mat3 m_sigma;

  ngl::Mat3 m_dR;
  ngl::Mat3 m_LH;  
  ngl::Mat3 m_identity;
  ngl::Mat3 m_dJFinvT;
  ngl::Mat3 m_JFinvT;
  ngl::Mat3 m_Ap;

};

#endif // GRID_H
