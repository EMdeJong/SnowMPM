//----------------------------------------------------------------------------------------------------------------------
/// @file Particle.h
/// @brief Container for the particles
/// @author Esther de Jong
/// @version 2.0
/// @date 20-08-2015
/// Revision History :
/// Added weight and gradient weight values to particles
/// Get and set functions
//----------------------------------------------------------------------------------------------------------------------

#ifndef PARTICLE_H
#define PARTICLE_H

#include <ngl/Vec3.h>
#include <ngl/Mat3.h>

class Particle
{
public:
  Particle(ngl::Vec3 _minPos, ngl::Vec3 _maxPos, ngl::Vec3 _centrepos);

  //set all the particle values that are calculated in the grid class
  inline void SetParticleDensity(float _pDensity) {m_pDensity=_pDensity;}
  inline void SetParticleVolume(float _pVolume) {m_pVolume=_pVolume;}
  inline void SetDeformationGradient(ngl::Mat3 _pDeformationGradient) {m_PDeformationGradient=_pDeformationGradient;}
  inline void SetVelocityGradient(ngl::Mat3 _pVelocityGradient) {m_pVelocityGradient=_pVelocityGradient;}
  inline void SetParticleVelocity(ngl::Vec3 _pVelocity) {m_pVelocity=_pVelocity;}
  inline void SetElasticDefG(ngl::Mat3 _elasticDefG) {m_elasticDefG=_elasticDefG;}
  inline void SetdR(ngl::Mat3 _dR){m_dR=_dR;}
  inline void SetdElasticDefG(ngl::Mat3 _dElasticDefG){m_dElasticDefG=_dElasticDefG;}
  inline void SetAp(ngl::Mat3 _Ap){ m_Ap=_Ap;}
  inline void SetMu(float _mu) {m_mu=_mu;}
  inline void SetLambda(float _lambda) {m_lambda=_lambda;}
  inline void SetFeHat(ngl::Mat3 _FeHat) {m_FeHat=_FeHat;}
  inline void SetReHat(ngl::Mat3 _ReHat) {m_ReHat=_ReHat;}
  inline void SetSeHat(ngl::Mat3 _SeHat) {m_SeHat=_SeHat;}
  inline void SetParticlePosition(ngl::Vec3 _pPosition){m_PPosition=_pPosition;}

  //set the particle range for the grid nodes influence
  void SetParticleRange(float _xMin, float _xMax, float _yMin,float _yMax, float _zMin, float _zMax);
  // set the elastic and plastic deformation gradient
  void SetElasticAndPlasticDefG(ngl::Mat3 _elasticDefG, ngl::Mat3 _plasticDefG);
  //particle collision
  void ParticleCollision();
  //update particle position
  void UpdateParticlePosition(float _timeStep);

private:
  bool m_collision;
public:
  ngl::Mat3 m_elasticDefG;
  ngl::Mat3 m_plasticDefG;
  ngl::Mat3 m_PDeformationGradient;  
  ngl::Mat3 m_pVelocityGradient;
  ngl::Mat3 m_pTempMat;
  ngl::Mat3 m_dElasticDefG;
  ngl::Mat3 m_dR;
  ngl::Mat3 m_Ap;
  ngl::Mat3 m_FeHat;
  ngl::Mat3 m_ReHat;
  ngl::Mat3 m_SeHat;

  ngl::Vec3 m_PPosition;
  ngl::Vec3 m_pVelocity;
  ngl::Vec3 m_minPos, m_maxPos;
  ngl::Vec3 m_collisionNormal;
  ngl::Vec3 m_velocityTang;
  ngl::Vec3 m_collisionVelocity;

  //test values needed for sphere distribution
  //float m_test;
  //float m_test2;
  //float m_testv;
  float m_timeStep;
  float m_pDensity;
  float m_pMass;
  float m_pVolume;
  float m_hardeningCoefficient;
  float m_weightX, m_weightY, m_weightZ, m_weight;
  float m_dWeightX, m_dWeightY, m_dWeightZ;
  float m_vn, m_vtLength;
  float m_mu, m_lambda;
  float m_xMin, m_xMax;
  float m_yMin, m_yMax;
  float m_zMin, m_zMax;
};

#endif // PARTICLE_H
