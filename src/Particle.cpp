//----------------------------------------------------------------------------------------------------------------------
/// @file Particle.cpp
/// @brief
/// @author Esther de Jong
/// @version 1.2
/// @date 15-07-2015
/// Revision History :
/// This is an initial version
/// Added particle position update
/// Added draw function
//----------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <ngl/Random.h>

#include "Particle.h"


Particle::Particle(ngl::Vec3 _minPos, ngl::Vec3 _maxPos, ngl::Vec3 _centrepos)
{
  /*circle distribution can be used if no obj is imported
  ngl::Random *rand=ngl::Random::instance();
  m_minPos=_centrepos;
  m_test=rand->randomPositiveNumber();
  m_test=0.01814*pow(m_test,1.0/3.0); //0.0310175 0.0181391 0.0152
  m_test2=rand->randomNumber();
  m_testv=2.0*ngl::PI*rand->randomPositiveNumber();
  m_PPosition.m_z=m_test*m_test2;
  m_test2=1-m_test2*m_test2;
  m_test2=sqrt(m_test2);
  m_PPosition.m_x=m_test*m_test2*cos(m_testv);
  m_PPosition.m_y=m_test*m_test2*sin(m_testv);
  m_PPosition=m_PPosition+m_minPos;
  */

  m_PPosition=_minPos;

  //variation in hardening coefficient within the sphere depending on the dinstance from the centre
  if(m_PPosition.length()>0.015)
  {
    m_hardeningCoefficient=5;
  }
  else
  {
    m_hardeningCoefficient=10;
  }
  m_pVelocity=_maxPos;
  m_pMass=0.00000244;
  m_pVelocityGradient.identity();
  m_PDeformationGradient.identity();
  m_elasticDefG.identity();
  m_plasticDefG.identity();
  m_collisionNormal=ngl::Vec3(0,1,0);
}

void Particle::SetElasticAndPlasticDefG(ngl::Mat3 _elasticDefG, ngl::Mat3 _plasticDefG)
{
  m_elasticDefG=_elasticDefG;
  m_plasticDefG=_plasticDefG;
}

void Particle::ParticleCollision()
{
  m_collision=false;
  if(m_PPosition.m_y<= 0)
  {
    m_collisionNormal=ngl::Vec3(0,1,0);
    m_collision=true;
  }
  else if(m_PPosition.m_y>= 1)
  {
    m_collisionNormal=ngl::Vec3(0,-1,0);
    m_collision=true;
  }
  if(m_PPosition.m_x>= 0.1)
  {
    m_collisionNormal=ngl::Vec3(-1,0,0);
    m_collision=true;
  }
  else if(m_PPosition.m_x<= -0.05)
  {
    m_collisionNormal=ngl::Vec3(1,0,0);
    m_collision=true;
  }
  if(m_PPosition.m_z>= 0.1)
  {
    m_collisionNormal=ngl::Vec3(0,0,-1);
    m_collision=true;
  }
  else if(m_PPosition.m_z<= -0.1)
  {
    m_collisionNormal=ngl::Vec3(0,0,1);
    m_collision=true;
  }
  if(m_collision==true)
  {
    m_vn=m_pVelocity.dot(m_collisionNormal);
    if(m_vn<0)
    {
      m_velocityTang=m_pVelocity-m_collisionNormal*m_vn;
      m_vtLength=m_velocityTang.length();
      if(m_vtLength<=(-0.1*m_vn))
      {
        m_collisionVelocity=0;
      }
      else
      {
        m_collisionVelocity=m_velocityTang+(0.1*m_vn*m_velocityTang)/m_vtLength;
      }
      m_pVelocity=m_collisionVelocity;
    }
  }
}

void Particle::UpdateParticlePosition(float _timeStep)
{
  m_timeStep=_timeStep;
  m_PPosition+=m_timeStep*m_pVelocity;
}

void Particle::SetParticleRange(float _xMin, float _xMax, float _yMin,float _yMax, float _zMin, float _zMax)
{
  m_xMin=_xMin;
  m_xMax=_xMax;
  m_yMin=_yMin;
  m_yMax=_yMax;
  m_zMin=_zMin;
  m_zMax=_zMax;
}
