//----------------------------------------------------------------------------------------------------------------------
/// @file MPM.cpp
/// @brief Calls all the MPM steps in the right order
/// @author Esther de Jong
/// @version 2.0
/// @date 23-08-2015
/// Revision History :
/// Added all the functions for the MPM
//----------------------------------------------------------------------------------------------------------------------

#include "MPM.h"

MPM::MPM(int _numParticles, float _timer)
{
  m_timer=_timer;
  m_numParticles=_numParticles;
  m_grid=new Grid(m_numParticles, ngl::Vec3 (-0.12,0.0,-0.12), ngl::Vec3 (0.12,0.1,0.12), 0.003, m_timer);
  m_firstStep=0;
}

void MPM::UpdateGridOperations()
{
  m_grid->ComputeGridMass();
  m_grid->ComputeGridVelocity();
  if(m_firstStep==0)
  {
    m_grid->ComputeParticleVolumes();
    m_firstStep=1;
  }
  m_grid->ComputeGridForces();
  m_grid->UpdateGridVelocities();
  m_grid->GridCollision();
  m_grid->ImplicitSolver();
  m_grid->CalculateDeformationGradient();
  m_grid->UpdateParticleVelocities();
  m_grid->ExplicitSolver();
}
