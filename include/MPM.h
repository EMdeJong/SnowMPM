//----------------------------------------------------------------------------------------------------------------------
/// @file MPM.h
/// @brief Calls all the functions for the material point method
/// @author Esther de Jong
/// @version 1.0
/// @date 23-06-2015
/// Revision History :
/// This is an initial version
//----------------------------------------------------------------------------------------------------------------------

#ifndef MPM_H
#define MPM_H

#include <vector>
#include "ngl/Vec3.h"
#include "ngl/Camera.h"
#include "InterpolationWeight.h"
#include "Grid.h"

class MPM
{
public:
  MPM(int _numParticles, float _timer);

  //initialize particles inside mesh (start with simple sphere)
  void CreateParticles(int _numParticles);
  //perform grid operations
  void UpdateGridOperations();
private:
  bool m_firstStep;
  int m_numParticles;
  float m_timer;  
  Grid* m_grid;
  InterpolationWeight *m_iPW;
};

#endif // MPM_H
