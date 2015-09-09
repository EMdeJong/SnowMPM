//----------------------------------------------------------------------------------------------------------------------
/// @file InterpolationWeight.h
/// @brief calculate the interpolated weight and its gradient
/// @author Esther de Jong
/// @version 1.2
/// @date 20-08-2015
/// Revision History :
/// Added extra function instead of constructor
/// Fixed the weight gradient function
//----------------------------------------------------------------------------------------------------------------------

#ifndef INTERPOLATIONWEIGHT_H
#define INTERPOLATIONWEIGHT_H

class InterpolationWeight
{
public:
  InterpolationWeight();
  //calculates the interpolation weight for the input x
  float CalInterpolationWeight(float _x);
  //calculates the derivative of the interpolation weight for the input x
  float IWGradient(float _x);
};

#endif // INTERPOLATIONWEIGHT_H
