//----------------------------------------------------------------------------------------------------------------------
/// @file InterpolationWeight.h
/// @brief calculate the interpolated weight and its gradient
/// @author Esther de Jong
/// @version 1.1
/// @date 24-06-2015
/// Revision History :
/// Added extra function instead of constructor
/// @class InterpolationWeight
/// @brief calculate the interpolated weight and its gradient
//----------------------------------------------------------------------------------------------------------------------

#include <math.h>
#include <iostream>
#include "InterpolationWeight.h"

InterpolationWeight::InterpolationWeight()
{

}

float InterpolationWeight::CalInterpolationWeight(float _x)
{
  float m_x=_x;
  float m_xAbs=fabs(m_x);
  float m_intWeight;
  if(m_xAbs>=0&&m_xAbs<1)
  {
    m_intWeight=0.5*pow(m_xAbs, 3.0)-pow(m_x,2.0)+2.0/3.0;
  }
  else if(m_xAbs>=1&&m_xAbs<2)
  {
    m_intWeight=(-1.0/6.0)*pow(m_xAbs, 3.0)+pow(m_x, 2.0)-2*m_xAbs+4.0/3.0;
  }
  else
  {
    m_intWeight=0;
  }
  return m_intWeight;
}

float InterpolationWeight::IWGradient(float _x)
{
  float m_xAbs=fabs(_x);
  float m_x=_x;
  float m_intWeight;
  if(m_xAbs>=0&&m_xAbs<1)
  {
    m_intWeight=(3.0/2.0)*pow(m_xAbs, 2.0)-2*m_xAbs;
  }
  else if(m_xAbs>=1&&m_xAbs<2)
  {
    m_intWeight=(-0.5)*pow(m_xAbs, 2.0)+2.0*m_xAbs-2;
  }
  else
  {
    m_intWeight=0;
  }
  if(m_x<0)
  {
    m_intWeight*=-1.0;
  }
  return m_intWeight;
}
