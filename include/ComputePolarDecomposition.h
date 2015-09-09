//----------------------------------------------------------------------------------------------------------------------
/// @file ComputePolarDecomposition.h
/// @brief Computes the polar and singular value decomposition
/// @author Esther de Jong
/// @version 1.3
/// @date 14-07-2015
/// Revision History :
/// Added the calculations to get the Rotation of the polar decomposition matrix
/// Used the Eigen library to compute the correct rotation from the polar decomposition
/// Added the deformation gradient update step which inclused updating and clamping the elastic part
//----------------------------------------------------------------------------------------------------------------------
#ifndef COMPUTEPOLARDECOMPOSITION_H
#define COMPUTEPOLARDECOMPOSITION_H

#include "ngl/Mat3.h"
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Dense>


class ComputePolarDecomposition
{
public:
  ComputePolarDecomposition();
  //compute the rotation matrix for the elastic deformation gradient
  ngl::Mat3 ComputePDRotation(ngl::Mat3 _elasticDG);
  //compute the rotation and transormation matrix for the deformation gradient
  ngl::Mat3 ComputeRandS(ngl::Mat3 _elasticDG);
  //compute the elastic and plastic deformation gradient with the clamped elastic matrix
  void ComputeElasticAndPlasticDefGrad(ngl::Mat3 _elasticDG, ngl::Mat3 _defGrad);
  //see if the elastic deformation gradient has to be clamped
  bool ComputeSVD(ngl::Mat3 _matrixA);

  inline ngl::Mat3 GetPDS(){return m_sPD;}
  inline ngl::Mat3 GetElasticDG(){return m_elasticDG;}
  inline ngl::Mat3 GetPlasticDG(){return m_plasticDG;}

private:
  bool m_outsideRange;
  float m_singularUpper;
  float m_singularLower;
  Eigen::Matrix3f m_matrixA, m_temp3Mat, m_singularMat;
  Eigen::MatrixXf m_rotationMat, m_matrixU, m_matrixV, m_tempMat, m_matrixS;
  Eigen::Vector3f m_singularValues;
  Eigen::JacobiSVD<Eigen::Matrix3f> m_svd;
  ngl::Mat3 m_elasticDG;
  ngl::Mat3 m_plasticDG;
  ngl::Mat3 m_diagonalPD;
  ngl::Mat3 m_rotationPD;
  ngl::Mat3 m_sPD;
  ngl::Mat3 m_ATranspose;
  ngl::Mat3 m_V;
  ngl::Mat3 m_U;
};

#endif // COMPUTEPOLARDECOMPOSITION_H
