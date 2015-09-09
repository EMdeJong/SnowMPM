#include <iostream>
#include "ComputePolarDecomposition.h"


ComputePolarDecomposition::ComputePolarDecomposition()
{

}

ngl::Mat3 ComputePolarDecomposition::ComputePDRotation(ngl::Mat3 _elasticDG)
{
  ngl::Mat3 m_elasticF=_elasticDG;
  Eigen::Matrix3f m_matrixPD;
  //convert ngl matrix to eigen matrix
  m_matrixPD(0,0)=m_elasticF.m_00;
  m_matrixPD(0,1)=m_elasticF.m_01;
  m_matrixPD(0,2)=m_elasticF.m_02;
  m_matrixPD(1,0)=m_elasticF.m_10;
  m_matrixPD(1,1)=m_elasticF.m_11;
  m_matrixPD(1,2)=m_elasticF.m_12;
  m_matrixPD(2,0)=m_elasticF.m_20;
  m_matrixPD(2,1)=m_elasticF.m_21;
  m_matrixPD(2,2)=m_elasticF.m_22;
  m_svd.compute(m_matrixPD, Eigen::ComputeFullU |Eigen::ComputeFullV);
  Eigen::MatrixXf m_matrixUt=m_svd.matrixU();
  Eigen::MatrixXf m_matrixVt=m_svd.matrixV();
  Eigen::MatrixXf m_rotationMat=m_matrixUt*m_matrixVt.transpose();

  //convert eigen matrix to ngl matrix
  ngl::Mat3 m_rotationR;
  m_rotationR.m_00=m_rotationMat(0,0);
  m_rotationR.m_01=m_rotationMat(0,1);
  m_rotationR.m_02=m_rotationMat(0,2);
  m_rotationR.m_10=m_rotationMat(1,0);
  m_rotationR.m_11=m_rotationMat(1,1);
  m_rotationR.m_12=m_rotationMat(1,2);
  m_rotationR.m_20=m_rotationMat(2,0);
  m_rotationR.m_21=m_rotationMat(2,1);
  m_rotationR.m_22=m_rotationMat(2,2);



  return m_rotationR;
}

ngl::Mat3 ComputePolarDecomposition::ComputeRandS(ngl::Mat3 _elasticDG)
{
  ngl::Mat3 m_elasticF=_elasticDG;
  Eigen::Matrix3f m_matrixPD;
  //convert ngl matrix to eigen matrix
  m_matrixPD(0,0)=m_elasticF.m_00;
  m_matrixPD(0,1)=m_elasticF.m_01;
  m_matrixPD(0,2)=m_elasticF.m_02;
  m_matrixPD(1,0)=m_elasticF.m_10;
  m_matrixPD(1,1)=m_elasticF.m_11;
  m_matrixPD(1,2)=m_elasticF.m_12;
  m_matrixPD(2,0)=m_elasticF.m_20;
  m_matrixPD(2,1)=m_elasticF.m_21;
  m_matrixPD(2,2)=m_elasticF.m_22;
  m_svd.compute(m_matrixPD, Eigen::ComputeFullU |Eigen::ComputeFullV);
  Eigen::MatrixXf m_matrixUt=m_svd.matrixU();
  Eigen::MatrixXf m_matrixVt=m_svd.matrixV();

  Eigen::MatrixXf m_rotationMat=m_matrixUt*(m_matrixVt.adjoint());

  Eigen::Vector3f m_singularValuest=m_svd.singularValues();
  Eigen::MatrixXf m_singularMatt=Eigen::Matrix3f::Zero();
  for(int i=0; i<3; ++i)
  {
    m_singularMatt(i,i)=m_singularValuest(i);
  }
  Eigen::MatrixXf m_matrixSt=m_matrixVt*m_singularMatt*(m_matrixVt.adjoint());

  //convert Eigen matrix to NGL matrix
  ngl::Mat3 m_rotationR;
  m_rotationR.m_00=m_rotationMat(0,0);
  m_rotationR.m_01=m_rotationMat(0,1);
  m_rotationR.m_02=m_rotationMat(0,2);
  m_rotationR.m_10=m_rotationMat(1,0);
  m_rotationR.m_11=m_rotationMat(1,1);
  m_rotationR.m_12=m_rotationMat(1,2);
  m_rotationR.m_20=m_rotationMat(2,0);
  m_rotationR.m_21=m_rotationMat(2,1);
  m_rotationR.m_22=m_rotationMat(2,2);

  m_sPD.m_00=m_matrixSt(0,0);
  m_sPD.m_01=m_matrixSt(0,1);
  m_sPD.m_02=m_matrixSt(0,2);
  m_sPD.m_10=m_matrixSt(1,0);
  m_sPD.m_11=m_matrixSt(1,1);
  m_sPD.m_12=m_matrixSt(1,2);
  m_sPD.m_20=m_matrixSt(2,0);
  m_sPD.m_21=m_matrixSt(2,1);
  m_sPD.m_22=m_matrixSt(2,2);

  return m_rotationR;
}

void ComputePolarDecomposition::ComputeElasticAndPlasticDefGrad(ngl::Mat3 _elasticDG, ngl::Mat3 _defGrad)
{
  m_singularMat=Eigen::Matrix3f::Zero();
  for(int i=0; i<3; ++i)
  {
    m_singularMat(i,i)=m_singularValues(i);
  }

  //recalculate FE
  m_matrixA=m_matrixU*m_singularMat*(m_matrixV.transpose().conjugate());

  m_elasticDG.m_00=m_matrixA(0,0);
  m_elasticDG.m_01=m_matrixA(0,1);
  m_elasticDG.m_02=m_matrixA(0,2);
  m_elasticDG.m_10=m_matrixA(1,0);
  m_elasticDG.m_11=m_matrixA(1,1);
  m_elasticDG.m_12=m_matrixA(1,2);
  m_elasticDG.m_20=m_matrixA(2,0);
  m_elasticDG.m_21=m_matrixA(2,1);
  m_elasticDG.m_22=m_matrixA(2,2);

  //calculate FP still needs to be multiplied with F
  m_matrixA=m_matrixV*m_singularMat.inverse()*m_matrixU.transpose();
  m_plasticDG.m_00=m_matrixA(0,0);
  m_plasticDG.m_01=m_matrixA(0,1);
  m_plasticDG.m_02=m_matrixA(0,2);
  m_plasticDG.m_10=m_matrixA(1,0);
  m_plasticDG.m_11=m_matrixA(1,1);
  m_plasticDG.m_12=m_matrixA(1,2);
  m_plasticDG.m_20=m_matrixA(2,0);
  m_plasticDG.m_21=m_matrixA(2,1);
  m_plasticDG.m_22=m_matrixA(2,2);
}

bool ComputePolarDecomposition::ComputeSVD(ngl::Mat3 _matrixA)
{
  //the upper and lower clamp value
  m_singularLower=1-0.025;
  m_singularUpper=1+0.0075;
  //from NGL to Eigen matrix
  m_matrixA(0,0)=_matrixA.m_00;
  m_matrixA(0,1)=_matrixA.m_01;
  m_matrixA(0,2)=_matrixA.m_02;
  m_matrixA(1,0)=_matrixA.m_10;
  m_matrixA(1,1)=_matrixA.m_11;
  m_matrixA(1,2)=_matrixA.m_12;
  m_matrixA(2,0)=_matrixA.m_20;
  m_matrixA(2,1)=_matrixA.m_21;
  m_matrixA(2,2)=_matrixA.m_22;

  m_svd.compute(m_matrixA, Eigen::ComputeFullU |Eigen::ComputeFullV);
  m_matrixU=m_svd.matrixU();
  m_matrixV=m_svd.matrixV();
  m_singularValues=m_svd.singularValues();
  m_singularMat=Eigen::Matrix3f::Zero();
  for(int i=0; i<3; ++i)
  {
    m_singularMat(i,i)=m_singularValues(i);
  }
  m_outsideRange=0;
  //see if any of the singular values needs to be clamped
  if(m_singularValues(0)>m_singularUpper)
  {
    m_outsideRange=1;
    m_singularValues(0)=m_singularUpper;
    if(m_singularValues(1)>m_singularUpper)
    {
      m_singularValues(1)=m_singularUpper;
      if(m_singularValues(2)>m_singularUpper)
      {
        m_singularValues(2)=m_singularUpper;
      }
    }
  }
  if(m_singularValues(2)<m_singularLower)
  {
    m_outsideRange=1;
    m_singularValues(2)=m_singularLower;
    if(m_singularValues(1)<m_singularLower)
    {
      m_singularValues(1)=m_singularLower;
      if(m_singularValues(0)<m_singularLower)
      {
        m_singularValues(0)=m_singularLower;
      }
    }
  }

  return m_outsideRange;
}
