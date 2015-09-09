//----------------------------------------------------------------------------------------------------------------------
/// @file Grid.cpp
/// @brief
/// @author Esther de Jong
/// @version 3.0
/// @date 23-08-2015
/// Revision History :
/// Made the Mass, Velocity, Volume and Draw functions
/// Grid Visualization
/// Working on node force function
/// Finished node force function and changed functions so weight has to be calculated once per timestep
/// included timestep input and grid velocity update
/// Added explicit time integration
/// Added the calculation of the new deformation gradient and the elastic and plastic parts
/// Added particle velocity update step
/// Implementing grid collision
/// Added particle collision
/// Changeble grid spacing
/// Added implicit update
//----------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <ngl/Obj.h>
#include "Grid.h"

Grid::Grid(int _numParticles, ngl::Vec3 _pMin, ngl::Vec3 _pMax, float _gridSpacing, float _timer)
{
  m_timeStep=_timer;
  m_timeStep/=1000;

  m_gridSpacing=_gridSpacing;
  m_gMin=_pMin;
  m_gMax=_pMax;
  m_centrePos=ngl::Vec3(0.0, 0.1, 0.0);

  //scene setup
  ngl::Obj mesh;
  mesh.load("Models/SmallSnowball.obj", false);

  std::vector<ngl::Vec3> pos(mesh.getNumVerts());
  for(size_t i = 0; i < pos.size(); ++i)
  {
    pos[i] = mesh.getVertexAtIndex(i);
    particles.push_back(Particle(pos[i], ngl::Vec3(-10.0,0.0,0.0), m_centrePos));
  }
  m_partamount=pos.size();
  mesh.load("Models/BigSnowball.obj", false);
  std::vector<ngl::Vec3> pos2(mesh.getNumVerts());
  for(size_t i = 0; i < pos2.size(); ++i)
  {
    pos2[i] = mesh.getVertexAtIndex(i);
    particles.push_back(Particle(pos2[i], ngl::Vec3(0.0,0.0,0.0), m_centrePos));
  }
  m_partamount+=pos2.size();

  m_gMin=m_gMin-ngl::Vec3 (3*m_gridSpacing,3*m_gridSpacing,3*m_gridSpacing);
  m_gMax=m_gMax+ngl::Vec3 (3*m_gridSpacing,3*m_gridSpacing,3*m_gridSpacing);
  m_gridSize=(m_gMax-m_gMin)/m_gridSpacing+ngl::Vec3 (1,1,1);
  m_newNodesAmount=m_nodesAmount=(m_gridSize.m_x)*(m_gridSize.m_y)*(m_gridSize.m_z);

  for(int i=0; i<m_nodesAmount; ++i)
  {
    m_nodes.push_back(GridNode());
  }

  m_iPWeight= new InterpolationWeight();
  m_polarDecomposition=new ComputePolarDecomposition();
  m_poissons=0.2;
  m_youngsmodulus=1.4e5;
  m_identity.identity();
  m_partamount=particles.size();
  m_gravity=ngl::Vec3(0.0,-9.8,0.0);
  m_weightDistance=2;
  m_expPart= new ExportParticleData();
  m_frameUpdate=100;
  m_frameNumber=0;
}

void Grid::ComputeGridMass()
{  
  if(m_newNodesAmount>m_nodesAmount)
  {

    for(int i=0; i<(m_newNodesAmount-m_nodesAmount); ++i)
    {
      m_nodes.push_back(GridNode());
    }
    m_gridSize=m_newGridSize;
    m_nodesAmount=(m_gridSize.m_x)*(m_gridSize.m_y)*(m_gridSize.m_z);
  }
  else if(m_newNodesAmount<m_nodesAmount)
  {

    for(int i=0; i<(m_nodesAmount-m_newNodesAmount); ++i)
    {
      m_nodes.pop_back();
    }
    m_gridSize=m_newGridSize;
    m_nodesAmount=m_newNodesAmount;
  }

    for(int k=0; k<m_gridSize.m_z; ++k)
    {
      for(int j=0; j<m_gridSize.m_y; ++j)
      {
        for(int i=0; i<m_gridSize.m_x; ++i)
        {
          int m_currentnode=(k*(m_gridSize.m_x)*m_gridSize.m_y)+(j*m_gridSize.m_x)+i;
          m_nodes[m_currentnode].m_nodePosition.m_x=m_gMin.m_x+i*m_gridSpacing;
          m_nodes[m_currentnode].m_nodePosition.m_y=m_gMin.m_y+j*m_gridSpacing;
          m_nodes[m_currentnode].m_nodePosition.m_z=m_gMin.m_z+k*m_gridSpacing;
        }
      }
    }

  //reset all values here so there is only one loop per step

  #pragma omp parallel for
    for(int i=0; i<m_nodesAmount; ++i)
    {
      m_nodes[i].m_nodeMass=0;
      m_nodes[i].m_nodeVelocity=0;
      m_nodes[i].m_newNodeVelocity=0;
      m_nodes[i].m_nodeForce=0;
      m_nodes[i].m_dFE=0;
      m_nodes[i].m_df=0.0;
    }

    #pragma omp parallel for
    for(int i=0; i<m_partamount; ++i)
    {
      ngl::Vec3 m_gridPosition=(particles[i].m_PPosition-m_gMin)/m_gridSpacing;
      float m_xMin=ceil(m_gridPosition.m_x-m_weightDistance);
      float m_xMax=floor(m_gridPosition.m_x+m_weightDistance);
      float m_yMin=ceil(m_gridPosition.m_y-m_weightDistance);
      float m_yMax=floor(m_gridPosition.m_y+m_weightDistance);
      float m_zMin=ceil(m_gridPosition.m_z-m_weightDistance);
      float m_zMax=floor(m_gridPosition.m_z+m_weightDistance);
      particles[i].SetParticleRange(m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax);

      for(int x=particles[i].m_xMin; x<=particles[i].m_xMax; ++x)
      {
        float m_partGridDist=(m_gridPosition.m_x-x);
        float m_weightX=m_iPWeight->CalInterpolationWeight(m_partGridDist);
        for(int y=particles[i].m_yMin; y<=particles[i].m_yMax; ++y)
        {
          m_partGridDist=(m_gridPosition.m_y-y);
          float m_weightY=m_iPWeight->CalInterpolationWeight(m_partGridDist);
          for(int z=particles[i].m_zMin; z<=particles[i].m_zMax; ++z)
          {
            m_partGridDist=(m_gridPosition.m_z-z);
            float m_weightZ=m_iPWeight->CalInterpolationWeight(m_partGridDist);
            float m_weight=m_weightX*m_weightY*m_weightZ;
            int m_currentnode=(z*m_gridSize.m_x*m_gridSize.m_y)+(y*m_gridSize.m_x)+x;
            #pragma omp critical(nodeupdate0)
            {
              m_nodes[m_currentnode].m_nodeMass += particles[i].m_pMass*m_weight;
            }
          }
        }
      }
    }
}

void Grid::ComputeGridVelocity()
{
  m_partamount=particles.size();
  #pragma omp parallel for
    for(int i=0; i<m_partamount; ++i)
    {
      ngl::Vec3 m_gridPosition=(particles[i].m_PPosition-m_gMin)/m_gridSpacing;
      for(int x=particles[i].m_xMin; x<particles[i].m_xMax; ++x)
      {
        float m_partGridDist=(m_gridPosition.m_x-x);
        float m_weightX=m_iPWeight->CalInterpolationWeight(m_partGridDist);
        for(int y=particles[i].m_yMin; y<particles[i].m_yMax; ++y)
        {
          m_partGridDist=(m_gridPosition.m_y-y);
          float m_weightY=m_iPWeight->CalInterpolationWeight(m_partGridDist);
          for(int z=particles[i].m_zMin; z<particles[i].m_zMax; ++z)
          {
            m_partGridDist=(m_gridPosition.m_z-z);
            float m_weightZ=m_iPWeight->CalInterpolationWeight(m_partGridDist);
            float m_weight=m_weightX*m_weightY*m_weightZ;
            int m_currentnode=(z*m_gridSize.m_x*m_gridSize.m_y)+(y*m_gridSize.m_x)+x;
            if(m_nodes[m_currentnode].m_nodeMass!=0)
            {
              float m_tempValue=(particles[i].m_pMass*m_weight)/m_nodes[m_currentnode].m_nodeMass;
              #pragma omp critical(nodeupdate1)
              {
                m_nodes[m_currentnode].m_nodeVelocity += particles[i].m_pVelocity.operator *(m_tempValue);
              }
            }
          }
        }
      }
    }
}

void Grid::ComputeParticleVolumes()
{
  m_gridSpacingCubed=pow(m_gridSpacing, 3.0);
  #pragma omp parallel for
    for(int i=0; i<m_partamount; ++i)
    {
      float m_pDensity=0;
      float m_pVolume=0;
      ngl::Vec3 m_gridPosition=(particles[i].m_PPosition-m_gMin)/m_gridSpacing;
      for(int x=particles[i].m_xMin; x<particles[i].m_xMax; ++x)
      {
        float m_partGridDist=(m_gridPosition.m_x-x);
        float m_weightX=m_iPWeight->CalInterpolationWeight(m_partGridDist);
        for(int y=particles[i].m_yMin; y<particles[i].m_yMax; ++y)
        {
          m_partGridDist=(m_gridPosition.m_y-y);
          float m_weightY=m_iPWeight->CalInterpolationWeight(m_partGridDist);
          for(int z=particles[i].m_zMin; z<particles[i].m_zMax; ++z)
          {
            m_partGridDist=(m_gridPosition.m_z-z);
            float m_weightZ=m_iPWeight->CalInterpolationWeight(m_partGridDist);
            float m_weight=m_weightX*m_weightY*m_weightZ;
            int m_currentnode=(z*m_gridSize.m_x*m_gridSize.m_y)+(y*m_gridSize.m_x)+x;
            m_pDensity+=(m_nodes[m_currentnode].m_nodeMass*m_weight)/m_gridSpacingCubed;
          }
        }
      }
      particles[i].SetParticleDensity(m_pDensity);
      m_pVolume=particles[i].m_pMass/particles[i].m_pDensity;
      assert(m_pVolume>0);
      particles[i].SetParticleVolume(m_pVolume);
    }
}

//compute grid velocities
void Grid::ComputeGridForces()
{
  #pragma omp parallel
    for(int i=0; i<m_partamount; ++i)
    {
      ngl::Mat3 m_gDeformationGradient=particles[i].m_PDeformationGradient;
      float m_jacobian=m_gDeformationGradient.determinant();

      ngl::Mat3 m_gPlasticDefG=particles[i].m_plasticDefG;
      float m_jPlastic=m_gPlasticDefG.determinant();
      ngl::Mat3 m_gElasticDefG=particles[i].m_elasticDefG;
      float m_jElastic=m_gElasticDefG.determinant();

      //lambda
      float m_lameFirstParameter=m_poissons*m_youngsmodulus/((1+m_poissons)*(1-2*m_poissons));
      float m_lambda=m_lameFirstParameter*exp(particles[i].m_hardeningCoefficient*(1-m_jPlastic));
      particles[i].SetLambda(m_lambda);

      //mu
      float m_lameSecondParameter=m_youngsmodulus/(2.0*(1.0+m_poissons));
      float m_mu=m_lameSecondParameter*exp(particles[i].m_hardeningCoefficient*(1-m_jPlastic));
      particles[i].SetMu(m_mu);
      ngl::Mat3 m_rotationPD;

      //calculate rotation of the polar decomposition
      #pragma omp critical(polardec)
      {
        m_rotationPD=m_polarDecomposition->ComputePDRotation(m_gElasticDefG);
      }
      //first part of sigma
      float m_tempValue=(2*m_mu)/m_jacobian;
      ngl::Mat3 m_tempMat=m_rotationPD.operator *(-1.0);
      m_tempMat=m_gElasticDefG.operator +(m_tempMat);
      //=2*m_mu/m_jacobian(F-R)
      ngl::Mat3 m_sigma=m_tempMat.operator *(m_tempValue);
      m_tempMat=m_gElasticDefG.transpose();
      m_sigma.operator *=(m_tempMat);

      //second part of sigma
      m_tempValue=(m_lambda/m_jacobian)*(m_jElastic-1)*m_jElastic;
      m_tempMat=m_identity.operator *(m_tempValue);
      m_sigma.operator +=(m_tempMat);

      //m_sigma=2*m_mu/m_jacobian*(m_gElasticDefG-m_rotationPD)*m_gElasticDefG.transpose()+m_lambda/m_jacobian*(m_jElastic-1)*m_jElastic*m_identity;
      float m_pVolume=-1.0*m_jacobian*particles[i].m_pVolume;

      ngl::Vec3 m_gridPosition=(particles[i].m_PPosition-m_gMin)/m_gridSpacing;
      for(int x=particles[i].m_xMin; x<particles[i].m_xMax; ++x)
      {
        float m_partGridDist=(m_gridPosition.m_x-x);
        float m_weightX=m_iPWeight->CalInterpolationWeight(m_partGridDist);
        float m_dWeightX=m_iPWeight->IWGradient(m_partGridDist);
        for(int y=particles[i].m_yMin; y<particles[i].m_yMax; ++y)
        {
          m_partGridDist=(m_gridPosition.m_y-y);
          float m_weightY=m_iPWeight->CalInterpolationWeight(m_partGridDist);
          float m_dWeightY=m_iPWeight->IWGradient(m_partGridDist);
          for(int z=particles[i].m_zMin; z<particles[i].m_zMax; ++z)
          {
            m_partGridDist=(m_gridPosition.m_z-z);
            float m_weightZ=m_iPWeight->CalInterpolationWeight(m_partGridDist);
            float m_dWeightZ=m_iPWeight->IWGradient(m_partGridDist);
            ngl::Vec3 m_weightGradient;
            m_weightGradient.m_x=m_dWeightX*m_weightY*m_weightZ;
            m_weightGradient.m_y=m_weightX*m_dWeightY*m_weightZ;
            m_weightGradient.m_z=m_weightX*m_weightY*m_dWeightZ;
            int m_currentnode=(z*m_gridSize.m_x*m_gridSize.m_y)+(y*m_gridSize.m_x)+x;
            m_tempMat=m_sigma.operator *(m_pVolume);
            #pragma omp critical(nodeupdate2)
            {
              m_nodes[m_currentnode].m_nodeForce+=m_tempMat.operator *(m_weightGradient);
            }
          }
        }
      }
    }
}

void Grid::UpdateGridVelocities()
{
  #pragma omp parallel for
    for(int i=0; i<m_nodesAmount; ++i)
    {
      if(m_nodes[i].m_nodeMass!=0){
        ngl::Vec3 m_tempVec=m_gravity.operator *(m_nodes[i].m_nodeMass);
        //f=f+g*m
        m_nodes[i].m_nodeForce.operator +=(m_tempVec);
        float m_tempValue=m_timeStep/m_nodes[i].m_nodeMass;
        m_tempVec=m_nodes[i].m_nodeForce.operator *(m_tempValue);
        //v=v+f*t/m
        m_nodes[i].m_newNodeVelocity=m_nodes[i].m_nodeVelocity.operator +(m_tempVec);
      }
    }
}

//grid collisions
void Grid::GridCollision()
{
  #pragma omp parallel for
    for(int i=0; i<m_nodesAmount; ++i)
    {
      if(m_nodes[i].m_nodeMass!=0){
        bool m_collision=false;
        ngl::Vec3 m_collisionNormal;
        ngl::Vec3 m_collisionVelocity;
        if(m_nodes[i].m_nodePosition.m_y< 0)
        {
          m_collisionNormal=ngl::Vec3(0,1,0);
          m_collision=true;
        }
        else if(m_nodes[i].m_nodePosition.m_y>= 1)
        {
          m_collisionNormal=ngl::Vec3(0,-1,0);
          m_collision=true;
        }
        if(m_nodes[i].m_nodePosition.m_x>= 0.1)
        {
          m_collisionNormal=ngl::Vec3(-1,0,0);
          m_collision=true;
        }
        else if(m_nodes[i].m_nodePosition.m_x<= -0.05)
        {
          m_collisionNormal=ngl::Vec3(1,0,0);
          m_collision=true;
        }
        if(m_nodes[i].m_nodePosition.m_z>= 0.1)
        {
          m_collisionNormal=ngl::Vec3(0,0,-1);
          m_collision=true;
        }
        else if(m_nodes[i].m_nodePosition.m_z<= -0.1)
        {
          m_collisionNormal=ngl::Vec3(0,0,1);
          m_collision=true;
        }
        if(m_collision==true)
        {
          float m_vn=m_nodes[i].m_newNodeVelocity.dot(m_collisionNormal);
          if(m_vn<0)
          {
            ngl::Vec3 m_velocityTang=m_nodes[i].m_newNodeVelocity-(m_collisionNormal*m_vn);
            float m_vtLength=m_velocityTang.length();
            if(m_vtLength<=(-0.1*m_vn))
            {
              m_collisionVelocity=0;
            }
            else
            {
              m_collisionVelocity=m_velocityTang+(0.1*m_vn*m_velocityTang)/m_vtLength;
            }
            m_nodes[i].m_newNodeVelocity=m_collisionVelocity;
          }
        }
      }
    }
}

void Grid::ExplicitSolver()
{
  #pragma omp parallel for
    for(int i=0; i<m_nodesAmount; ++i)
    {
      m_nodes[i].m_nodeVelocity=m_nodes[i].m_newNodeVelocity;
    }
}

//solve the linear system for semi-implicit integration
void Grid::ImplicitSolver()
{
  //initial values
  ComputeFeReSeHat();
  #pragma omp parallel for
    for(int i=0; i<m_nodesAmount; ++i)
    {
      if(m_nodes[i].m_nodeMass!=0)
      {
        m_nodes[i].m_v=m_nodes[i].m_newNodeVelocity;
        m_nodes[i].m_r=m_nodes[i].m_v;
        m_nodes[i].m_implicitSolved=false;
      }
    }
  ComputeAr();
  #pragma omp parallel for
    for(int i=0; i<m_nodesAmount; ++i)
    {
      if(m_nodes[i].m_nodeMass!=0)
      {
        m_nodes[i].m_r=m_nodes[i].m_v-m_nodes[i].m_Ar;
        m_nodes[i].m_p=m_nodes[i].m_r;
      }
    }
  ComputeAr();
  #pragma omp parallel for
    for(int i=0; i<m_nodesAmount; ++i)
    {
      if(m_nodes[i].m_nodeMass!=0)
      {
        m_nodes[i].m_Ap=m_nodes[i].m_Ar;
      }
    }

  //while all residuals are too big or certain amount of loops
  int j=0;
  do
  {
    m_implicitDone=true;
    #pragma omp parallel for
      for(int i=0; i<m_nodesAmount; ++i)
      {
        if(m_nodes[i].m_nodeMass!=0&&m_nodes[i].m_implicitSolved==false)
        {
          //alpha
          m_nodes[i].m_beta=m_nodes[i].m_r.inner(m_nodes[i].m_Ar);
          m_nodes[i].m_alpha=m_nodes[i].m_Ap.inner(m_nodes[i].m_Ap);
          if(m_nodes[i].m_alpha!=0)
          {
            m_nodes[i].m_alpha=m_nodes[i].m_beta/m_nodes[i].m_alpha;
          }
          else
          {
            m_nodes[i].m_alpha=0;
          }
          //v=v+alpha*p
          m_nodes[i].m_v=m_nodes[i].m_v+m_nodes[i].m_alpha*m_nodes[i].m_p;
          m_nodes[i].m_r=m_nodes[i].m_r-m_nodes[i].m_alpha*m_nodes[i].m_Ap;
          m_nodes[i].m_residual=m_nodes[i].m_r.dot(m_nodes[i].m_r);
          if(m_nodes[i].m_residual<0.0001)
          {
            m_nodes[i].m_implicitSolved=true;
          }
        }
      }

    ComputeAr();
    #pragma omp parallel for
      for(int i=0; i<m_nodesAmount; ++i)
      {
        if(m_nodes[i].m_nodeMass!=0&&m_nodes[i].m_implicitSolved==false)
        {
          //beta
          float m_tempValue=m_nodes[i].m_r.inner(m_nodes[i].m_Ar);
          if(m_nodes[i].m_beta!=0)
          {
            m_nodes[i].m_beta=m_tempValue/m_nodes[i].m_beta;
          }
          else
          {
            m_nodes[i].m_beta=0;
          }
          m_nodes[i].m_p=m_nodes[i].m_r+m_nodes[i].m_beta*m_nodes[i].m_p;
          m_nodes[i].m_Ap=m_nodes[i].m_Ar+m_nodes[i].m_beta*m_nodes[i].m_Ap;
          m_implicitDone=false;
        }
      }
  }while((++j<30&&(m_implicitDone==false)));
  #pragma omp parallel for
    for(int i=0; i<m_nodesAmount; ++i)
    {
      if(m_nodes[i].m_nodeMass!=0)
      {
        m_nodes[i].m_newNodeVelocity=m_nodes[i].m_v;
      }
    }
}

void Grid::ComputeFeReSeHat()
{
    for(int i=0; i<m_partamount; ++i)
    {
      ngl::Mat3 m_pVelocityGradient=0;
      ngl::Mat3 m_tempMat;
      ngl::Vec3 m_gridPosition=(particles[i].m_PPosition-m_gMin)/m_gridSpacing;

      for(int x=particles[i].m_xMin; x<particles[i].m_xMax; ++x)
      {
        float m_partGridDist=(m_gridPosition.m_x-x);
        float m_weightX=m_iPWeight->CalInterpolationWeight(m_partGridDist);
        float m_dWeightX=m_iPWeight->IWGradient(m_partGridDist);
        for(int y=particles[i].m_yMin; y<particles[i].m_yMax; ++y)
        {
          m_partGridDist=(m_gridPosition.m_y-y);
          float m_weightY=m_iPWeight->CalInterpolationWeight(m_partGridDist);
          float m_dWeightY=m_iPWeight->IWGradient(m_partGridDist);
          for(int z=particles[i].m_zMin; z<particles[i].m_zMax; ++z)
          {
            m_partGridDist=(m_gridPosition.m_z-z);
            float m_weightZ=m_iPWeight->CalInterpolationWeight(m_partGridDist);
            float m_dWeightZ=m_iPWeight->IWGradient(m_partGridDist);
            ngl::Vec3 m_weightGradient;
            m_weightGradient.m_x=m_dWeightX*m_weightY*m_weightZ;
            m_weightGradient.m_y=m_weightX*m_dWeightY*m_weightZ;
            m_weightGradient.m_z=m_weightX*m_weightY*m_dWeightZ;
            int m_currentnode=(z*m_gridSize.m_x*m_gridSize.m_y)+(y*m_gridSize.m_x)+x;

            m_tempMat.m_00=m_timeStep*m_nodes[m_currentnode].m_newNodeVelocity.m_x*m_weightGradient.m_x;
            m_tempMat.m_01=m_timeStep*m_nodes[m_currentnode].m_newNodeVelocity.m_x*m_weightGradient.m_y;
            m_tempMat.m_02=m_timeStep*m_nodes[m_currentnode].m_newNodeVelocity.m_x*m_weightGradient.m_z;
            m_tempMat.m_10=m_timeStep*m_nodes[m_currentnode].m_newNodeVelocity.m_y*m_weightGradient.m_x;
            m_tempMat.m_11=m_timeStep*m_nodes[m_currentnode].m_newNodeVelocity.m_y*m_weightGradient.m_y;
            m_tempMat.m_12=m_timeStep*m_nodes[m_currentnode].m_newNodeVelocity.m_y*m_weightGradient.m_z;
            m_tempMat.m_20=m_timeStep*m_nodes[m_currentnode].m_newNodeVelocity.m_z*m_weightGradient.m_x;
            m_tempMat.m_21=m_timeStep*m_nodes[m_currentnode].m_newNodeVelocity.m_z*m_weightGradient.m_y;
            m_tempMat.m_22=m_timeStep*m_nodes[m_currentnode].m_newNodeVelocity.m_z*m_weightGradient.m_z;

            m_pVelocityGradient.operator +=(m_tempMat);
          }
        }
      }
      particles[i].SetVelocityGradient(m_pVelocityGradient);
      m_tempMat=m_pVelocityGradient;
      m_tempMat=m_tempMat.operator +=(m_tempMat.identity());
      ngl::Mat3 m_gElasticDefG=m_tempMat.operator *(particles[i].m_elasticDefG);
      particles[i].SetFeHat(m_gElasticDefG);
      ngl::Mat3 m_rotationPD=m_polarDecomposition->ComputeRandS(particles[i].m_FeHat);
      ngl::Mat3 m_sPD=m_polarDecomposition->GetPDS();
      particles[i].SetReHat(m_rotationPD);
      particles[i].SetSeHat(m_sPD);
    }
}

void Grid::ComputeAr()
{
  //df=sum(Vp*Ap*FEp^T*weight gradient)
  //Ap=2*mu(dFE-dRE)+lambda*JF^-T(JF^-T:dF)+lambda*(J-1)d(JF^-T)
  #pragma omp parallel for
  //compute dF
    for(int i=0; i<m_partamount; ++i)
    {
      ngl::Mat3 m_dElasticDefG=0;
      ngl::Vec3 m_gridPosition=(particles[i].m_PPosition-m_gMin)/m_gridSpacing;
      for(int x=particles[i].m_xMin; x<particles[i].m_xMax; ++x)
      {
        float m_partGridDist=(m_gridPosition.m_x-x);
        float m_weightX=m_iPWeight->CalInterpolationWeight(m_partGridDist);
        float m_dWeightX=m_iPWeight->IWGradient(m_partGridDist);
        for(int y=particles[i].m_yMin; y<particles[i].m_yMax; ++y)
        {
          m_partGridDist=(m_gridPosition.m_y-y);
          float m_weightY=m_iPWeight->CalInterpolationWeight(m_partGridDist);
          float m_dWeightY=m_iPWeight->IWGradient(m_partGridDist);
          for(int z=particles[i].m_zMin; z<particles[i].m_zMax; ++z)
          {
            m_partGridDist=(m_gridPosition.m_z-z);
            float m_weightZ=m_iPWeight->CalInterpolationWeight(m_partGridDist);
            float m_dWeightZ=m_iPWeight->IWGradient(m_partGridDist);
            ngl::Vec3 m_weightGradient;
            m_weightGradient.m_x=m_dWeightX*m_weightY*m_weightZ;
            m_weightGradient.m_y=m_weightX*m_dWeightY*m_weightZ;
            m_weightGradient.m_z=m_weightX*m_weightY*m_dWeightZ;
            int m_currentnode=(z*m_gridSize.m_x*m_gridSize.m_y)+(y*m_gridSize.m_x)+x;
            ngl::Mat3 m_tempMat;
            m_tempMat.m_00=m_nodes[m_currentnode].m_r.m_x*m_weightGradient.m_x;
            m_tempMat.m_01=m_nodes[m_currentnode].m_r.m_x*m_weightGradient.m_y;
            m_tempMat.m_02=m_nodes[m_currentnode].m_r.m_x*m_weightGradient.m_z;
            m_tempMat.m_10=m_nodes[m_currentnode].m_r.m_y*m_weightGradient.m_x;
            m_tempMat.m_11=m_nodes[m_currentnode].m_r.m_y*m_weightGradient.m_y;
            m_tempMat.m_12=m_nodes[m_currentnode].m_r.m_y*m_weightGradient.m_z;
            m_tempMat.m_20=m_nodes[m_currentnode].m_r.m_z*m_weightGradient.m_x;
            m_tempMat.m_21=m_nodes[m_currentnode].m_r.m_z*m_weightGradient.m_y;
            m_tempMat.m_22=m_nodes[m_currentnode].m_r.m_z*m_weightGradient.m_z;
            m_tempMat.operator *=(m_timeStep);
            m_dElasticDefG+=m_tempMat;
          }
        }
      }
      m_dElasticDefG.operator *=(particles[i].m_elasticDefG);
      particles[i].SetdElasticDefG(m_dElasticDefG);
    }
  //calculate df
  for(int i=0; i<m_partamount; ++i)
  {
    ngl::Mat3 m_rotationPD=particles[i].m_ReHat;
    ngl::Mat3 m_sPD=particles[i].m_SeHat;
    //calculate dR=R^T*dF-dF^TR
    m_LH=m_rotationPD;

    m_LH.transpose();
    m_LH.operator *=(particles[i].m_dElasticDefG);
    ngl::Mat3 m_tempMat;
    m_tempMat=particles[i].m_dElasticDefG;
    m_tempMat.transpose();
    m_tempMat.operator *=(m_rotationPD);
    m_tempMat.operator *=(-1.0);
    m_LH.operator +=(m_tempMat);

    //the S matrix
    m_tempMat.m_00=m_sPD.m_00+m_sPD.m_11;
    m_tempMat.m_01=m_sPD.m_21;
    m_tempMat.m_02=-1.0*m_sPD.m_20;
    m_tempMat.m_10=m_sPD.m_21;
    m_tempMat.m_11=m_sPD.m_00+m_sPD.m_22;
    m_tempMat.m_12=m_sPD.m_10;
    m_tempMat.m_20=-1.0*m_sPD.m_20;
    m_tempMat.m_21=m_sPD.m_10;
    m_tempMat.m_22=m_sPD.m_11+m_sPD.m_22;


    ngl::Vec3 m_tempVec;
    m_tempVec.m_x=m_LH.m_01;
    m_tempVec.m_y=m_LH.m_02;
    m_tempVec.m_z=m_LH.m_12;

    m_tempMat.inverse();
    m_tempVec=m_tempMat.operator *(m_tempVec);

    m_tempMat=0.0;
    m_tempMat.m_01=m_tempVec.m_x;
    m_tempMat.m_02=m_tempVec.m_y;
    m_tempMat.m_10=-1.0*m_tempVec.m_x;
    m_tempMat.m_12=m_tempVec.m_z;
    m_tempMat.m_20=-1.0*m_tempVec.m_y;
    m_tempMat.m_21=-1.0*m_tempVec.m_z;

    m_dR=m_rotationPD.operator *(m_tempMat);

    particles[i].SetdR(m_dR);
    //dJF_invTrans
    ngl::Mat3 m_dElasticDefG=particles[i].m_dElasticDefG;
    ngl::Mat3 m_gElasticDefG=particles[i].m_FeHat;
    m_dJFinvT.m_00=m_gElasticDefG.m_11*m_dElasticDefG.m_22-m_gElasticDefG.m_21*m_dElasticDefG.m_12-m_gElasticDefG.m_12*m_dElasticDefG.m_21+m_gElasticDefG.m_22*m_dElasticDefG.m_11;
    m_dJFinvT.m_01=m_gElasticDefG.m_20*m_dElasticDefG.m_12-m_gElasticDefG.m_10*m_dElasticDefG.m_22+m_gElasticDefG.m_12*m_dElasticDefG.m_20-m_gElasticDefG.m_22*m_dElasticDefG.m_10;
    m_dJFinvT.m_02=m_gElasticDefG.m_10*m_dElasticDefG.m_21-m_gElasticDefG.m_20*m_dElasticDefG.m_11-m_gElasticDefG.m_11*m_dElasticDefG.m_20+m_gElasticDefG.m_21*m_dElasticDefG.m_10;
    m_dJFinvT.m_10=m_gElasticDefG.m_21*m_dElasticDefG.m_02-m_gElasticDefG.m_01*m_dElasticDefG.m_22+m_gElasticDefG.m_02*m_dElasticDefG.m_21-m_gElasticDefG.m_22*m_dElasticDefG.m_01;
    m_dJFinvT.m_11=m_gElasticDefG.m_00*m_dElasticDefG.m_22-m_gElasticDefG.m_20*m_dElasticDefG.m_02-m_gElasticDefG.m_02*m_dElasticDefG.m_20+m_gElasticDefG.m_22*m_dElasticDefG.m_00;
    m_dJFinvT.m_12=m_gElasticDefG.m_20*m_dElasticDefG.m_01-m_gElasticDefG.m_00*m_dElasticDefG.m_21+m_gElasticDefG.m_01*m_dElasticDefG.m_20-m_gElasticDefG.m_21*m_dElasticDefG.m_00;
    m_dJFinvT.m_20=m_gElasticDefG.m_01*m_dElasticDefG.m_12-m_gElasticDefG.m_11*m_dElasticDefG.m_02-m_gElasticDefG.m_02*m_dElasticDefG.m_11+m_gElasticDefG.m_12*m_dElasticDefG.m_01;
    m_dJFinvT.m_21=m_gElasticDefG.m_10*m_dElasticDefG.m_02-m_gElasticDefG.m_00*m_dElasticDefG.m_12+m_gElasticDefG.m_02*m_dElasticDefG.m_10-m_gElasticDefG.m_12*m_dElasticDefG.m_00;
    m_dJFinvT.m_22=m_gElasticDefG.m_00*m_dElasticDefG.m_11-m_gElasticDefG.m_10*m_dElasticDefG.m_01-m_gElasticDefG.m_01*m_dElasticDefG.m_10+m_gElasticDefG.m_11*m_dElasticDefG.m_00;

    //mat3 JFe_invTrans = mat3::cofactor( Fe );
    m_JFinvT.m_00=m_gElasticDefG.m_11*m_gElasticDefG.m_22-m_gElasticDefG.m_21*m_gElasticDefG.m_12;
    m_JFinvT.m_01=m_gElasticDefG.m_20*m_dElasticDefG.m_12-m_gElasticDefG.m_10*m_dElasticDefG.m_22;
    m_JFinvT.m_02=m_gElasticDefG.m_10*m_dElasticDefG.m_21-m_gElasticDefG.m_20*m_dElasticDefG.m_11;
    m_JFinvT.m_10=m_gElasticDefG.m_21*m_dElasticDefG.m_02-m_gElasticDefG.m_01*m_dElasticDefG.m_22;
    m_JFinvT.m_11=m_gElasticDefG.m_00*m_dElasticDefG.m_22-m_gElasticDefG.m_20*m_dElasticDefG.m_02;
    m_JFinvT.m_12=m_gElasticDefG.m_20*m_dElasticDefG.m_01-m_gElasticDefG.m_00*m_dElasticDefG.m_21;
    m_JFinvT.m_20=m_gElasticDefG.m_01*m_dElasticDefG.m_12-m_gElasticDefG.m_11*m_dElasticDefG.m_02;
    m_JFinvT.m_21=m_gElasticDefG.m_10*m_dElasticDefG.m_02-m_gElasticDefG.m_00*m_dElasticDefG.m_12;
    m_JFinvT.m_22=m_gElasticDefG.m_00*m_dElasticDefG.m_11-m_gElasticDefG.m_10*m_dElasticDefG.m_01;

    //Ap=2*mu(dFE-dRE)+lambda*JF^-T(JF^-T:dF)+lambda*(J-1)d(JF^-T)
    m_Ap=m_dR.operator *(-1.0);
    m_Ap.operator +=(m_dElasticDefG);
    float m_mu=2.0*particles[i].m_mu;
    m_Ap.operator *=(m_mu);
    float m_lambda=particles[i].m_lambda;

    float m_tempValue=m_JFinvT.m_00*m_dElasticDefG.m_00+m_JFinvT.m_01*m_dElasticDefG.m_01+m_JFinvT.m_02*m_dElasticDefG.m_02
            +m_JFinvT.m_10*m_dElasticDefG.m_10+m_JFinvT.m_11*m_dElasticDefG.m_11+m_JFinvT.m_12*m_dElasticDefG.m_12
            +m_JFinvT.m_20*m_dElasticDefG.m_20+m_JFinvT.m_21*m_dElasticDefG.m_21+m_JFinvT.m_22*m_dElasticDefG.m_22;   
    m_tempValue*=m_lambda;
    m_tempMat=m_JFinvT.operator *(m_tempValue);
    m_Ap.operator +=(m_tempMat);
    //lambda*(J-1)d(JF^-T)
    float m_jElastic=m_gElasticDefG.determinant();
    m_tempValue=m_lambda*(m_jElastic-1);
    m_tempMat=m_dJFinvT.operator *(m_tempValue);
    m_Ap.operator +=(m_tempMat);
    particles[i].SetAp(m_Ap);

    //df=sum(Vp*Ap*FEp^T*weight gradient)
    ngl::Vec3 m_gridPosition=(particles[i].m_PPosition-m_gMin)/m_gridSpacing;

    for(int x=particles[i].m_xMin; x<particles[i].m_xMax; ++x)
    {
      float m_partGridDist=(m_gridPosition.m_x-x);
      float m_weightX=m_iPWeight->CalInterpolationWeight(m_partGridDist);
      float m_dWeightX=m_iPWeight->IWGradient(m_partGridDist);
      for(int y=particles[i].m_yMin; y<particles[i].m_yMax; ++y)
      {
        m_partGridDist=(m_gridPosition.m_y-y);
        float m_weightY=m_iPWeight->CalInterpolationWeight(m_partGridDist);
        float m_dWeightY=m_iPWeight->IWGradient(m_partGridDist);
        for(int z=particles[i].m_zMin; z<particles[i].m_zMax; ++z)
        {
          int m_currentnode=(z*m_gridSize.m_x*m_gridSize.m_y)+(y*m_gridSize.m_x)+x;
          if(m_nodes[m_currentnode].m_implicitSolved==false)
          {
            m_partGridDist=(m_gridPosition.m_z-z);
            float m_weightZ=m_iPWeight->CalInterpolationWeight(m_partGridDist);
            float m_dWeightZ=m_iPWeight->IWGradient(m_partGridDist);

            ngl::Vec3 m_weightGradient;
            m_weightGradient.m_x=m_dWeightX*m_weightY*m_weightZ;
            m_weightGradient.m_y=m_weightX*m_dWeightY*m_weightZ;
            m_weightGradient.m_z=m_weightX*m_weightY*m_dWeightZ;
            m_tempMat=m_gElasticDefG.transpose();
            m_tempMat=m_Ap.operator *(m_tempMat);

            m_tempMat.operator *=(particles[i].m_pVolume);
            m_tempVec=m_tempMat.operator *(m_weightGradient);
            m_tempVec.operator *=(-1.0);
            #pragma omp critical(nodeupdate3)
            {
              m_nodes[m_currentnode].m_df+=m_tempVec;
            }
          }
        }
      }
    }
  }
  //Ar
  #pragma omp parallel for
    for(int i=0; i<m_nodesAmount; ++i)
    {
      if(m_nodes[i].m_nodeMass!=0&&m_nodes[i].m_implicitSolved==false)
      {
        m_nodes[i].m_Ar=m_nodes[i].m_r-m_nodes[i].m_df.operator *(0.5*m_timeStep/m_nodes[i].m_nodeMass);
      }
    }
}

void Grid::CalculateDeformationGradient()
{
  for(int i=0; i<m_partamount; ++i)
  {
    ngl::Mat3 m_pVelocityGradient=0;
    ngl::Mat3 m_tempMat;
    ngl::Vec3 m_gridPosition=(particles[i].m_PPosition-m_gMin)/m_gridSpacing;
    for(int x=particles[i].m_xMin; x<particles[i].m_xMax; ++x)
    {
      float m_partGridDist=(m_gridPosition.m_x-x);
      float m_weightX=m_iPWeight->CalInterpolationWeight(m_partGridDist);
      float m_dWeightX=m_iPWeight->IWGradient(m_partGridDist);
      for(int y=particles[i].m_yMin; y<particles[i].m_yMax; ++y)
      {
        m_partGridDist=(m_gridPosition.m_y-y);
        float m_weightY=m_iPWeight->CalInterpolationWeight(m_partGridDist);
        float m_dWeightY=m_iPWeight->IWGradient(m_partGridDist);
        for(int z=particles[i].m_zMin; z<particles[i].m_zMax; ++z)
        {
          m_partGridDist=(m_gridPosition.m_z-z);
          float m_weightZ=m_iPWeight->CalInterpolationWeight(m_partGridDist);
          float m_dWeightZ=m_iPWeight->IWGradient(m_partGridDist);

          ngl::Vec3 m_weightGradient;
          m_weightGradient.m_x=m_dWeightX*m_weightY*m_weightZ;
          m_weightGradient.m_y=m_weightX*m_dWeightY*m_weightZ;
          m_weightGradient.m_z=m_weightX*m_weightY*m_dWeightZ;
          int m_currentnode=(z*m_gridSize.m_x*m_gridSize.m_y)+(y*m_gridSize.m_x)+x;

          m_tempMat.m_00=m_nodes[m_currentnode].m_newNodeVelocity.m_x*m_weightGradient.m_x;
          m_tempMat.m_01=m_nodes[m_currentnode].m_newNodeVelocity.m_x*m_weightGradient.m_y;
          m_tempMat.m_02=m_nodes[m_currentnode].m_newNodeVelocity.m_x*m_weightGradient.m_z;
          m_tempMat.m_10=m_nodes[m_currentnode].m_newNodeVelocity.m_y*m_weightGradient.m_x;
          m_tempMat.m_11=m_nodes[m_currentnode].m_newNodeVelocity.m_y*m_weightGradient.m_y;
          m_tempMat.m_12=m_nodes[m_currentnode].m_newNodeVelocity.m_y*m_weightGradient.m_z;
          m_tempMat.m_20=m_nodes[m_currentnode].m_newNodeVelocity.m_z*m_weightGradient.m_x;
          m_tempMat.m_21=m_nodes[m_currentnode].m_newNodeVelocity.m_z*m_weightGradient.m_y;
          m_tempMat.m_22=m_nodes[m_currentnode].m_newNodeVelocity.m_z*m_weightGradient.m_z;

          m_pVelocityGradient.operator +=(m_tempMat);
        }
      }
    }
    particles[i].SetVelocityGradient(m_pVelocityGradient);
    m_tempMat=m_pVelocityGradient.operator *(m_timeStep);
    m_tempMat.operator +=(m_identity);

    ngl::Mat3 m_gElasticDefG=m_tempMat.operator *(particles[i].m_elasticDefG);
    ngl::Mat3 m_gPlasticDefG=particles[i].m_plasticDefG;
    ngl::Mat3 m_gDeformationGradient=m_gElasticDefG.operator *(m_gPlasticDefG);
    particles[i].SetDeformationGradient(m_gDeformationGradient);
    m_outsideRange=m_polarDecomposition->ComputeSVD(m_gElasticDefG);
    //clamp sigma and get Fep and Fpp, input is m_gElasic and pdefmoration gradient
    if(m_outsideRange==1)
    {
       m_polarDecomposition->ComputeElasticAndPlasticDefGrad(m_gElasticDefG, m_gDeformationGradient);
       m_gElasticDefG=m_polarDecomposition->GetElasticDG();
       m_gPlasticDefG=m_polarDecomposition->GetPlasticDG();
       m_gPlasticDefG.operator *=(m_gDeformationGradient);
       particles[i].SetElasticAndPlasticDefG(m_gElasticDefG, m_gPlasticDefG);
    }
    else
    {
      particles[i].SetElasticAndPlasticDefG(m_gElasticDefG, m_gPlasticDefG);
    }
  }
}

void Grid::UpdateParticleVelocities()
{
  //calculate vpic
  //weighted sum of the new velocity of all grid nodes
  //calculate vflip
  //current particle velocity + the weighted sum of the new grid velocity minus the old grid velocity
  //calculate and set new
  m_newGridSize=m_gridSize;
  m_newGMin= m_newGMax=particles[0].m_PPosition;

    for(int i=0; i<m_partamount; ++i)
    {
      ngl::Vec3 m_picVelocity=0;
      ngl::Vec3 m_flipVelocity=0;
      ngl::Vec3 m_gridPosition=(particles[i].m_PPosition-m_gMin)/m_gridSpacing;
      for(int x=particles[i].m_xMin; x<particles[i].m_xMax; ++x)
      {
        float m_partGridDist=(m_gridPosition.m_x-x);//*m_gridSpacing;
        float m_weightX=m_iPWeight->CalInterpolationWeight(m_partGridDist);
        for(int y=particles[i].m_yMin; y<particles[i].m_yMax; ++y)
        {
          m_partGridDist=(m_gridPosition.m_y-y);//*m_gridSpacing;
          float m_weightY=m_iPWeight->CalInterpolationWeight(m_partGridDist);
          for(int z=particles[i].m_zMin; z<particles[i].m_zMax; ++z)
          {
            m_partGridDist=(m_gridPosition.m_z-z);//*m_gridSpacing;
            float m_weightZ=m_iPWeight->CalInterpolationWeight(m_partGridDist);
            float m_weight=m_weightX*m_weightY*m_weightZ;
            float m_currentnode=(z*m_gridSize.m_x*m_gridSize.m_y)+(y*m_gridSize.m_x)+x;
            m_picVelocity+=m_nodes[m_currentnode].m_newNodeVelocity*m_weight;
            m_flipVelocity+=(m_nodes[m_currentnode].m_newNodeVelocity-m_nodes[m_currentnode].m_nodeVelocity)*m_weight;
          }
        }
      }
      m_flipVelocity+=particles[i].m_pVelocity;
      m_pVelocity=(1-0.95)*m_picVelocity+0.95*m_flipVelocity;

      particles[i].SetParticleVelocity(m_pVelocity);
      particles[i].ParticleCollision();
      particles[i].UpdateParticlePosition(m_timeStep);
      if(particles[i].m_PPosition.m_x<m_newGMin.m_x&&particles[i].m_PPosition.m_x>-1)
      {
        m_newGMin.m_x=particles[i].m_PPosition.m_x;
      }
      else if(particles[i].m_PPosition.m_x>m_newGMax.m_x&&particles[i].m_PPosition.m_x<1)
      {
        m_newGMax.m_x=particles[i].m_PPosition.m_x;
      }
      if(particles[i].m_PPosition.m_y<m_newGMin.m_y&&particles[i].m_PPosition.m_y>-0.1)
      {
        m_newGMin.m_y=particles[i].m_PPosition.m_y;
      }
      else if(particles[i].m_PPosition.m_y>m_newGMax.m_y&&particles[i].m_PPosition.m_y<1)
      {
        m_newGMax.m_y=particles[i].m_PPosition.m_y;
      }
      if(particles[i].m_PPosition.m_z<m_newGMin.m_z&&particles[i].m_PPosition.m_z>-1)
      {
        m_newGMin.m_z=particles[i].m_PPosition.m_z;
      }
      else if(particles[i].m_PPosition.m_z>m_newGMax.m_z&&particles[i].m_PPosition.m_z<1)
      {
        m_newGMax.m_z=particles[i].m_PPosition.m_z;
      }
      if(particles[i].m_PPosition.m_x<-1||particles[i].m_PPosition.m_x>1||
         particles[i].m_PPosition.m_y<-0.1||particles[i].m_PPosition.m_y>1||
         particles[i].m_PPosition.m_z<-1||particles[i].m_PPosition.m_z>1)
      {
        particles[i].SetParticlePosition(ngl::Vec3(0.0,0.0,0.0));
      }
    }
  m_newGMax=m_newGMax+ngl::Vec3(2*m_gridSpacing,2*m_gridSpacing,2*m_gridSpacing);
  m_gMin=m_newGMin=m_newGMin-ngl::Vec3(2*m_gridSpacing,2*m_gridSpacing,2*m_gridSpacing);
  m_newGridSize=(m_newGMax-m_newGMin)/m_gridSpacing+ngl::Vec3 (1,1,1);
  m_newGridSize.m_x=ceil(m_newGridSize.m_x);
  m_newGridSize.m_y=ceil(m_newGridSize.m_y);
  m_newGridSize.m_z=ceil(m_newGridSize.m_z);
  m_newNodesAmount=(m_newGridSize.m_x)*(m_newGridSize.m_y)*(m_newGridSize.m_z);
  //write out the pointcloud file every 100 timesteps
  if(m_frameUpdate==100)
  {
    m_expPart->MakeFile(m_frameNumber);
    for(int i=0; i<m_partamount; ++i)
    {
      m_expPart->InputParticleData(particles[i].m_PPosition);
    }
    m_expPart->CloseFile();
    m_frameNumber+=1;
    m_frameUpdate=0;
  }
  m_frameUpdate+=1;
}
