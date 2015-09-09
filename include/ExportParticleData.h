#ifndef EXPORTPARTICLEDATA_H
#define EXPORTPARTICLEDATA_H

#include <sstream>
#include <iomanip>
#include <iostream>
#include <pointcloud.h>

#include <ngl/Vec3.h>


struct xyz	 { float x, y, z; };


class ExportParticleData
{
public:
  ExportParticleData();
  //creates the new file with the current frameNumber
  void MakeFile(int _frameNumber);
  //adds the particle position of the current particle to the pointcloud file
  void InputParticleData(ngl::Vec3 _partPos);
  //closes the current pointcloud file when all particles are added
  void CloseFile();

private:
  struct xyz position;
  struct xyz normal;
  int m_framenumber;
  int m_nvars;
  int m_test;
  PtcPointCloud outptc;

  const char* m_outputfileNameCharTemp;
  char *m_outputfileNameChar;
  char **m_vartypes;
  char **m_varnames;
  float w2e[16];
  float w2n[16];
  float m_format[3];
  ngl::Vec3 m_curPartPosition;
  std::stringstream m_outputfileName;
  std::string m_temp;
};

#endif // EXPORTPARTICLEDATA_H
