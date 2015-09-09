#include "ExportParticleData.h"

ExportParticleData::ExportParticleData()
{
  m_nvars=0;
  m_vartypes= NULL;
  m_varnames= NULL;
}

void ExportParticleData::MakeFile(int _frameNumber)
{
	m_framenumber=_frameNumber;
	//reset the filename
	m_outputfileName.str("");
	m_outputfileName << "SnowF" <<std::setw(4)<<std::setfill('0')<<m_framenumber;
	m_temp=m_outputfileName.str()+".ptc";
	m_outputfileNameCharTemp=m_temp.c_str();
	m_outputfileNameChar=(char*)m_outputfileNameCharTemp;
	outptc = PtcCreatePointCloudFile(m_outputfileNameChar,
							 m_nvars, m_vartypes, m_varnames,
							 w2e, w2n, m_format);
}

void ExportParticleData::InputParticleData(ngl::Vec3 _partPos)
{
  position.x=_partPos.m_x;
  position.y=_partPos.m_y;
  position.z=_partPos.m_z;

  PtcWriteDataPoint(outptc,
          (float*)&position,
          (float*)&normal,
          0.2,
          0);
}

void ExportParticleData::CloseFile()
{
  PtcFinishPointCloudFile(outptc);
}
