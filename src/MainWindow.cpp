#include "MainWindow.h"
#include "ui_MainWindow.h"

//----------------------------------------------------------------------------------------------------------------------
MainWindow::MainWindow(QWidget *parent) :QMainWindow(parent), m_ui(new Ui::MainWindow)
{
  // create an openGL format and pass to the new GLWidget
  QGLFormat format;
  format.setVersion(3,3);
  format.setProfile( QGLFormat::CoreProfile);
  m_timer=0.2;

	// setup the user interface
	m_ui->setupUi(this);
	// create our GL window for drawing the spring

	// add the glWindow to the UI

  m_MPM=new MPM(1340, m_timer);//m_ui->m_pAmount->value());
  m_gl=new  GLWindow(format,m_timer,this, m_MPM);
  m_ui->s_mainWindowGridLayout->addWidget(m_gl,0,0,4,1);
}

//----------------------------------------------------------------------------------------------------------------------
MainWindow::~MainWindow()
{
  delete m_ui;
  delete m_MPM;
}

//----------------------------------------------------------------------------------------------------------------------
void MainWindow::toggleSim(bool _s)
{
	if(_s == true)
	{
		m_gl->startSimTimer();
	}
	else
	{
		m_gl->stopSimTimer();
	}
}
