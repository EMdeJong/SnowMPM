#include "GLWindow.h"
#include <iostream>
#include <ngl/Camera.h>
#include <ngl/Colour.h>
#include <ngl/Light.h>
#include <ngl/Mat4.h>
#include <ngl/Transformation.h>
#include <ngl/Material.h>
#include <ngl/NGLInit.h>
#include <ngl/Obj.h>
#include <ngl/Random.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/ShaderLib.h>
#include <boost/foreach.hpp>



//----------------------------------------------------------------------------------------------------------------------
GLWindow::GLWindow(const QGLFormat _format, float _timer, QWidget *_parent, MPM *_mpm) : QGLWidget( _format, _parent )
{
  // set this widget to have the initial keyboard focus
  setFocus();
  // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
  this->resize(_parent->size());
  // Now set the initial GLWindow attributes to default values
  // Roate is false
  m_rotate=false;
  // mouse rotation values set to 0
  m_spinXFace=0;
  m_spinYFace=0;
	m_timerValue=_timer;
	m_MPM=_mpm;
}

GLWindow::~GLWindow()
{
	ngl::NGLInit *init = ngl::NGLInit::instance();
	init->NGLQuit();
}
// This virtual function is called once before the first call to paintGL() or resizeGL(),
//and then once whenever the widget has been assigned a new QGLContext.
// This function should set up any required OpenGL context rendering flags, defining display lists, etc.

//----------------------------------------------------------------------------------------------------------------------
void GLWindow::initializeGL()
{
  ngl::NGLInit::instance();
  glClearColor(0.4f, 0.4f, 0.4f, 1.0f);			   // Grey Background
  // enable depth testing for drawing
 glEnable(GL_DEPTH_TEST);
 // Now we will create a basic Camera from the graphics library
 // This is a static camera so it only needs to be set once
 // First create Values for the camera position
 ngl::Vec3 from(0,2.5,15);
 ngl::Vec3 to(0,2,0);
 ngl::Vec3 up(0,1,0);
 ngl::NGLInit::instance();
 m_cam= new ngl::Camera(from,to,up);
 // set the shape using FOV 45 Aspect Ratio based on Width and Height
 // The final two are near and far clipping planes of 0.5 and 10
 m_cam->setShape(45,(float)720.0/576.0,0.5,150);




 // now to load the shader andstuck on different screen set the values
 // grab an instance of shader manager
 ngl::ShaderLib *shader=ngl::ShaderLib::instance();
 // we are creating a shader called Phong
 shader->createShaderProgram("Phong");
 // now we are going to create empty shaders for Frag and Vert
 shader->attachShader("PhongVertex",ngl::VERTEX);
 shader->attachShader("PhongFragment",ngl::FRAGMENT);
 // attach the source
 shader->loadShaderSource("PhongVertex","shaders/PhongVertex.glsl");
 shader->loadShaderSource("PhongFragment","shaders/PhongFragment.glsl");
 // compile the shaders
 shader->compileShader("PhongVertex");
 shader->compileShader("PhongFragment");
 // add them to the program
 shader->attachShaderToProgram("Phong","PhongVertex");
 shader->attachShaderToProgram("Phong","PhongFragment");

 // now we have associated this data we can link the shader
 shader->linkProgramObject("Phong");
 // and make it active ready to load values
 (*shader)["Phong"]->use();
 shader->setShaderParam1i("Normalize",1);

 // now pass the modelView and projection values to the shader
 // the shader will use the currently active material and light0 so set them
 ngl::Material m(ngl::SILVER);
 m.loadToShader("material");
 ngl::Light light(ngl::Vec3(0,0,2),ngl::Colour(1,1,1,1),ngl::Colour(1,1,1,1),ngl::POINTLIGHT);
 // now create our light this is done after the camera so we can pass the
 // transpose of the projection matrix to the light to do correct eye space
 // transformations
 ngl::Mat4 iv=m_cam->getViewMatrix();
 iv.transpose();
 light.setTransform(iv);
 light.setAttenuation(1,0,0);
 light.enable();
 // load these values to the shader as well
 light.loadToShader("light");
 shader->createShaderProgram("Colour");

 shader->attachShader("ColourVertex",ngl::VERTEX);
 shader->attachShader("ColourFragment",ngl::FRAGMENT);
 shader->loadShaderSource("ColourVertex","shaders/ColourVertex.glsl");
 shader->loadShaderSource("ColourFragment","shaders/ColourFragment.glsl");

 shader->compileShader("ColourVertex");
 shader->compileShader("ColourFragment");
 shader->attachShaderToProgram("Colour","ColourVertex");
 shader->attachShaderToProgram("Colour","ColourFragment");

 shader->bindAttribute("Colour",0,"inVert");

 shader->linkProgramObject("Colour");

 ngl::VAOPrimitives *prim=ngl::VAOPrimitives::instance();
 prim->createSphere("sphere",0.05,10);
 prim->createTrianglePlane("ground", 10, 10, 1, 1,(0,1,0));
 startSimTimer();

}

//----------------------------------------------------------------------------------------------------------------------
//This virtual function is called whenever the widget has been resized.
// The new size is passed in width and height.
void GLWindow::resizeGL( int _w, int _h)
{
  glViewport(0,0,_w,_h);
  m_cam->setShape(45,(float)_w/_h,0.5,150);
}

void GLWindow::loadMatricesToShader()
{
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();
  (*shader)["Phong"]->use();

  ngl::Mat4 MV;
  ngl::Mat4 MVP;
  ngl::Mat3 normalMatrix;
  ngl::Mat4 M;
  M=m_transform.getMatrix();
  MV=M*m_mouseGlobalTX*m_cam->getViewMatrix() ;
  MVP=MV*m_cam->getProjectionMatrix();
  normalMatrix=MV;
  normalMatrix.inverse();
  shader->setShaderParamFromMat4("MV",MV);
  shader->setShaderParamFromMat4("MVP",MVP);
  shader->setShaderParamFromMat3("normalMatrix",normalMatrix);
  shader->setShaderParamFromMat4("M",M);
}

void GLWindow::loadMatricesToColourShader()
{
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();
  (*shader)["Colour"]->use();
  ngl::Mat4 MVP;
  MVP=m_transform.getMatrix()*m_mouseGlobalTX*m_cam->getVPMatrix() ;
  shader->setShaderParamFromMat4("MVP",MVP);

}


//----------------------------------------------------------------------------------------------------------------------
//This virtual function is called whenever the widget needs to be painted.
// this is our main drawing routine
void GLWindow::paintGL()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void GLWindow::timerEvent( QTimerEvent *_event)
{
	updateGL();
	m_MPM->UpdateGridOperations();
}

void GLWindow::startSimTimer()
{
	m_timer=startTimer(m_timerValue);
}

void GLWindow::stopSimTimer()
{
 killTimer(m_timer);
}


