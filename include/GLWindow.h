#ifndef __GL_WINDOW_H__
#define __GL_WINDOW_H__


#include <ngl/Camera.h>
#include <ngl/Colour.h>
#include <ngl/Transformation.h>
#include <ngl/Types.h>
#include "MPM.h"

// must be included after our stuff becuase GLEW needs to be first
#include <QtOpenGL>

//----------------------------------------------------------------------------------------------------------------------
/// @file GLWindow.h
/// @brief a basic Qt GL window class for ngl demos
/// @author Jonathan Macey
/// @version 1.0
/// @date 10/10/10
/// Revision History :
/// Initial Version 10/10/10 (Binary day ;-0 )
/// @class GLWindow
/// @brief our main glwindow widget for NGL applications all drawing elements are
/// put in this file
//----------------------------------------------------------------------------------------------------------------------
class GLWindow : public QGLWidget
{
Q_OBJECT        // must include this if you use Qt signals/slots
public :
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief Constructor for GLWindow
	/// @param[in] _timer the time value for simulation updates
  /// @param [in] _parent the parent window to create the GL context in
	//----------------------------------------------------------------------------------------------------------------------
		GLWindow(const QGLFormat _format, float _timer,QWidget *_parent, MPM*_mpm);
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief dtor
	//----------------------------------------------------------------------------------------------------------------------
	~GLWindow();
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief set the spring to use
	/// @param[in] _s the spring
	//----------------------------------------------------------------------------------------------------------------------
	//inline void setSpring(RK4Spring *_s){m_spring=_s;}
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief start the simulation timer
	//----------------------------------------------------------------------------------------------------------------------
	void startSimTimer();
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief stop the simulation timer
	//----------------------------------------------------------------------------------------------------------------------
	void stopSimTimer();

public slots :
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief set the timer duration value
	/// @param[in] _v the timer value in ms
	//----------------------------------------------------------------------------------------------------------------------
	inline void setTimerDuration(int _v){m_timerValue=_v;}



private :
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief the timer
	//----------------------------------------------------------------------------------------------------------------------
	float m_timer;
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief used to store the x rotation mouse value
	//----------------------------------------------------------------------------------------------------------------------
  int m_spinXFace;
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief used to store the y rotation mouse value
	//----------------------------------------------------------------------------------------------------------------------
  int m_spinYFace;
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief flag to indicate if the mouse button is pressed when dragging
	//----------------------------------------------------------------------------------------------------------------------
  bool m_rotate;
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief the previous x mouse value
	//----------------------------------------------------------------------------------------------------------------------
  int m_origX;
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief the previous y mouse value
	//----------------------------------------------------------------------------------------------------------------------
  int m_origY;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief used to store the global mouse transforms
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Mat4 m_mouseGlobalTX;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Our Camera
	//----------------------------------------------------------------------------------------------------------------------
  ngl::Camera *m_cam;
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief transformation stack for the gl transformations etc
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief our transformation stack used for drawing
	//----------------------------------------------------------------------------------------------------------------------
	ngl::Transformation m_transform;
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief our spring pointer
	//----------------------------------------------------------------------------------------------------------------------
	//RK4Spring *m_spring;
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief the timer value in ms
	//----------------------------------------------------------------------------------------------------------------------
	float m_timerValue;

	MPM *m_MPM;

protected:

	//----------------------------------------------------------------------------------------------------------------------
  /// @brief  The following methods must be implimented in the sub class
  /// this is called when the window is created
	//----------------------------------------------------------------------------------------------------------------------
  void initializeGL();

	//----------------------------------------------------------------------------------------------------------------------
  /// @brief this is called whenever the window is re-sized
  /// @param[in] _w the width of the resized window
  /// @param[in] _h the height of the resized window
	//----------------------------------------------------------------------------------------------------------------------
	void resizeGL(const int _w, const int _h);
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief this is the main gl drawing routine which is called whenever the window needs to
  /// be re-drawn
	//----------------------------------------------------------------------------------------------------------------------
  void paintGL();

private :
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief this method is called every time a mouse is moved
  /// @param _event the Qt Event structure
	//----------------------------------------------------------------------------------------------------------------------
	//void mouseMoveEvent ( QMouseEvent * _event );
	//----------------------------------------------------------------------------------------------------------------------
  /// @brief this method is called everytime the mouse button is pressed
  /// inherited from QObject and overridden here.
  /// @param _event the Qt Event structure
	//----------------------------------------------------------------------------------------------------------------------

	//void mousePressEvent ( QMouseEvent *_event );

	//----------------------------------------------------------------------------------------------------------------------
  /// @brief this method is called everytime the mouse button is released
  /// inherited from QObject and overridden here.
  /// @param _event the Qt Event structure
	//----------------------------------------------------------------------------------------------------------------------
	//void mouseReleaseEvent (QMouseEvent *_event );
	//----------------------------------------------------------------------------------------------------------------------
	/// @brief timer event trigered by startTimer
	//----------------------------------------------------------------------------------------------------------------------
	void timerEvent(QTimerEvent *_event	);
	void loadMatricesToShader();
	void loadMatricesToColourShader();

};

#endif
