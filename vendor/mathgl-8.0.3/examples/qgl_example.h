#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <mgl2/mgl.h>
#include <QtGlobal>

#if (QT_VERSION_MAJOR>5)
#include <QOpenGLWidget>
class MainWindow : public QOpenGLWidget
#else
#include <QGLWidget>
class MainWindow : public QGLWidget
#endif
{
	Q_OBJECT
protected:
	mglGraph *gr;			// pointer to MathGL core class
	void resizeGL(int nWidth, int nHeight);	// Method called after each window resize
	void paintGL();			// Method to display the image on the screen
	void initializeGL();	// Method to initialize OpenGL
public:
	MainWindow(QWidget *parent = 0);
	~MainWindow();
};

#endif // MAINWINDOW_H
