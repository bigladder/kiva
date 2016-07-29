/***************************************************************************
 * Fl_MathGL.h is part of Math Graphic Library
 * Copyright (C) 2007-2016 Alexey Balakin <mathgl.abalakin@gmail.ru>       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library General Public License as       *
 *   published by the Free Software Foundation; either version 3 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef _MGL_FL_MATHGL_H_
#define _MGL_FL_MATHGL_H_

#ifdef __MWERKS__
# define FL_DLL
#endif

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Counter.H>
#include <FL/Fl_Menu_Bar.H>
#include <mgl2/fltk.h>
class mglCanvas;
//-----------------------------------------------------------------------------
/// Class is FLTK widget which display MathGL graphics
class MGL_EXPORT Fl_MathGL : public Fl_Widget
{
	friend class Fl_MGLView;
public:
	Fl_Valuator	*tet_val;	///< pointer to external tet-angle validator
	Fl_Valuator	*phi_val;	///< pointer to external phi-angle validator
	mglCanvas *gr;			///< Built-in mglCanvas instance (mglCanvasFLTK is used by default)

	Fl_MathGL(int x, int y, int w, int h, const char *label=0);
	virtual ~Fl_MathGL();

	/// Update (redraw) plot
	virtual void update();
	/// Set angles for additional plot rotation
	inline void set_angle(mreal t, mreal p){	tet = t;	phi = p;	}
	/// Set bitwise flags for general state (1-Alpha, 2-Light)
	inline void set_flag(int f)	{	flag = f;	}
	/// Set flags for handling mouse
	void set_graph(HMGL gr);	///< Set grapher object
	inline void set_graph(mglGraph *Gr)
	{	set_graph(Gr->Self());	}
	/// Get pointer to grapher
	inline HMGL get_graph()	{	return (HMGL)gr;	}
	/// Get mglDraw pointer or NULL
	inline mglDraw *get_class()
	{	mglDraw *d=0;
		if(draw_func==mgl_draw_class)	d = (mglDraw*)draw_par;
		if(draw_cl)	d = draw_cl;	return d;	}
	/// Set drawing functions and its parameter
	inline void set_draw(int (*func)(mglBase *gr, void *par), void *par)
	{	if(draw_cl)	delete draw_cl;	draw_cl=0;	draw_func=func;	draw_par=par;	}
	inline void set_draw(mglDraw *dr)	{	if(draw_cl)	delete draw_cl;	draw_cl=dr;	draw_func=0;	}
	inline void set_draw(int (*dr)(mglGraph *gr))
	{	set_draw(dr?mgl_draw_graph:0,(void*)dr);	}
	void set_state(bool z, bool r)	{	zoom = z;	rotate = r;	}
	/// Set zoom in/out region
	inline void set_zoom(mreal X1, mreal Y1, mreal X2, mreal Y2)
	{	x1 = X1;	x2 = X2;	y1 = Y1;	y2 = Y2;	update();	}
	/// Get zoom region
	inline void get_zoom(mreal *X1, mreal *Y1, mreal *X2, mreal *Y2)
	{	*X1 = x1;	*X2 = x2;	*Y1 = y1;	*Y2 = y2;	}
	/// Set popup menu pointer
	inline void set_popup(const Fl_Menu_Item *pmenu, Fl_Widget *wdg, void *v)
	{	popup = pmenu;	wpar = wdg;	vpar = v;	}
	inline void zoom_region(mreal xx1,mreal xx2,mreal yy1, mreal yy2)
	{	x1=xx1;	y1=yy1;	x2=xx2;	y2=yy2;	}
	void stop(bool stop=true);

protected:
	void *draw_par;		///< Parameters for drawing function mglCanvasWnd::DrawFunc.
	/// Drawing function for window procedure. It should return the number of frames.
	int (*draw_func)(mglBase *gr, void *par);
	mglDraw *draw_cl;

	const Fl_Menu_Item *popup;	///< pointer to popup menu items
	Fl_Widget *wpar;			///< widget for popup menu
	void *vpar;					///< parameter for popup menu
	mreal tet,phi;				///< rotation angles
	bool rotate;				///< flag for handle mouse
	bool zoom;					///< flag for zoom by mouse
	bool wire;
	mreal x1,x2,y1,y2;			///< zoom region
	int flag;					///< bitwise flag for general state (1-Alpha, 2-Light)
	int x0,y0,xe,ye;			///< mouse position
	char pos[128];

	virtual void draw();		///< quick drawing function
	int handle(int code);		///< handle mouse events
	void resize(int x, int y, int w, int h);	///< resize control
};
//-----------------------------------------------------------------------------
class MGL_EXPORT Fl_MGLView : public Fl_Window
{
public:
	Fl_MathGL *FMGL;		///< Control which draw graphics
	Fl_Scroll *scroll;
	Fl_Menu_Bar	*menu;

	void *par;				///< Parameter for handling animation
	void (*next)(void*);	///< Callback function for next frame
	void (*prev)(void*);	///< Callback function for prev frame
	mreal (*delay)(void*);	///< Callback function for delay
	void (*reload)(void*);	///< Callback function for reloading

	void toggle_alpha()	{	toggle(alpha, alpha_bt, "Graphics/Alpha");	}
	void toggle_light()	{	toggle(light, light_bt, "Graphics/Light");	}
	void toggle_sshow()	{	toggle(sshow, anim_bt, "Graphics/Slideshow");	}
	void toggle_grid()	{	toggle(grid, grid_bt, "Graphics/Grid");	}
	void toggle_zoom()	{	toggle(zoom, zoom_bt);	}
	void toggle_rotate(){	toggle(rotate, rotate_bt);	}
	void setoff_zoom()	{	setoff(zoom, zoom_bt);	}
	void setoff_rotate(){	setoff(rotate, rotate_bt);	}
	bool is_sshow()		{	return sshow;	}
	void toggle_pause()	{	toggle(pauseC, pause_bt, "Graphics/Pause calc");	exec_pause();	}
	void adjust()
	{	mgl_set_size(FMGL->get_graph(),scroll->w(),scroll->h());	FMGL->size(scroll->w(),scroll->h());	update();	}

	Fl_MGLView(int x, int y, int w, int h, const char *label=0);
	virtual ~Fl_MGLView();
	void update();			///< Update picture by calling user drawing function
protected:
	Fl_Button *alpha_bt, *light_bt, *rotate_bt, *anim_bt, *zoom_bt, *grid_bt, *pause_bt;
//	Fl_Counter *tet, *phi;

	int grid, alpha, light;	///< Current states of wire, alpha, light switches (toggle buttons)
	int sshow, rotate, zoom;///< Current states of slideshow, rotate, zoom switches (toggle buttons)
	int pauseC;	///< Current state of pause for calculations

	void toggle(int &val, Fl_Button *b, const char *txt=NULL);
	void setoff(int &val, Fl_Button *b, const char *txt=NULL);
	void exec_pause();
};
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_makemenu_fltk(Fl_Menu_ *m, Fl_MGLView *w);
//-----------------------------------------------------------------------------
#endif
