/***************************************************************************
 * evalc.cpp is part of Math Graphic Library
 * Copyright (C) 2007-2016 Alexey Balakin <mathgl.abalakin@gmail.ru>       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License  as       *
 *   published by the Free Software Foundation; either version 3 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public     *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <time.h>
#include "mgl2/datac.h"
#include "mgl2/evalc.h"
#if MGL_HAVE_GSL
#include <gsl/gsl_sf.h>
#endif
//-----------------------------------------------------------------------------
//	constants for expression parsing
enum{
EQ_NUM=0,	// a variable substitution
EQ_RND,		// random number
EQ_A,		// numeric constant
EQ_S,		// subdata for given variable
// normal functions of 2 arguments
EQ_LT,		// comparison x<y			!!! MUST BE FIRST 2-PLACE FUNCTION
EQ_GT,		// comparison x>y
EQ_EQ,		// comparison x=y
EQ_ADD,		// addition x+y
EQ_SUB,		// substraction x-y
EQ_MUL,		// multiplication x*y
EQ_DIV,		// division x/y
EQ_IPOW,	// power x^n for integer n
EQ_POW,		// power x^y
EQ_LOG,		// logarithm of x on base a, log_a(x) = ln(x)/ln(a)
EQ_CMPLX,	// return a+i*b
EQ_HYPOT,	// return sqrt(a*a+b*b)
// normal functions of 1 argument
EQ_SIN,		// sine function \sin(x).			!!! MUST BE FIRST 1-PLACE FUNCTION
EQ_COS,		// cosine function \cos(x).
EQ_TAN,		// tangent function \tan(x).
EQ_ASIN,	// inverse sine function \asin(x).
EQ_ACOS,	// inverse cosine function \acos(x).
EQ_ATAN,	// inverse tangent function \atan(x).
EQ_SINH,	// hyperbolic sine function \sinh(x).
EQ_COSH,	// hyperbolic cosine function \cosh(x).
EQ_TANH,	// hyperbolic tangent function \tanh(x).
EQ_ASINH,	// inverse hyperbolic sine function \asinh(x).
EQ_ACOSH,	// inverse hyperbolic cosine function \acosh(x).
EQ_ATANH,	// inverse hyperbolic tangent function \atanh(x).
EQ_SQRT,	// square root function \sqrt(x)
EQ_EXP,		// exponential function \exp(x)
EQ_EXPI,	// exponential function \exp(i*x)
EQ_LN,		// logarithm of x, ln(x)
EQ_LG,		// decimal logarithm of x, lg(x) = ln(x)/ln(10)
EQ_ABS,		// absolute value
EQ_ARG,		// argument (or phase) of complex number
EQ_CONJ,	// complex conjugate
EQ_REAL,	// real part
EQ_IMAG,	// imaginary part
EQ_NORM,	// square of absolute value |u|^2
EQ_LAST		// id of last entry
};
//-----------------------------------------------------------------------------
int mglFormulaC::Error=0;
bool MGL_LOCAL_PURE mglCheck(char *str,int n);
long MGL_LOCAL_PURE mglFindInText(const char *str, const char *lst);
//-----------------------------------------------------------------------------
mglFormulaC::~mglFormulaC()
{
	if(tmp)		delete tmp;
	if(Left)	delete Left;
	if(Right)	delete Right;
}
//-----------------------------------------------------------------------------
// Formula constructor (automatically parse and "compile" formula)
mglFormulaC::mglFormulaC(const char *string):Res(0)
{
	dat = tmp = NULL;
	dx1=dy1=dz1=0;	dx2=dy2=dz2=1;
	Error=Kod=0;
	Left=Right=0;
	if(!string)	{	Kod = EQ_NUM;	Res = 0;	return;	}
	char *str = new char[strlen(string)+1];
	strcpy(str,string);
	long n,len;
	mgl_strtrim(str);
	mgl_strlwr(str);
	len=strlen(str);
	if(str[0]==0) {	delete []str;	return;	}
	if(str[0]=='(' && mglCheck(&(str[1]),len-2))	// remove braces
	{
		memmove(str,str+1,len);
		len-=2;	str[len]=0;
	}
	len=strlen(str);
	if(str[0]==':' && str[1]!=0)		//	this data file for interpolation
	{
		double sx1,sx2,sy1,sy2,sz1,sz2;
		char *buf = strchr(str+1,':');
		if(buf && *buf)
		{
			*buf = 0;
			int r = sscanf(buf+1,"%lg:%lg:%lg:%lg:%lg:%lg",&sx1,&sx2,&sy1,&sy2,&sz1,&sz2);
			if(r>1 && sx1!=sx2)	{	dx1=sx1;	dx2=sx2;	}
			if(r>3 && sy1!=sy2)	{	dy1=sy1;	dy2=sy2;	}
			if(r>5 && sz1!=sz2)	{	dz1=sz1;	dz2=sz2;	}
		}
		dat = tmp = new mglDataC(str+1);
		delete []str;	return;
	}
	n=mglFindInText(str,"<>=");				// low priority -- conditions
	if(n>=0)
	{
		if(str[n]=='<') Kod=EQ_LT;
		else if(str[n]=='>') Kod=EQ_GT;
		else Kod=EQ_EQ;
		str[n]=0;
		Left=new mglFormulaC(str);
		Right=new mglFormulaC(str+n+1);
		delete []str;	return;
	}
	n=mglFindInText(str,"+-");				// normal priority -- additions
	if(n>=0 && (n<2 || !strchr("eE",str[n-1]) || (str[n-2]!='.' && !isdigit(str[n-2]))))
	{
		if(str[n]=='+') Kod=EQ_ADD; else Kod=EQ_SUB;
		str[n]=0;
		Left=new mglFormulaC(str);
		Right=new mglFormulaC(str+n+1);
		delete []str;	return;
	}
	n=mglFindInText(str,"*/");				// high priority -- multiplications
	if(n>=0)
	{
		if(str[n]=='*') Kod=EQ_MUL; else Kod=EQ_DIV;
		str[n]=0;
		Left=new mglFormulaC(str);
		Right=new mglFormulaC(str+n+1);
		delete []str;	return;
	}
	n=mglFindInText(str,"^");				// highest priority -- power
	if(n>=0)
	{
		Kod=EQ_IPOW;		str[n]=0;
		Left=new mglFormulaC(str);
		Right=new mglFormulaC(str+n+1);
		delete []str;	return;
	}

	for(n=0;n<len;n++)	if(str[n]=='(')	break;
	if(n>=len)								// this is number or variable
	{
		Kod = EQ_NUM;
//		Left = Right = 0;
		if(str[1]==0 && str[0]>='a' && str[0]<='z')	// available variables
		{	Kod=EQ_A;	Res = str[0]-'a';	}
		else if(!strcmp(str,"rnd"))	Kod=EQ_RND;
		else if(!strcmp(str,"pi"))	Res=M_PI;
		else if(!strcmp(str,"inf"))	Res=INFINITY;
		else if(str[0]=='i')	Res = dual(0,atof(str+1));
		else	Res = (str[len-1]=='i') ? dual(0,atof(str)) : atof(str);
	}
	else
	{
		char name[128];
		mgl_strncpy(name,str,128);	name[127]=name[n]=0;
		memmove(str,str+n+1,len-n);
		len=strlen(str);		str[--len]=0;
		if(strlen(name)==1)
		{	Kod=EQ_S;	Res = name[0]-'a';	}
		else
		{
			if(!strcmp(name,"sin")) Kod=EQ_SIN;
			else if(!strcmp(name,"cos")) Kod=EQ_COS;
			else if(!strcmp(name,"tg")) Kod=EQ_TAN;
			else if(!strcmp(name,"tan")) Kod=EQ_TAN;
			else if(!strcmp(name,"asin")) Kod=EQ_ASIN;
			else if(!strcmp(name,"acos")) Kod=EQ_ACOS;
			else if(!strcmp(name,"atan")) Kod=EQ_ATAN;
			else if(!strcmp(name,"sinh")) Kod=EQ_SINH;
			else if(!strcmp(name,"cosh")) Kod=EQ_COSH;
			else if(!strcmp(name,"tanh")) Kod=EQ_TANH;
			else if(!strcmp(name,"sh")) Kod=EQ_SINH;
			else if(!strcmp(name,"ch")) Kod=EQ_COSH;
			else if(!strcmp(name,"th")) Kod=EQ_TANH;
			else if(!strcmp(name,"sqrt")) Kod=EQ_SQRT;
			else if(!strcmp(name,"log")) Kod=EQ_LOG;
			else if(!strcmp(name,"pow")) Kod=EQ_POW;
			else if(!strcmp(name,"exp")) Kod=EQ_EXP;
			else if(!strcmp(name,"lg")) Kod=EQ_LG;
			else if(!strcmp(name,"ln")) Kod=EQ_LN;
			else if(!strcmp(name,"abs")) Kod=EQ_ABS;
			else if(!strcmp(name,"arg")) Kod=EQ_ARG;
			else if(!strcmp(name,"conj")) Kod=EQ_CONJ;
			else if(!strcmp(name,"real")) Kod=EQ_REAL;
			else if(!strcmp(name,"imag")) Kod=EQ_IMAG;
			else if(!strcmp(name,"norm")) Kod=EQ_NORM;
			else if(!strcmp(name,"cmplx")) Kod=EQ_CMPLX;
			else if(!strcmp(name,"hypot")) Kod=EQ_HYPOT;
			else {	delete []str;	return;	}	// unknown function
		}
		n=mglFindInText(str,",");
		if(n>=0)
		{
			str[n]=0;
			Left=new mglFormulaC(str);
			Right=new mglFormulaC(str+n+1);
		}
		else
			Left=new mglFormulaC(str);
	}
	delete []str;
}
//-----------------------------------------------------------------------------
// evaluate formula for 'x'='r', 'y'='n'='v', 't'='z', 'u'='a' variables
dual mglFormulaC::Calc(dual x,dual y,dual t,dual u) const
{
	Error=0;
	dual a1[MGL_VS];	memset((void*)a1,0,MGL_VS*sizeof(dual));
	a1['a'-'a'] = a1['u'-'a'] = u;
	a1['x'-'a'] = a1['r'-'a'] = x;
	a1['y'-'a'] = a1['n'-'a'] = a1['v'-'a'] = y;
	a1['z'-'a'] = a1['t'-'a'] = t;
	a1['i'-'a'] = dual(0,1);
	dual b = CalcIn(a1);
	return mgl_isfin(b) ? b : NAN;
}
//-----------------------------------------------------------------------------
// evaluate formula for 'x'='r', 'y'='n', 't'='z', 'u'='a', 'v'='b', 'w'='c' variables
dual mglFormulaC::Calc(dual x,dual y,dual t,dual u,dual v,dual w) const
{
	Error=0;
	dual a1[MGL_VS];	memset((void*)a1,0,MGL_VS*sizeof(dual));
	a1['c'-'a'] = a1['w'-'a'] = w;
	a1['b'-'a'] = a1['v'-'a'] = v;
	a1['a'-'a'] = a1['u'-'a'] = u;
	a1['x'-'a'] = a1['r'-'a'] = x;
	a1['y'-'a'] = a1['n'-'a'] = y;
	a1['z'-'a'] = a1['t'-'a'] = t;
	a1['i'-'a'] = dual(0,1);
	dual b = CalcIn(a1);
	return mgl_isfin(b) ? b : NAN;
}
//-----------------------------------------------------------------------------
// evaluate formula for arbitrary set of variables
dual mglFormulaC::Calc(const dual var[MGL_VS]) const
{
	Error=0;
	dual b = CalcIn(var);
	return mgl_isfin(b) ? b : NAN;
}
//-----------------------------------------------------------------------------
dual MGL_LOCAL_CONST ceqc(const dual &a,const dual &b)	{return a==b?1:0;}
dual MGL_LOCAL_CONST cltc(const dual &a,const dual &b)	{return real(a-b)<0?1:0;}
dual MGL_LOCAL_CONST cgtc(const dual &a,const dual &b)	{return real(a-b)>0?1:0;}
dual MGL_LOCAL_CONST addc(const dual &a,const dual &b)	{return a+b;}
dual MGL_LOCAL_CONST subc(const dual &a,const dual &b)	{return a-b;}
dual MGL_LOCAL_CONST mulc(const dual &a,const dual &b)	{return a*b;}
dual MGL_LOCAL_CONST divc(const dual &a,const dual &b)	{return a/b;}
dual MGL_LOCAL_CONST ipwc(const dual &a,const dual &b)	{return mgl_ipowc(a,int(b.real()));}
dual MGL_LOCAL_CONST powc(const dual &a,const dual &b)	{return exp(b*log(a));	}
dual MGL_LOCAL_CONST llgc(const dual &a,const dual &b)	{return log(a)/log(b);	}
dual MGL_LOCAL_CONST cmplxc(const dual &a,const dual &b)	{return a+dual(0,1)*b;	}
dual MGL_LOCAL_CONST expi(const dual &a)	{	return exp(dual(0,1)*a);	}
dual MGL_LOCAL_CONST expi(const double &a)	{	return dual(cos(a),sin(a));	}
//-----------------------------------------------------------------------------
const dual MGL_NO_EXPORT ic = dual(0,1);
dual MGL_LOCAL_CONST hypotc(const dual &x, const dual &y)	{	return sqrt(x*x+y*y);	}
dual MGL_LOCAL_CONST asinhc(const dual &x)	{	return log(x+sqrt(x*x+mreal(1)));	}
dual MGL_LOCAL_CONST acoshc(const dual &x)	{	return log(x+sqrt(x*x-mreal(1)));	}
dual MGL_LOCAL_CONST atanhc(const dual &x)	{	return log((mreal(1)+x)/(mreal(1)-x))/mreal(2);	}
dual MGL_LOCAL_CONST conjc(const dual &x)	{	return dual(real(x),-imag(x));	}
dual MGL_LOCAL_CONST sinc(const dual &x)	{	return sin(x);	}
dual MGL_LOCAL_CONST cosc(const dual &x)	{	return cos(x);	}
dual MGL_LOCAL_CONST tanc(const dual &x)	{	return tan(x);	}
dual MGL_LOCAL_CONST sinhc(const dual &x)	{	return sinh(x);	}
dual MGL_LOCAL_CONST coshc(const dual &x)	{	return cosh(x);	}
dual MGL_LOCAL_CONST tanhc(const dual &x)	{	return tanh(x);	}
dual MGL_LOCAL_CONST asinc(const dual &x)	{	return log(ic*x+sqrt(mreal(1)-x*x))/ic;	}
dual MGL_LOCAL_CONST acosc(const dual &x)	{	return log(x+sqrt(x*x-mreal(1)))/ic;	}
dual MGL_LOCAL_CONST atanc(const dual &x)	{	return log((ic-x)/(ic+x))/(mreal(2)*ic);	}
dual MGL_LOCAL_CONST expc(const dual &x)	{	return exp(x);	}
dual MGL_LOCAL_CONST sqrtc(const dual &x)	{	return sqrt(x);	}
dual MGL_LOCAL_CONST logc(const dual &x)	{	return log(x);	}
dual MGL_LOCAL_CONST absc(const dual &x)	{	return abs(x);	}
dual MGL_LOCAL_CONST argc(const dual &x)	{	return arg(x);	}
dual MGL_LOCAL_CONST lgc(const dual &x)	{	return log10(x);}
dual MGL_LOCAL_CONST realc(const dual &x)	{	return real(x);	}
dual MGL_LOCAL_CONST imagc(const dual &x)	{	return imag(x);	}
dual MGL_LOCAL_CONST normc(const dual &x)	{	return norm(x);	}
//-----------------------------------------------------------------------------
typedef dual (*func_1)(const dual&);
typedef dual (*func_2)(const dual&, const dual&);
static const func_2 f2[EQ_SIN-EQ_LT] = {cltc,cgtc,ceqc,addc,subc,mulc,divc,ipwc,powc,llgc,cmplxc,hypotc};
static const func_1 f1[EQ_LAST-EQ_SIN] = {sinc,cosc,tanc,asinc,acosc,atanc,sinhc,coshc,tanhc,
					asinhc,acoshc,atanhc,sqrtc,expc,expi,logc,lgc,absc,argc,conjc,realc,imagc,normc};
// evaluation of embedded (included) expressions
dual mglFormulaC::CalcIn(const dual *a1) const
{
	if(dat)
	{
		mreal x = (real(a1['x'-'a'])-dx1)*(dat->GetNx()-1)/(dx2-dx1);
		mreal y = (real(a1['y'-'a'])-dy1)*(dat->GetNy()-1)/(dy2-dy1);
		mreal z = (real(a1['z'-'a'])-dz1)*(dat->GetNz()-1)/(dz2-dz1);
		return mgl_datac_spline(dat,x,y,z);
	}
	if(Kod<EQ_LT)
	{
		if(Kod==EQ_RND)	return mgl_rnd();
		else	return (Kod==EQ_A) ? a1[int(Res.real())] : Res;
	}

	dual a = Left->CalcIn(a1);
	if(mgl_isfin(a))
	{
		if(Kod<EQ_SIN)
			return Right?f2[Kod-EQ_LT](a,Right->CalcIn(a1)):NAN;
		else
			return f1[Kod-EQ_SIN](a);
	}
	return NAN;
}
//-----------------------------------------------------------------------------
void mglFormulaC::CalcV(HADT res, HCDT var[MGL_VS]) const
{
	const long nn = res->GetNN();
	if(dat)
	{
		HCDT X = var['x'-'a'], Y = var['y'-'a'], Z = var['z'-'a'];
		long nx = X?X->GetNN():0, ny = Y?Y->GetNN():0, nz = Z?Z->GetNN():0;
		if(nx==nn && ny==nn && nz==nn)
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,X->vthr(i),Y->vthr(i),Z->vthr(i));
		else if(nx==nn && ny==nn)
		{
			double a = Z?Z->v(0):0;
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,X->vthr(i),Y->vthr(i),a);
		}
		else if(nx==nn && nz==nn)
		{
			double a = Y?Y->v(0):0;
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,X->vthr(i),a,Z->vthr(i));
		}
		else if(nz==nn && ny==nn)
		{
			double a = X?X->v(0):0;
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,a,Y->vthr(i),Z->vthr(i));
		}
		else if(nx==nn)
		{
			double a = Y?Y->v(0):0, b = Z?Z->v(0):0;
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,X->vthr(i),a,b);
		}
		else if(ny==nn)
		{
			double a = X?X->v(0):0, b = Z?Z->v(0):0;
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,a,Y->vthr(i),b);
		}
		else if(nz==nn)
		{
			double a = X?X->v(0):0, b = Y?Y->v(0):0;
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,a,b,Z->vthr(i));
		}
		else
		{
			double a = X?X->v(0):0, b = Y?Y->v(0):0, c = Z?Z->v(0):0;
			double v = mgl_data_spline(dat,a,b,c);
			for(long i=0;i<nn;i++)	res->a[i] = v;
		}
		return;
	}
	if(Kod<EQ_LT)
	{
		int k = int(Res.real()+0.5);
		if(Kod==EQ_NUM)
			for(long i=0;i<nn;i++)	res->a[i] = Res;
		else if(Kod==EQ_RND)	
			for(long i=0;i<nn;i++)	res->a[i] = mgl_rnd();
		else if(Kod==EQ_A)
		{
			HCDT a = var[k];
			const mglDataC *c = dynamic_cast<const mglDataC *>(a);
			if(a)
			{
				if(a->GetNN() < nn)
				{
					const dual val = c?c->a[0]:a->v(0);
					for(long i=0;i<nn;i++)	res->a[i] = val;
				}
				else
				{
					if(c)	for(long i=0;i<nn;i++)	res->a[i] = c->a[i];
					else	for(long i=0;i<nn;i++)	res->a[i] = a->vthr(i);
				}
			}
		}
		else if(Kod==EQ_S)
		{
			Left->CalcV(res,var);
			HCDT a = var[k];
			const mglDataC *c = dynamic_cast<const mglDataC *>(a);
			const long na = a->GetNN();
			if(c)
				for(long i=0;i<nn;i++)
				{
					const long kk = long(res->a[i].real()+0.5);
					if(kk>=0 && kk<na)	res->a[i] = c->a[kk];
					else	res->a[i] = 0.;
				}
			else if(a)
				for(long i=0;i<nn;i++)
				{
					const long kk = long(res->a[i].real()+0.5);
					if(kk>=0 && kk<na)	res->a[i] = a->vthr(kk);
					else	res->a[i] = NAN;
				}
		}
		return;
	}
	Left->CalcV(res,var);
	if(Kod<EQ_SIN)
	{
		const func_2 ff=f2[Kod-EQ_LT];
		const dual &v = Right->Res;
		if(Right->Kod==EQ_NUM)	switch(Kod)
		{
		case EQ_ADD:
			for(long i=0;i<nn;i++)	res->a[i] += v;
			break;
		case EQ_SUB:
			for(long i=0;i<nn;i++)	res->a[i] -= v;
			break;
		case EQ_MUL:
			for(long i=0;i<nn;i++)	res->a[i] *= v;
			break;
		case EQ_DIV:
			for(long i=0;i<nn;i++)	res->a[i] /= v;
			break;
		default:
			for(long i=0;i<nn;i++)	res->a[i] = ff(res->a[i], v);
		}
		else if(Right->Kod==EQ_RND)	switch(Kod)
		{
		case EQ_ADD:
			for(long i=0;i<nn;i++)	res->a[i] += mgl_rnd();
			break;
		case EQ_SUB:
			for(long i=0;i<nn;i++)	res->a[i] -= mgl_rnd();
			break;
		case EQ_MUL:
			for(long i=0;i<nn;i++)	res->a[i] *= mgl_rnd();
			break;
		case EQ_DIV:
			for(long i=0;i<nn;i++)	res->a[i] /= mgl_rnd();
			break;
		default:
			for(long i=0;i<nn;i++)	res->a[i] = ff(res->a[i], mgl_rnd());
		}
		else
		{
			mglDataC tmp(res->nx,res->ny,res->nz);
			Right->CalcVomp(&tmp, var);
			switch(Kod)
			{
			case EQ_ADD:
				for(long i=0;i<nn;i++)	res->a[i] += tmp.a[i];
				break;
			case EQ_SUB:
				for(long i=0;i<nn;i++)	res->a[i] -= tmp.a[i];
				break;
			case EQ_MUL:
				for(long i=0;i<nn;i++)	res->a[i] *= tmp.a[i];
				break;
			case EQ_DIV:
				for(long i=0;i<nn;i++)	res->a[i] /= tmp.a[i];
				break;
			default:
				for(long i=0;i<nn;i++)	res->a[i] = ff(res->a[i], tmp.a[i]);
			}
		}
	}
	else if(Kod<EQ_LAST)
	{
		const func_1 ff=f1[Kod-EQ_SIN];
		for(long i=0;i<nn;i++)
			res->a[i] = ff(res->a[i]);
	}
	else
		for(long i=0;i<nn;i++)	res->a[i] = NAN;
}
//-----------------------------------------------------------------------------
void mglFormulaC::CalcVomp(HADT res, HCDT var[MGL_VS]) const
{
	const long nn = res->GetNN();
	if(dat)
	{
		HCDT X = var['x'-'a'], Y = var['y'-'a'], Z = var['z'-'a'];
		long nx = X?X->GetNN():0, ny = Y?Y->GetNN():0, nz = Z?Z->GetNN():0;
		if(nx==nn && ny==nn && nz==nn)
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,X->vthr(i),Y->vthr(i),Z->vthr(i));
		else if(nx==nn && ny==nn)
		{
			double a = Z?Z->v(0):0;
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,X->vthr(i),Y->vthr(i),a);
		}
		else if(nx==nn && nz==nn)
		{
			double a = Y?Y->v(0):0;
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,X->vthr(i),a,Z->vthr(i));
		}
		else if(nz==nn && ny==nn)
		{
			double a = X?X->v(0):0;
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,a,Y->vthr(i),Z->vthr(i));
		}
		else if(nx==nn)
		{
			double a = Y?Y->v(0):0, b = Z?Z->v(0):0;
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,X->vthr(i),a,b);
		}
		else if(ny==nn)
		{
			double a = X?X->v(0):0, b = Z?Z->v(0):0;
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,a,Y->vthr(i),b);
		}
		else if(nz==nn)
		{
			double a = X?X->v(0):0, b = Y?Y->v(0):0;
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = mgl_datac_spline(dat,a,b,Z->vthr(i));
		}
		else
		{
			double a = X?X->v(0):0, b = Y?Y->v(0):0, c = Z?Z->v(0):0;
			double v = mgl_data_spline(dat,a,b,c);
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = v;
		}
		return;
	}
	if(Kod<EQ_LT)
	{
		int k = int(Res.real()+0.5);
		if(Kod==EQ_NUM)
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = Res;
		else if(Kod==EQ_RND)	
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = mgl_rnd();
		else if(Kod==EQ_A)
		{
			HCDT a = var[k];
			const mglDataC *c = dynamic_cast<const mglDataC *>(a);
			if(a)
			{
				if(a->GetNN() < nn)
				{
					const dual val = c?c->a[0]:a->v(0);
					for(long i=0;i<nn;i++)	res->a[i] = val;
				}
				else
				{
					if(c)	for(long i=0;i<nn;i++)	res->a[i] = c->a[i];
					else	for(long i=0;i<nn;i++)	res->a[i] = a->vthr(i);
				}
			}
		}
		else if(Kod==EQ_S)
		{
			Left->CalcVomp(res,var);
			HCDT a = var[k];
			const mglDataC *c = dynamic_cast<const mglDataC *>(a);
			const long na = a?a->GetNN():0;
			if(c)
#pragma omp parallel for schedule(static)
				for(long i=0;i<nn;i++)
				{
					const long kk = long(res->a[i].real()+0.5);
					if(kk>=0 && kk<na)	res->a[i] = c->a[kk];
					else	res->a[i] = 0.;
				}
			else if(a)
#pragma omp parallel for schedule(static)
				for(long i=0;i<nn;i++)
				{
					const long kk = long(res->a[i].real()+0.5);
					if(kk>=0 && kk<na)	res->a[i] = a->vthr(kk);
					else	res->a[i] = NAN;
				}
		}
		return;
	}
	Left->CalcVomp(res,var);
	if(Kod<EQ_SIN)
	{
		const func_2 ff=f2[Kod-EQ_LT];
		const dual &v = Right->Res;
		if(Right->Kod==EQ_NUM)	switch(Kod)
		{
		case EQ_ADD:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] += v;
			break;
		case EQ_SUB:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] -= v;
			break;
		case EQ_MUL:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] *= v;
			break;
		case EQ_DIV:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] /= v;
			break;
		default:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = ff(res->a[i], v);
		}
		else if(Right->Kod==EQ_RND)	switch(Kod)
		{
		case EQ_ADD:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] += mgl_rnd();
			break;
		case EQ_SUB:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] -= mgl_rnd();
			break;
		case EQ_MUL:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] *= mgl_rnd();
			break;
		case EQ_DIV:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] /= mgl_rnd();
			break;
		default:
#pragma omp parallel for schedule(static)
			for(long i=0;i<nn;i++)	res->a[i] = ff(res->a[i], mgl_rnd());
		}
		else
		{
			mglDataC tmp(res->nx,res->ny,res->nz);
			Right->CalcVomp(&tmp, var);
			switch(Kod)
			{
			case EQ_ADD:
#pragma omp parallel for schedule(static)
				for(long i=0;i<nn;i++)	res->a[i] += tmp.a[i];
				break;
			case EQ_SUB:
#pragma omp parallel for schedule(static)
				for(long i=0;i<nn;i++)	res->a[i] -= tmp.a[i];
				break;
			case EQ_MUL:
#pragma omp parallel for schedule(static)
				for(long i=0;i<nn;i++)	res->a[i] *= tmp.a[i];
				break;
			case EQ_DIV:
#pragma omp parallel for schedule(static)
				for(long i=0;i<nn;i++)	res->a[i] /= tmp.a[i];
				break;
			default:
#pragma omp parallel for schedule(static)
				for(long i=0;i<nn;i++)	res->a[i] = ff(res->a[i], tmp.a[i]);
			}
		}
	}
	else if(Kod<EQ_LAST)
	{
		const func_1 ff=f1[Kod-EQ_SIN];
#pragma omp parallel for schedule(static)
		for(long i=0;i<nn;i++)
			res->a[i] = ff(res->a[i]);
	}
	else
#pragma omp parallel for schedule(static)
		for(long i=0;i<nn;i++)	res->a[i] = NAN;
}
//-----------------------------------------------------------------------------
dual MGL_LOCAL_CONST mgl_ipowc_c(dual x,int n)
{
	dual t;
	if(n==2)	t = x*x;
	else if(n==1)	t = x;
	else if(n<0)	t = mreal(1)/mgl_ipowc_c(x,-n);
	else if(n==0)	t = mreal(1);
	else
	{
		t = mgl_ipowc_c(x,n/2);	t = t*t;
		if(n%2==1)	t *= x;
	}
	return t;
}
cmdual MGL_EXPORT_CONST mgl_ipowc(mdual x,int n)
{	return mdual(mgl_ipowc_c(x,n));	}
cmdual MGL_EXPORT_PURE mgl_ipowc_(mdual *x,int *n)	{	return mgl_ipowc(*x,*n);	}
//-----------------------------------------------------------------------------
HAEX MGL_EXPORT mgl_create_cexpr(const char *expr)	{	return new mglFormulaC(expr);	}
uintptr_t MGL_EXPORT mgl_create_cexpr_(const char *expr, int l)
{	char *s=new char[l+1];	memcpy(s,expr,l);	s[l]=0;
	uintptr_t res = uintptr_t(mgl_create_cexpr(s));
	delete []s;	return res;	}
void MGL_EXPORT mgl_delete_cexpr(HAEX ex)	{	if(ex)	delete ex;	}
void MGL_EXPORT mgl_delete_cexpr_(uintptr_t *ex)	{	mgl_delete_cexpr((HAEX)ex);	}
cmdual MGL_EXPORT mgl_cexpr_eval(HAEX ex, mdual x, mdual y, mdual z)
{	return mdual(ex->Calc(x,y,z));	}
cmdual MGL_EXPORT mgl_cexpr_eval_(uintptr_t *ex, mdual *x, mdual *y, mdual *z)
{	return mgl_cexpr_eval((HAEX) ex, *x,*y,*z);		}
cmdual MGL_EXPORT mgl_cexpr_eval_v(HAEX ex, mdual *var)
{	return mdual(ex->Calc(reinterpret_cast<dual*>(var)));	}
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_cexpr_eval_dat(HAEX ex, HADT res, HCDT vars[MGL_VS])
{	ex->CalcV(res,vars);	}
void MGL_EXPORT mgl_cexpr_eval_dat_(uintptr_t *ex, uintptr_t *res, uintptr_t vars[MGL_VS])
{
	HCDT vs[MGL_VS];
	for(int i=0;i<MGL_VS;i++)	vs[i] = (HCDT)(vars[i]);
	mgl_cexpr_eval_dat((HAEX) *ex, (HADT) *res, vs);
}
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_cexpr_eval_omp(HAEX ex, HADT res, HCDT vars[MGL_VS])
{	ex->CalcVomp(res,vars);	}
void MGL_EXPORT mgl_cexpr_eval_omp_(uintptr_t *ex, uintptr_t *res, uintptr_t vars[MGL_VS])
{
	HCDT vs[MGL_VS];
	for(int i=0;i<MGL_VS;i++)	vs[i] = (HCDT)(vars[i]);
	mgl_cexpr_eval_omp((HAEX) *ex, (HADT) *res, vs);
}
//-----------------------------------------------------------------------------
