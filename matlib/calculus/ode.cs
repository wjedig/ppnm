using System;
using static System.Console;
using System.IO;
using static System.Math;
using linalg;
using genlist;

namespace calculus{
public static class rkint{
	public static (vector,vector) rkstep12(
			Func<double,vector,vector> f, /* the f from dy/dx=f(x,y) */
			double x,                    /* the current value of the variable */
			vector y,                    /* the current value y(x) of the sought function */
			double h                     /* the step to be taken */
			)
	{
		vector k0 = f(x,y);			/* embedded lower order formula (Euler) */
		vector k1 = f(x+h/2, y+k0*(h/2));	/* higher order formula (midpoint) */
		vector yh = y+k1*h;			/* y(x+h) estimate */
		vector er = (k1-k0)*h;			/* error estimate */
		return (yh,er);
	}
	public static vector driver(
			Func<double,vector,vector> f,	/* the f from dy/dx=f(x,y) */
			double a,			/* the start-point a */
			vector ya,			/* y(a) */
			double b,			/* the end-point of the integration */
			//ref double count,
			double h=0.01,			/* initial step-size */
			double acc=0.01,		/* absolute accuracy goal */
			double eps=0.01,		/* relative accuracy goal */
			genlist<double> xlist = null, 
			genlist<vector> ylist = null,
			double maxstep = 0.0;
			)
	{
		if(a>b) throw new ArgumentException("driver: a>b");
		//count ++;
		double x=a; vector y=ya.copy();
		bool isxlist = (xlist!=null);
		bool isylist = (ylist!=null);
		if(isxlist) xlist.add(x);
		if(isylist) ylist.add(y);
		//WriteLine($"xlist exists: {isxlist}, ylist exists: {isylist}");
		do{
			if(x>=b) return y; /* job done! */
			if(x+h>b) h=b-x; /* last step should end at b */
			var (yh,erv) = rkstep12(f,x,y,h);
			vector tol = new vector(y.size); 
			bool accept = true;
			for(int i=0; i<y.size; i++){
				tol[i] = (acc+eps*Abs(yh[i]))*Sqrt(h/(b-a));
				if(!(erv[i]<tol[i])) accept=false;
			}
			if(accept){
				x+=h; y=yh; 
				if(isxlist) xlist.add(x);
				if(isylist) ylist.add(y);
			}
			double factor = tol[0]/Abs(erv[0]);
			for(int i=1 ; i<y.size ; i++) factor=Min(factor,tol[i]/Abs(erv[i]));
			h *= Min( Pow(factor,0.25)*0.95 ,2);
			if(maxstep>0 && h>maxstep) h = maxstep;
		} while(true);
	} //driver
	//public static vector driver_count(
	//		Func<double,vector,vector> f,	/* the f from dy/dx=f(x,y) */
	//		double a,			/* the start-point a */
	//		vector ya,			/* y(a) */
	//		double b,			/* the end-point of the integration */
	//		double h=0.01,			/* initial step-size */
	//		double acc=0.01,		/* absolute accuracy goal */
	//		double eps=0.01,		/* relative accuracy goal */
	//		genlist<double> xlist = null, genlist<vector> ylist = null
	//		)
	//{
	//	double c = 0;
	//	return driver(f,a,ya,b,h,acc,eps,xlist,ylist);
	//}
} //rkint
} //calculus
