using System;
using static System.Console;
using static System.Math;
using linalg;

static class rkint{
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
	public static (genlist<double>, genlist<vector>) driver(
			Func<double,vector,vector> f,	/* the f from dy/dx=f(x,y) */
			double a,			/* the start-point a */
			vector ya,			/* y(a) */
			double b,			/* the end-point of the integration */
			double h=0.01,			/* initial step-size */
			double acc=0.01,		/* absolute accuracy goal */
			double eps=0.01			/* relative accuracy goal */
			)
	{
		if(a>b) throw new ArgumentException("driver: a>b");
		double x=a; vector y=ya.copy();
		var xlist=new genlist<double>(); xlist.add(x);
		var ylist=new genlist<vector>(); ylist.add(y);
		do{
			if(x>=b) return (xlist, ylist); /* job done! */
			if(x+h>b) h=b-x; /* last step should end at b */
			var (yh,erv) = rkstep12(f,x,y,h);
			double tol = (acc+eps*yh.norm()) * Sqrt(h/(b-a));
			double err = erv.norm();
			if(err<=tol){ // accept step
				x+=h; y=yh;
				xlist.add(x);
				ylist.add(y);
			}
			h *= Min(Pow(tol/err,0.25)*0.95 , 2); // reajust stepsize
		} while(true);
	} //driver
} //rkint

class main{
	public static void Main(string[] args){
		foreach(string run in args){
			if(run=="-A2"){
				Func<double, vector, vector> f = (double x, vector y) => { /* Plotting y'' = -y */
					vector yprime = new vector(2);
					yprime[0] = y[1];
					yprime[1] = -y[0];
					return yprime;
				};
				var (xlist, ylist) = rkint.driver(f,0.0,new vector(2.0,0.0),10.0);
				for(int i=0; i<ylist.size; i++){
					vector y = ylist[i];
					double x = xlist[i];
					WriteLine($"{x} {y[0]} {y[1]}");
				}
			}
			if(run=="-A1"){
				Func<double, vector, vector> f = (double x, vector y) => { /* plotting y' = -y */
					vector yprime = new vector(1);
					yprime[0] = -y[0];
					return yprime;
				};
				var (xlist, ylist) = rkint.driver(f,0,new vector(4.0),10.0);
				for(int i=0; i<ylist.size; i++){
					vector y = ylist[i];
					double x = xlist[i];
					WriteLine($"{x} {y[0]}");

				}
			}
			if(run=="-A3"){
				Func<double, vector, vector> f = (double x, vector y) => { /* plotting damp oscil */
					vector yprime = new vector(2);
					yprime[0] = y[1];
					yprime[1] = -0.6*y[1] - 7*y[0];
					return yprime;
				};
				var (xlist, ylist) = rkint.driver(f,0,new vector(4.0,0),10);
				for(int i=0; i<ylist.size; i++){
					vector y = ylist[i];
					double x = xlist[i];
					WriteLine($"{x} {y[0]} {y[1]}");
				}
			}
		}
	}
}
