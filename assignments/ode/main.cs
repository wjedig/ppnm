using System;
using static System.Console;
using System.IO;
using static System.Math;
using linalg;
using genlist;

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
	public static vector driver(
			Func<double,vector,vector> f,	/* the f from dy/dx=f(x,y) */
			double a,			/* the start-point a */
			vector ya,			/* y(a) */
			double b,			/* the end-point of the integration */
			double h=0.01,			/* initial step-size */
			double acc=0.01,		/* absolute accuracy goal */
			double eps=0.01,		/* relative accuracy goal */
			genlist<double> xlist = null, genlist<vector> ylist = null
			,double maxstep = 0
			)
	{
		if(a>b) throw new ArgumentException("driver: a>b");
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
			if(maxstep>0 && h>maxstep)h=maxstep;
		} while(true);
	} //driver
} //rkint

class main{
	public static void Main(string[] args){
		var xlist = new genlist<double>();
		var ylist = new genlist<vector>();
		foreach(string run in args){
			if(run=="-A2"){
				Func<double, vector, vector> f = (double x, vector y) => { /* Plotting y'' = -y */
					vector yprime = new vector(2);
					yprime[0] = y[1];
					yprime[1] = -y[0];
					return yprime;
				};
				rkint.driver(
					f: f, 
					a: 0.0,
					ya: new vector(2.0,0.0),
					b: 10.0,
					xlist: xlist, 
					ylist: ylist);
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
				rkint.driver(
					f: f, 
					a: 0.0,
					ya: new vector(4.0),
					b: 10.0,
					xlist: xlist, 
					ylist: ylist);
				WriteLine($"xlist size: {xlist.size}");
				for(int i=0; i<ylist.size; i++){
					vector y = ylist[i];
					double x = xlist[i];
					WriteLine($"{x} {y[0]}");
				}
			}
			if(run=="-A3"){
				var outstream = new StreamWriter("Out.A.txt", append:false);
				outstream.WriteLine("--- PART A ---");
				outstream.WriteLine("running for function y'=-y (should be exponential decay)");
				outstream.WriteLine("running for function y''=-y (should be sine/cosine)");
				outstream.WriteLine("Both y'=-y and y''=-y are plotted in ODE_test.svg");
				outstream.WriteLine("running for Dampened Harmonic Oscilator...");
				outstream.WriteLine("This is plotted in ODE_damp.svg");
				outstream.Close();
				Func<double, vector, vector> f = (double x, vector y) => { /* plotting damp oscil */
					vector yprime = new vector(2);
					yprime[0] = y[1];
					yprime[1] = -0.6*y[1] - 7*y[0];
					return yprime;
				};
				rkint.driver(
					f: f, 
					a: 0.0,
					ya: new vector(4.0,0.0),
					b: 10.0,
					xlist: xlist, 
					ylist: ylist);
				for(int i=0; i<ylist.size; i++){
					vector y = ylist[i];
					double x = xlist[i];
					WriteLine($"{x} {y[0]} {y[1]}");
				}
			}
			if(run=="-B"){
				Func<double, vector, vector> f = (double x, vector y) => {
					/* u'' + u' = 1 + eps*u^2 */
					vector yprime = new vector(2);
					yprime[0] = y[1];
					yprime[1] = 1 - y[0];
					return yprime;
				};
				Func<double, vector, vector> feps = (double x, vector y) => {
					vector yprime = new vector(2);
					yprime[0] = y[1];
					yprime[1] = 1 - y[0] + 0.01*y[0]*y[0];
					return yprime;
				};
				genlist<double> xlist1, xlist2, xlist3;
				genlist<vector> ylist1, ylist2, ylist3;
				xlist1 = new genlist<double>(); xlist2 = new genlist<double>(); xlist3 = new genlist<double>();
				ylist1 = new genlist<vector>(); ylist2 = new genlist<vector>(); ylist3 = new genlist<vector>();
				rkint.driver(
					f: f,
					a: 0.0,
					ya: new vector(1.0,0.0),
					b: 100.0,
					xlist: xlist1,
					ylist: ylist1,
					maxstep:0.1);
				rkint.driver(
					f: f,
					a: 0.0,
					ya: new vector(0.5,-0.5),
					b: 100.0,
					xlist: xlist2,
					ylist: ylist2);
				rkint.driver(
					f: feps,
					a: 0.0,
					ya: new vector(0.5,-0.5),
					b: 100.0,
					xlist: xlist3,
					ylist: ylist3);
				var outstream = new StreamWriter("Out.B.txt",append:false);
				outstream.WriteLine("--- Solved ODEs for planetary motion ---");
				outstream.WriteLine($"No. of steps for Newtonian Circular Motion: {xlist1.size}");
				outstream.WriteLine("NOTE: due to constant function, max step was set to 0.1 to increase no. of plotting points");
				outstream.WriteLine($"No. of steps for Newtonian Elliptical Motion: {xlist2.size}");
				outstream.WriteLine($"No. of steps for relativistic precession: {xlist3.size}");
				outstream.WriteLine("These planetary motions are plotted in Planets.svg");
				outstream.Close();
				var outputstream = new StreamWriter("out.B1.data",append:false);
				for(int i=0; i<ylist1.size; i++){
					vector y = ylist1[i];
					double x = xlist1[i];
					outputstream.WriteLine($"{x} {y[0]}");
				}
				outputstream.Close();
				outputstream = new StreamWriter("out.B2.data",append:false);
				for(int i=0; i<ylist2.size; i++){
					vector y = ylist2[i];
					double x = xlist2[i];
					outputstream.WriteLine($"{x} {y[0]}");
				}
				outputstream.Close();
				outputstream = new StreamWriter("out.B3.data",append:false);
				for(int i=0; i<ylist3.size; i++){
					vector y = ylist3[i];
					double x = xlist3[i];
					outputstream.WriteLine($"{x} {y[0]}");
				}
				outputstream.Close();
				
			}
		}
	}
}
