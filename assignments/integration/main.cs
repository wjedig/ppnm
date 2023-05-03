using System;
using System.IO;
using static System.Console;
using static System.Math;
using static System.Double;
//using linalg;


public static class numint{
	public static double integrate(
			Func<double,double> f,	/* Function to be integrated */
			double a, 		/* start point */
			double b,		/* end point */
			ref double count,
			double d = 0.001,
			double eps = 0.001,
			double f2 = NaN,
			double f3 = NaN
			)
	{
		double h = b-a;
		if(!(IsNaN(count))) count+=1;
		if(IsNaN(f2)){f2 = f(a+2.0*h/6.0); f3=f(a+4.0*h/6.0);}
		double f1 = f(a+h/6), f4 = f(a+5*h/6);
		double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // higher order rule
		double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule
		double err = Abs(Q-q);
		if(err <= d+eps*Abs(Q)) return Q;
		else return integrate(f, a, (a+b)/2.0, ref count, d/Sqrt(2.0),eps,f1,f2) + 
			integrate(f, (a+b)/2.0, b, ref count, d/Sqrt(2),eps,f3,f4);
	}
	public static double integrate(			//overloaded method without counter
			Func<double,double> f,
			double a,
			double b,
			double d = 0.001,
			double eps = 0.001){
		double c = 0;
		return integrate(f,a,b,ref c,d,eps);
	}
			
	public static double erf(double z){
		if(z<0) return -erf(-z);
		if(1.0<z){
			Func<double,double> f = (double t) => Exp(-Pow(z+(1-t)/t,2))/t/t;
			return 1 - (2.0/Sqrt(PI))*integrate(f,0.0,1.0);
		}
		Func<double,double> f1 = (double x) => Exp(-x*x);
		return (2.0/Sqrt(PI))*integrate(f1,0.0,z);
	}
	public static double openccint(
			Func<double,double> f,	/* Function to be integrated */
			double a, 		/* start point */
			double b,		/* end point */
			ref double count,	
			double d = 0.001,
			double eps = 0.001
			){
		Func<double,double> F = (double theta) => f((a+b)/2 + (b-a)*Cos(theta)/2)*(b-a)*Sin(theta)/2;
		return integrate(F,0,PI,ref count,d,eps);
	}
	public static double openccint(			//overloaded method without counter
			Func<double,double> f,
			double a,
			double b,
			double d = 0.001,
			double eps = 0.001){
		double c = 0;
		return openccint(f,a,b,ref c,d,eps);
	}
}

class main{
	public static void Main(string[] args){
		foreach(string run in args){
			if(run=="-A"){
				var outstream = new StreamWriter("Out.A.txt",append:false);
				double eps = 0.001;
				Func<double,double> f = (double x) => Pow(x,0.5);
				double a = numint.integrate(f,0.0,1.0);
				string isAccurate = (Abs(a-(2.0/3.0)) < eps) ? "Yes!" : "No.";
				outstream.WriteLine($"integral 1 = {a} (expected: 2/3). Within accuracy? {isAccurate}");
				f = (double x) => 1.0/(Sqrt(x));
				a = numint.integrate(f,0.0,1.0);
				isAccurate = (Abs(a-2.0) < eps) ? "Yes!" : "No.";
				outstream.WriteLine($"integral 2 = {a} (expected: 2). Within accuracy? {isAccurate}");
				f = (double x) => 4*Sqrt(1-x*x);
				a = numint.integrate(f,0.0,1.0);
				isAccurate = (Abs(a-PI) < eps) ? "Yes!" : "No.";
				outstream.WriteLine($"integral 3 = {a} (expected: pi). Within accuracy? {isAccurate}");
				f = (double x) => Log(x)/(Sqrt(x));
				a = numint.integrate(f,0.0,1.0);
				isAccurate = (Abs(a+4.0) < eps) ? "Yes!" : "No.";
				outstream.WriteLine($"integral 4 = {a} (expected: -4). Within accuracy? {isAccurate}");
				outstream.WriteLine("\n--- Comparing calculated and tabulated values ---");
				outstream.WriteLine("\nx  erf(x)[tabulated]  erf(x)[calculated]  difference");
				var instream = new StreamReader("erf_tab.data");
				double err = 0;
				for(string line=instream.ReadLine(); line != null; line=instream.ReadLine()){
					string[] xy = line.Split(" ");
					double x = double.Parse(xy[0]);
					double tab_y = double.Parse(xy[1]);
					double calc_y = numint.erf(x);
					double diff = Abs(tab_y-calc_y);
					//outstream.WriteLine($"{diff}");
					if(diff>err) err = diff;
					outstream.WriteLine($"{x} {tab_y} {calc_y} {diff}");
				}
				instream.Close();
				outstream.WriteLine($"Biggest error: {err}");
				outstream.Close();
				double start = -5.0;
				double stop = 5.0;
				double n = 100;
				double spacing = (stop-start)/n;
				outstream = new StreamWriter("out.a.data",append:false);
				do{
					outstream.WriteLine($"{start} {numint.erf(start)}");
					start += spacing;
				}while(start<=stop);
				outstream.Close();
			}//A
			if(run=="-B"){
				Func<double,double> f = (double x) => 1/Sqrt(x);
				double c = 0; 
				var val = numint.integrate(f:f,a:0,b:1,count: ref c);
				WriteLine("--- Integral of (1/sqrt(x) from 0 to 1 ---)");
				WriteLine($"Integral count (classic): {c}");
				WriteLine($"Integral evaluation (classic): {val}");
				c = 0;
				val = numint.openccint(f:f,a:0,b:1,count: ref c);
				WriteLine($"Integral count (Clenshaw–Curtis): {c}");
				WriteLine($"Integral evaluation (Clenshaw–Curtis): {val}");
				WriteLine("For my scipy evaluation, integral count was 231");
				f = (double x) => Log(x)/Sqrt(x);
				c = 0; 
				val = numint.integrate(f:f,a:0,b:1,count: ref c);
				WriteLine("--- Integral of ln(x)/sqrt(x) from 0 to 1 ---)");
				WriteLine($"Integral count (classic): {c}");
				WriteLine($"Integral evaluation (classic): {val}");
				c = 0;
				val = numint.openccint(f:f,a:0,b:1,count: ref c);
				WriteLine($"Integral count (Clenshaw–Curtis): {c}");
				WriteLine($"Integral evaluation (Clenshaw–Curtis): {val}");
				WriteLine("For my scipy evaluation, integral count was 315");

			}
		}
	}
}
