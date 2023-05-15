using System;
using static System.Console;
using static System.Math;
using System.IO;
using static System.Random;
using linalg;
using genlist;
using calculus;

class main{
	public static void Main(string[] args){
		foreach(string run in args){
			if(run=="-A"){
				double acc = 1e-6;
				bool trace = true;
				var ann1 = new ann(5,trace:trace,acc:acc);
				int n = 200;
				var rand = new Random();
				vector x = new vector(n), y = new vector(n);
				var outstream = new StreamWriter("out.atrain.data",append:false);
				Func<double,double> g = (double a) => Cos(5.0*a-1.0)*Exp(-a*a);
				for(int i=0; i<n; i++){
					x[i] = rand.NextDouble()*2.0 - 1.0; /* n samples in [-1,1] */
					y[i] = g(x[i]);
					outstream.WriteLine($"{x[i]} {y[i]}");
				}
				outstream.Close();
				ann1.fit(x,y);
				vector p = ann1.get_params();
				p.print("params given by: ");
				int count = ann1.get_count();
				WriteLine($"Convergence took {count} iterations");
				outstream = new StreamWriter("out.afit.data",append:false);
				int npoints = 100;
				vector xfit = new vector(npoints), yfit = new vector(npoints);
				for(int i=0; i<npoints; i++){
					xfit[i] = (2.0*i-npoints)/npoints; /* n samples in [-1,1] */
					yfit[i] = ann1.predict(xfit[i]);
					outstream.WriteLine($"{xfit[i]} {yfit[i]}");
				}
				outstream.Close();
			}
		}
	}
}
