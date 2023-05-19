using System;
using static System.Console;
using static System.Math;
using System.IO;
using static System.Random;
using linalg;
using genlist;
using calculus;

class main{
	public static bool trace = false;
	public static double acc = 1e-2;
	public static int nx = 50;
	public static void Main(string[] args){
		foreach(string arg in args){
			string[] words = arg.Split(':');
			if(words[0] == "-trace") trace = bool.Parse(words[1]);
			if(words[0] == "-acc")	acc = double.Parse(words[1]);
			if(words[0] == "-nx")	nx = int.Parse(words[1]);
		}
		vector x = new vector(nx), y = new vector(nx);
		var outstream = new StreamWriter("out.atrain.data",append:false);
		Func<double,double> g = (double a) => Cos(5.0*a-1.0)*Exp(-a*a);
		for(int i=0; i<nx; i++){
			x[i] = -1 + i*2.0/(nx-1); /* n samples in [-1,1] */
			y[i] = g(x[i]);
			outstream.WriteLine($"{x[i]} {y[i]}");
		}
		outstream.Close();
		foreach(string run in args){
			if(run=="-A"){
				var ann1 = new ann(1,trace:trace,acc:acc);
				ann1.fit(x,y);
				vector p = ann1.get_params();
				WriteLine("\n--- ANN with 1 neurons ---");
				WriteLine("Parameters given by:");
				for(int i=0; i<p.size/3; i++){
					Write($"a[{i}] = {p[i]}, ");
					Write($"b[{i}] = {p[i+1]}, ");
					WriteLine($"w[{i}] = {p[i+2]}");
				}
				int count = ann1.get_count();
				WriteLine($"Convergence took {count} iterations for 1 neuron");
				outstream = new StreamWriter("out.afit1.data",append:false);
				int npoints = 100;
				vector xfit = new vector(npoints), yfit = new vector(npoints);
				for(int i=0; i<npoints; i++){
					xfit[i] = (2.0*i-npoints)/npoints; /* n samples in [-1,1] */
					yfit[i] = ann1.predict(xfit[i]);
					outstream.WriteLine($"{xfit[i]} {yfit[i]}");
				}
				outstream.Close();

				var ann10 = new ann(10,trace:trace,acc:acc);
				ann10.fit(x,y);
				p = ann10.get_params();
				WriteLine("\n--- ANN with 10 neurons ---");
				WriteLine("Parameters given by:");
				for(int i=0; i<p.size/3; i++){
					Write($"a[{i}] = {p[i]}, ");
					Write($"b[{i}] = {p[i+1]}, ");
					WriteLine($"w[{i}] = {p[i+2]}");
				}
				count = ann10.get_count();
				WriteLine($"Convergence took {count} iterations for 10 neurons");
				outstream = new StreamWriter("out.afit2.data",append:false);
				for(int i=0; i<npoints; i++){
					yfit[i] = ann10.predict(xfit[i]);
					outstream.WriteLine($"{xfit[i]} {yfit[i]}");
				}
				outstream.Close();
			}//A
			if(run=="-B"){
				WriteLine("--- PART B ---");
				WriteLine("\nTesting for g(x) = 4x*e^(-x^2)");
				int neurons = 5;
				var ann5 = new ann(neurons,trace:trace,acc:acc);
				outstream = new StreamWriter("out.btrain1.data",append:false);
				nx = 100;
				x = new vector(nx); y = new vector(nx);
				Func<double,double> f = (double a) => 4*a*Exp(-a*a);
				for(int i=0; i<nx; i++){
					x[i] = -1 + i*2.0/(nx-1); /* n samples in [-1,1] */
					y[i] = f(x[i]);
					outstream.WriteLine($"{x[i]} {y[i]}");
				}
				outstream.Close();
				ann5.fit(x,y);
				Func<double,double> Gwave = (double a) => 4*Exp(-a*a);
				Func<double,double> fp = (double a) => (1-2*a*a)*Gwave(a);
				Func<double,double> fpp = (double a) => a*(4*a*a - 6)*Gwave(a);
				Func<double,double> F = (double a) => (-1.0/2)*Gwave(a);
				Func<double,double> fa = (double a) => F(a) - F(0);
				int count = ann5.get_count();
				WriteLine($"Convergence took {count} iterations for {neurons} neurons");
				outstream = new StreamWriter("out.bfit1.data",append:false);
				int npoints = 100;
				double xfit, yfit, ypfit, ypexp, yppfit, yppexp, yafit, yaexp;
				for(int i=0; i<npoints; i++){
					xfit = (2.0*i-npoints)/npoints; /* n samples in [-1,1] */
					yfit = ann5.predict(xfit,"f");
					ypfit = ann5.predict(xfit,"fp");
					ypexp = fp(xfit);
					yppfit = ann5.predict(xfit,"fpp");
					yppexp = fpp(xfit);
					yafit = ann5.predict(xfit,"fa");
					yaexp = fa(xfit);
					outstream.WriteLine($"{xfit} {yfit} {ypfit} {ypexp} {yppfit} {yppexp} {yafit} {yaexp}");
				}
				outstream.Close();
				f = (double a) => a*a;
				WriteLine("\nTesting for g(x) = x^2");
				WriteLine("This should make the g'(x) = 2x, g''(x) = 2, and G(x) = 1/2 x^3...");
				WriteLine("... which should be easy to confirm visually on ANN_2.svg");
				outstream = new StreamWriter("out.btrain2.data",append:false);
				nx = 100;
				x = new vector(nx); y = new vector(nx);
				for(int i=0; i<nx; i++){
					x[i] = -1 + i*2.0/(nx-1); /* n samples in [-1,1] */
					y[i] = f(x[i]);
					outstream.WriteLine($"{x[i]} {y[i]}");
				}
				outstream.Close();
				ann5 = new ann(neurons,trace:trace,acc:acc);
				ann5.fit(x,y);
				count = ann5.get_count();
				WriteLine($"Convergence took {count} iterations for {neurons} neurons");
				outstream = new StreamWriter("out.bfit2.data",append:false);
				npoints = 100;
				fp = (double a) => 2*a;
				fpp = (double a) => 2.0;
				F = (double a) => 0.5*a*a*a;
				for(int i=0; i<npoints; i++){
					xfit = (2.0*i-npoints)/npoints; /* n samples in [-1,1] */
					yfit = ann5.predict(xfit,"f");
					ypfit = ann5.predict(xfit,"fp");
					ypexp = fp(xfit);
					yppfit = ann5.predict(xfit,"fpp");
					yppexp = fpp(xfit);
					yafit = ann5.predict(xfit,"fa");
					yaexp = F(xfit);
					outstream.WriteLine($"{xfit} {yfit} {ypfit} {ypexp} {yppfit} {yppexp} {yafit} {yaexp}");
				}
				outstream.Close();
			}
		}//foreach
	}//Main
}//main
