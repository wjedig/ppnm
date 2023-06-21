using System;
using System.IO;
using static System.Console;
using static System.Math;
using static System.Random;
using linalg;

public static class mc{
	public static (double,double) plainmc(Func<vector,double> f,vector a,vector b,int N){
		int dim = a.size;
		double V = 1;
		for(int i=0; i<dim; i++) V*=b[i]-a[i];
		double sum = 0, sum2 = 0;
		vector x = new vector(dim);
		var rnd = new Random();
		for(int i=0; i<N; i++){
			for(int k=0; k<dim; k++) x[k] = a[k]+rnd.NextDouble()*(b[k]-a[k]);
			double fx = f(x);
			sum += fx;
			sum2 += fx*fx;
		}
		double mean = sum/N;
		double sigma = Sqrt(sum2/N - mean*mean);
		var result = (mean*V,sigma*V/Sqrt(N));
		return result;
	}
	static double corput(int n, int b){
		double q=0, bk = 1.0/b;
		while(n>0){
			q += (n % b) * bk;
			n /= b;
			bk /= b;
		} return q;
	}
	public static void halton(int n, vector x, int offset=0){
		int dim = x.size;
		int[] bases = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
		if(dim + offset > bases.Length) throw new ArgumentException("dimension too large to handle");
		for(int i=0; i<dim; i++) x[i] = corput(n,bases[i+offset]);
	}
	public static (double, double) quasimc(Func<vector,double> f, vector a, vector b, int N){
		int dim = a.size;
		double V = 1;
		for(int i=0; i<dim; i++) V*=b[i]-a[i];
		double sum = 0, sum2 = 0;
		vector x = new vector(dim);
		for(int i=0; i<N; i++){
			halton(i,x);
			double fx = f(x);
			sum += fx;
			halton(i,x,dim);
			fx = f(x);
			sum2 += fx;
		}
		double mean = sum/N, mean2 = sum2/N;
		var result = (mean*V,Abs(mean-mean2)*V);
		return result;
	}
}

class main{
	public static void Main(string[] args){
		foreach(string run in args){
			if(run=="-A"){
				/* Checking the area of a circle with radius 1 */
				Func<vector,double> f = (vector v) => (v.norm() <= 1.0) ? 1.0 : 0.0;
				vector a = new vector(-1.0,-1.0);
				vector b = new vector(1.0,1.0);
				double acc_area = PI;
				double area, pred_err, acc_err;
				var outstream = new StreamWriter("out.a.data",append:false);
				int n = 1;
				double stop = 2e5;
				int spacing = 1;
				do{
					(area, pred_err) = mc.plainmc(f,a,b,n);
					acc_err = Abs(acc_area-area);
					outstream.WriteLine($"{n} {pred_err} {acc_err}");
					spacing += 2;
					n += spacing;
				}while(n <= stop);
				outstream.Close();
				WriteLine($"--- Estimating area of circle with 1e7 points ---");
				WriteLine($"\nCalculation by 'throwing'");
				n = Convert.ToInt32(1e7);
				(area, pred_err) = mc.plainmc(f,a,b,n);
				WriteLine($"Area calculated to be : {area}");
				WriteLine($"Actual area: {acc_area}");
				WriteLine($"Error estimated to be : {pred_err}");
				WriteLine($"Actual error: {Abs(acc_area-area)}");
				WriteLine($"\nCalculation by evaluating function value with spherical coordinates");
				WriteLine($"This means evaluating f(r,phi) = r from r = phi = 0, to r = 1, phi = 2pi");
				f = (vector v) => v[0];
				a = new vector(0.0,0.0);
				b = new vector(1.0,2.0*PI);
				(area, pred_err) = mc.plainmc(f,a,b,n);
				WriteLine($"Area calculated to be : {area}");
				WriteLine($"Actual area: {acc_area}");
				WriteLine($"Error estimated to be : {pred_err}");
				WriteLine($"Actual error: {Abs(acc_area-area)}");
				WriteLine($"\n\n--- Estimating volume of unit sphere ---");
				f = (vector v) => v[0]*v[0]*Sin(v[2]);
				a = new vector(0.0,0.0,0.0);
				b = new vector(1.0,2.0*PI,PI);
				acc_area = (4.0/3.0)*PI;
				(area, pred_err) = mc.plainmc(f,a,b,n);
				WriteLine($"Area calculated to be : {area}");
				WriteLine($"Actual area: {acc_area} (4/3 pi)");
				WriteLine($"Error estimated to be : {pred_err}");
				WriteLine($"Actual error: {Abs(acc_area-area)}");
				WriteLine($"\n\n--- Estimating area of crazy function from hw with 1e7 points ---");
				f = (vector v) => {
					double fval = Pow(PI,-3)*(1.0/(1.0-Cos(v[0])*Cos(v[1])*Cos(v[2])));
					return fval;
				};
				a = new vector(0.0,0.0,0.0);
				b = new vector(PI,PI,PI);
				(area, pred_err) = mc.plainmc(f,a,b,n);
				acc_area = 1.3932039296856768591842462603255;
				WriteLine($"Area calculated to be : {area}");
				WriteLine($"Actual area: {acc_area}");
				WriteLine($"Error estimated to be : {pred_err}");
				WriteLine($"Actual error: {Abs(acc_area-area)}");
			}//A
			if(run=="-B"){
				WriteLine("--- Estimating the area of a circle with radius 1, with Quasirandom sequences ---");
				WriteLine("(I'm again using 1e7 points)");
				Func<vector,double> f = (vector v) => v[0];
				vector a = new vector(0.0,0.0);
				vector b = new vector(1.0,2.0*PI);
				double acc_area = PI;
				double area, pred_err;
				int n = Convert.ToInt32(1e7);
				(area, pred_err) = mc.quasimc(f,a,b,n);
				WriteLine($"Area calculated to be : {area}");
				WriteLine($"Actual area: {acc_area}");
				WriteLine($"Error estimated to be : {pred_err}");
				WriteLine($"Actual error: {Abs(acc_area-area)}");
				var outstream = new StreamWriter("out.b.data",append:false);
				n = 1;
				double stop = 2e5;
				int spacing = 1;
				int sspacing = 1;
				double area2, acc_err1, acc_err2;
				do{
					(area, pred_err) = mc.plainmc(f,a,b,n);
					(area2, _) = mc.quasimc(f,a,b,n);
					acc_err1 = Abs(acc_area-area);
					acc_err2 = Abs(acc_area-area2);
					outstream.WriteLine($"{n} {acc_err1} {acc_err2}");
					spacing += sspacing;
					sspacing += 1;
					n += spacing;
				}while(n <= stop);
				outstream.Close();
			}
		}
	}
}
