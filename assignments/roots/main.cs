using System;
using static System.Console;
using static System.Math;
using System.IO;
using linalg;
using calculus;
using genlist;

public static class root{
	public static vector newton(Func<vector,vector> f, vector x, double eps=1e-2, bool trace=false){
		Func<vector,double,double> thresh = (vector x1, double l1) => (1-(l1/2.0))*(f(x1).norm());
		Func<vector,vector,double,double> val = (vector x2, vector delx, double l2) => f(x2+(l2*delx)).norm();
		while(f(x).norm() > eps){
			//WriteLine("evaluating J");
			matrix J = jacobi.numderiv(f,x);
			vector b = -f(x);
			//WriteLine("evaluating delta x");
			vector deltax = QRGS.solve(J,b);
			double lambda = 1.0;
			//WriteLine($"Current value of f: {f(x).norm()}");
			do{
				//WriteLine($"Evaluating at lambda = {lambda}");
				lambda /= 2.0;
			} while( (val(x,deltax,lambda) > thresh(x,lambda)) && lambda >= (1.0/1024.0) );
			x += lambda*deltax;
		}
		return x;
	}
	public static vector newton(Func<vector,vector> f, double x, double eps=1e-2){
		vector xvec = new vector(1);
		xvec[0] = x;
		return newton(f,xvec,eps);
	}
}

class main{
	public static int counter;
	public static int maxrep = 1000;
	public static double rmax = 10;
	public static double rmin = 0.0001;
	public static double eps = 0.01;
	public static double acc = 0.01;

	public static void Main(string[] args){
		foreach(string run in args){
			if(run == "-test"){
				Func<vector,vector> f = (vector x) =>{
					vector y = new vector(3);
					y[0] = x[0]+x[1]; y[1] = x[0]*x[0]; y[2] = x[1]*x[1];
					return y;
				};
				WriteLine("--- Testing Jacobi Matrix for function f(x1,x2) ---\n");
				WriteLine("f1(x) = x1+x2");
				WriteLine("f2(x) = x1^2");
				WriteLine("f3(x) = x2^2");
				vector x1 = new vector (2.0,4.0);
				x1.print("\nTesting for x = ");
				matrix A = jacobi.numderiv(f,x1);
				A.print("\nJ = ");
			}//test
			if(run == "-A"){		
				A();
			}//A
			if(run=="-B:test"){
				int c = 0;
				double Emin = B();
				WriteLine($"Emin found to be {Emin} in {c} iterations");
			}
			if(run=="-B:plot"){
				var outstream = new StreamWriter("out.b.data", append:false);
				double Emin;
				outstream.WriteLine("Emin rmax rmin eps acc");
				rmax = 8.0; rmin = 0.01; eps = 0.01; acc = 0.01;
				for(double i=0; i<10; i++){
					WriteLine($"Progress: {i*10}%");
					rmax = 3.0 + (0.5*i);
					Emin = B();
					outstream.Write($"{Emin} {rmax} ");
					rmax = 8.0; 
					rmin = 0.02 + (i*0.02);
					Emin = B();
					outstream.Write($"{Emin} {rmin} ");
					rmin = 0.01;
					eps = 0.1 + (i*0.1);
					Emin = B();
					outstream.Write($"{Emin} {eps} ");
					eps = 0.01;
					acc = 0.1 + (i*0.1);
					Emin = B();
					outstream.Write($"{Emin} {acc}\n");
					acc = 0.1;
				}
				outstream.Close();
			}
		}
	}//Main
	
	public static void A(){
		Func<vector,vector> f = (vector x) =>{		
			vector y = new vector(2);
			y[0] = x[0]*x[0]-16;
			y[1] = x[0]+4;
			return y;
		};
		vector xguess = new vector(-10.0);
		vector x1 = root.newton(f,xguess);
		WriteLine("--- Testing Jacobi Matrix for function f(x) ---\n");
		WriteLine("f1(x) = x^2-16");
		WriteLine("f2(x) = x+4");
		xguess.print("\nInitial guess: ");
		x1.print("Root found at: ");
		WriteLine("Actual root: (-4)");
		f = (vector x) =>{
			vector y = new vector(3);
			y[0] = x[0]+x[1]-1;
			y[1] = x[0]*x[0]-4;
			y[2] = x[1]*x[1]-9;
			return y;
		};
		xguess = new vector(-1.0,1.0);
		x1 = root.newton(f,xguess);
		WriteLine("\n--- Testing Jacobi Matrix for function f(x1,x2) ---\n");
		WriteLine("f1(x) = x1+x2-1");
		WriteLine("f2(x) = x1^2-4");
		WriteLine("f3(x) = x2^2-9");
		xguess.print("\nInitial guess: ");
		x1.print("Root found at: ");
		WriteLine("Actual root: (-2,3)");

		f = (vector x) =>{
			vector y = new vector(2);
			y[0] = x[2]*x[2]+x[1]-x[0];
			y[1] = Sqrt(x[0])-x[1]+3.0;
			return y;
		};
		//xguess = new vector(15.0,6.0,10.0);
		//x1 = root.newton(f,xguess);
		//WriteLine("\n--- Testing Jacobi Matrix for function f(x1,x2) ---\n");
		//WriteLine("f1(x) = x3^2+x2-x1");
		//WriteLine("f2(x) = Sqrt(x1)-x2+3");
		//xguess.print("\nInitial guess: ");
		//x1.print("Root found at: ");
		//WriteLine("Actual root: (16,7,3)");
		WriteLine("\n--- Finding extrema of Rosenbruck's valley function---");
		WriteLine("Derivative found to be: [-2(1-x)-400x(y-x^2)  ,  200(y-x^2)]");
		f = (vector x) => {
			vector y = new vector(2);
			y[0] = -2.0*(1.0-x[0]) - 400.0*x[0]*(x[1]-x[0]*x[0]);
			y[1] = 200.0*(x[1]-x[0]*x[0]);
			return y;
		};
		xguess = new vector(0.5,0.5);
		eps = 0.001; acc = 0.001;
		x1 = root.newton(f,xguess);
		xguess.print("\nInitial guess: ");
		x1.print("Root found at: ");
		WriteLine("Actual root: (1.0,1.0)");
		xguess = new vector(-0.9,0.9);
		eps = 0.001; acc = 0.001;
		x1 = root.newton(f,xguess);
		xguess.print("\nInitial guess: ");
		x1.print("Root found at: ");
		WriteLine("Actual root: (1.0,1.0)");
	}

	public static vector M(vector Evec){
		counter ++;
		double E = Evec[0];
		vector yval = new vector(1);
		vector finit = new vector(rmin-rmin*rmin, 1.0 - 2.0*rmin);
		vector yfinal = new vector(1); yfinal[0] = finit[0];
		if(counter > maxrep){
			yval[0] = yfinal[0];
			WriteLine($"Not converged for rmax:{rmax}, rmin:{rmin}, eps:{eps}, acc:{acc}");
			return yval;
		}
		Func<double,vector,vector> f = (double r, vector y) => {
			vector yprime = new vector(2);
			yprime[0] = y[1];
			yprime[1] = -2.0*(E + 1.0/r)*y[0];
			//WriteLine($"y = {y[0]}");
			return yprime;
		};
		yfinal = rkint.driver(
				f: f,
				a: rmin,
				ya: finit,
				b: rmax,
				acc: acc,
				eps: eps);
		yval[0] = yfinal[0];
		return yval;
	}

	public static double B(bool trace = false){
		counter = 0;
		vector Eguess = new vector(-0.7);
		vector Emin = root.newton(M,Eguess);
		if(trace) WriteLine($"Converged in {counter} iterations");
		return Emin[0];
	}
}//main

