using System;
using System.IO;
using static System.Console;
using static System.Math;
using linalg;
using genlist;
using calculus;

class main{
	public static genlist<double> energy;
	public static genlist<double> signal;
	public static genlist<double> error;

	public static void Main(string[] args){
		foreach(string run in args){
			if(run=="-A"){
				Func<vector,double> f = (vector x) => x[0]*x[0] - 10*x[0]; /* x^2 - 10x */
				vector start = new vector(-10.0);
				WriteLine("Finding minimum for x^2-10");
				start.print("Initial guess: ");
				vector stop = min.qnewton(f,start);
				stop.print("Minimum found at: ");
				WriteLine("Expected minimum at x = 5");
				WriteLine($"Minimization took {min.count1} iterations ({min.count2} backtracking iterations).");
				
				WriteLine("\nFinding a minimum for Rosenbrock's valley function");
				f = (vector x) => Pow(1.0-x[0],2) + 100.0*Pow(x[1]-x[0]*x[0],2);
				start = new vector(0.0,0.0);
				stop = min.qnewton(f,start);
				start.print("Initial guess: ");
				stop.print("Minimum found at: ");
				WriteLine("Expected minimum at (1,1)");
				WriteLine($"Minimization took {min.count1} iterations ({min.count2} backtracking iterations).");
				
				WriteLine("\nFinding a minimum for Himmelblau's function");
				f = (vector x) => Pow(x[0]*x[0]+x[1]-11.0,2) + Pow(x[0]+x[1]*x[1]-7,2);
				start = new vector(0.0,0.0);
				stop = min.qnewton(f,start);
				start.print("Initial guess: ");
				stop.print("Minimum found at: ");
				WriteLine("Expected minimum at (3,2) (there are 4 total minima)");
				WriteLine($"Minimization took {min.count1} iterations ({min.count2} backtracking iterations).");

			}//A
			if(run=="-B"){
				min.minlambda = 1.0/1024.0;
				energy = new genlist<double>();
				signal = new genlist<double>();
				error  = new genlist<double>();
				var separators = new char[] {' ','\t'};
				var options = StringSplitOptions.RemoveEmptyEntries;
				Console.In.ReadLine(); /* Removing columns from data */
				do{
					string line=Console.In.ReadLine();
					if(line==null)break;
					string[] words=line.Split(separators,options);
					energy.add(double.Parse(words[0]));
					signal.add(double.Parse(words[1]));
					error.add(double.Parse(words[2]));
				}while(true);
				//WriteLine($"EList size: {energy.size}, Slist size: {signal.size}, errList size: {error.size}");

				vector guess = new vector (120.0,10.0,10.0);
				guess.print("Initial guess: ");
				WriteLine("Originally, this was run with (10.0,10.0,10.0), which worked...");
				WriteLine("But took a little too long (~45 seconds)");
				vector fit = min.qnewton(D,guess);
				fit.print("Parameters found to be: ");
				WriteLine($"Minimization took {min.count1} iterations ({min.count2} backtracking iterations).");
				var outstream = new StreamWriter("out.bparams.par",append:false);
				outstream.WriteLine($"m = {fit[0]}");
				outstream.WriteLine($"G = {fit[1]}");
				outstream.WriteLine($"A = {fit[2]}");
				outstream.Close();
			}
		}
	}
	public static double D(vector parameters){
		Func<double,vector,double> F = (double E, vector pars) => {
			double m = pars[0], G = pars[1], A = pars[2];
			return A/(Pow(E-m,2)+(G*G)/4);
		};
		double sum = 0.0;
		for(int i=0; i<energy.size; i++){
			sum += Pow((F(energy[i],parameters)-signal[i])/error[i],2);
		}
		//WriteLine($"D(E) = {sum}");
		return sum;
	}
}
