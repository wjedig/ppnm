using System;
using System.IO;
using static System.Console;
using static System.Math;
using System.Collections.Generic;
using linalg;

public static class Utils{
	public static void writeFitParams(vector a){
		string[] coeffs = {"a", "b", "c"};
		for(int i=0 ; i<a.size ; i++){
			WriteLine($"{coeffs[i]} = {a[i]}");
		}
	}
	public static void printBlank(){
		Write("\n *************************** \n\n");
	}
	public static double poly2fit(vector c, double x){
		var fs = new Func<double,double>[] { z => 1.0, z => z, z => z*z };
		return c[0]*fs[0](x) + c[1]*fs[1](x) + c[2]*fs[2](x);
	}
	public static (vector xs, vector ys, vector dys) readFile(string infile){
		var inputstream = new StreamReader(infile);
		List<double> xsl = new List<double>(), ysl = new List<double>(), dysl = new List<double>();
		for(string line=inputstream.ReadLine(); line != null; line=inputstream.ReadLine()){
			string[] words = line.Split(" ");
			xsl.Add(double.Parse(words[0]));
			ysl.Add(double.Parse(words[1]));
			dysl.Add(double.Parse(words[2]));
		}
		inputstream.Close();
		vector xs, ys, dys;
		xs = new vector(xsl.ToArray());
		ys = new vector(ysl.ToArray());
		dys = new vector(dysl.ToArray());
		return (xs, ys, dys);
	}
}


public static class Test{
	public static void testA(bool writeFit, bool writeData, bool trace, int n = 10){
		var fs = new Func<double,double>[] { z => 1.0, z => z, z => z*z };
		var rnd = new System.Random(1);
		vector xs = new vector(n);
		vector ys = new vector(n);
		vector dy = new vector(n);
		for(double i=0 ; i<n ; i++){
			xs[(int)i] = i;
			ys[(int)i] = i*i + i*rnd.NextDouble() - rnd.NextDouble();
			dy[(int)i] = rnd.NextDouble();
		}	
		(vector c, _) = LS.lsfit(fs,xs,ys,dy);
		if(writeFit){
			double start = xs[0];
			double stop = xs[n-1];
			double spacing = (stop-start)/100.0;
			for(double i=start ; i<=stop ; i+=spacing){
				WriteLine($"{i} {Utils.poly2fit(c,i)}");			
			}
		}
		if(writeData){
			for(int i=0; i<n; i++){
				WriteLine($"{xs[i]} {ys[i]} {dy[i]}");
			}
		}
		if(trace){
			xs.print("x = ");
			ys.print("y = ");
			dy.print("dy = ");
			c.print("best fit for a + bx + cx^2");
		}
	}
	public static (vector, matrix) fitDecay(vector xs, vector ys, vector dys, bool writeFit, bool writeData){
		var fs = new Func<double,double>[] { z => 1.0, z => -z};
		int n = xs.size;
		vector lny = new vector(n), dlny = new vector(n);
		for(int i=0 ; i<n ; i++){
			lny[i] = Log(ys[i]);
			dlny[i] = dys[i]/ys[i];
		}
		(vector c, matrix S) = LS.lsfit(fs,xs,lny,dlny);
		if(writeFit){
			Utils.writeFitParams(c);	
		}
		if(writeData){
			for(int i=0 ; i<n ; i++){
				WriteLine($"{xs[i]} {lny[i]} {dlny[i]}");
			}
		}
		return (c,S);
		
	}
}

class main{
	public static void Main(string[] args){
		bool writeFit = false;
		bool writeData = false;
		bool trace = false;
		string infile = null;
		bool writeCov = false;
		foreach(string arg in args){
			string[] words = arg.Split(":");
			if(words[0] == "-writeFit"){writeFit = true;}
			if(words[0] == "-writeData"){writeData = true;}
			if(words[0] == "-trace"){trace = true;}
			if(words[0]=="-input") infile = words[1];
			if(words[0] == "-testB"){writeCov = true;}
		}
		if(infile!=null){
			(vector xs, vector ys, vector dys) = Utils.readFile(infile);	
			(vector c, matrix S) = Test.fitDecay(xs,ys,dys,writeFit,writeData);
			if(writeCov){
				WriteLine("Fit parameters given by: ");
				Utils.writeFitParams(c);
				WriteLine("(Here b refers to lambda)");
				S.print("Covariance matrix given by:");
				Utils.printBlank();
				WriteLine("Half life given by:");
				double err = S[1][1];
				WriteLine($"{Log(2)/c[1]} +/- {err*Log(2)/(Pow(c[1],2))} days");
				WriteLine("Compared to the theoretical value of 3.6 days");
			}
			//xs.print("xs = ");
			//ys.print("ys = ");
			//dys.print("dys = ");
		}
		foreach(string run in args){
			if(run=="-test"){Test.testA(writeFit,writeData,trace);}
		}
	}
}
