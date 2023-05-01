using System;
using static System.Console;
using static System.Math;
using linalg;

public static class utils{
	public static int binsearch(vector x, double z){
		/* locates the interval for z by bisection */
		if(!(x[0]<=z && z<=x[x.size-1])) throw new Exception("Binsearch: z outside of x range");
		int i = 0, j = x.size-1;
		while(j-i>1){
			int mid = (i+j)/2;
			if(z>x[mid]) i = mid; else j = mid;
		}
		return i;
	}
}


public static class linspline{
	public static double linterp(vector x, vector y, double z, int i = -1){
		if(i == -1) i = utils.binsearch(x, z);
		double dx = x[i+1]-x[i];
		if(!(dx>0)) throw new Exception("dx must be greater than 0");
		double dy = y[i+1]-y[i];
		return y[i] + (dy/dx)*(z-x[i]);
	}
	public static double linterpInteg(vector x, vector y, double z){
		double sum = 0;
		int i = utils.binsearch(x, z);
		double y_z = linterp(x,y,z,i);
		for(int j=0; j<i; j++){
			sum += (x[j+1]-x[j])*(y[j]+y[j+1])/2;
		}
		sum += (z-x[i])*(y_z+y[i])/2;
		return sum;
	}
}

class qspline {
	vector x,y,b,c;
	public qspline(vector xs,vector ys){
		x=xs.copy(); y=ys.copy(); /* calculate b and c */
		int n = x.size;
		c = new vector(n-1);
		b = new vector(n-1);
		double[] p = new double[n-1], h = new double[n-1]; //p[n-1], h[n-1];
		for(int i=0 ; i<n-1 ; i++){
			h[i] = x[i+1] - x[i];
			p[i] = (y[i+1] - y[i])/h[i];
		}
		c[0] = 0;
		for(int i=0 ; i<n-2 ; i++){
			c[i+1] = (p[i+1]-p[i]-c[i]*h[i])/h[i+1];
		}
		c[n-2] /= 2;
		for(int i=n-3 ; i>=0 ; i--){
			c[i] = (p[i+1]-p[i]-c[i+1]*h[i+1])/h[i];
		}
		for(int i=0 ; i<n-1 ; i++){
			b[i] = p[i] - (c[i]*h[i]);
		}
	}
	public double evaluate(double z){/* evaluate the spline */
		int i = utils.binsearch(x,z);
		double s = y[i] + b[i]*(z-x[i]) + c[i]*Pow(z-x[i],2);
		return s;
	}
	public double derivative(double z){/* evaluate the derivative */	
		int i = utils.binsearch(x,z);
		double deriv = b[i] + 2*c[i]*(z-x[i]);
		return deriv;
	}
	public double integral(double z){/* evaluate the integral */	
		int i = utils.binsearch(x,z);
		double sum = 0;
		for(int j=0 ; j<i ; j++){
			double dx = x[i+1]-x[i];
			sum += y[j]*dx + b[i]*Pow(dx,2)/2 + c[i]*Pow(dx,3)/6;
		}	
		return sum;
	}
	public vector get_b(){return b;}
	public vector get_c(){return c;}
}

class main{
	public static void Main(string[] args){
		foreach(string run in args){
			if(run=="-A"){
				vector xsvec = new vector(11), ysvec = new vector(11);
				int i = 0;
				for(string line=In.ReadLine(); line != null; line=In.ReadLine()){
					string[] vals = line.Split(" ");
					xsvec[i] = double.Parse(vals[0]);
					ysvec[i] = double.Parse(vals[1]);
					//WriteLine($"{xsvec[i]} {ysvec[i]}");
					i++;
				}
				double integ = linspline.linterpInteg(xsvec, ysvec, 1.7);
				WriteLine($"integral = {integ}");
			}
			if(run=="-plot"){
				vector xsvec = new vector(11), ysvec = new vector(11);
				int i = 0;
				for(string line=In.ReadLine(); line != null; line=In.ReadLine()){
					string[] vals = line.Split(" ");
					xsvec[i] = double.Parse(vals[0]);
					ysvec[i] = double.Parse(vals[1]);
					//WriteLine($"{xsvec[i]} {ysvec[i]}");
					i++;
				}
				int n = 100;
				vector z_xs = new vector(n), z_ys = new vector(n);
				double dxz = (xsvec[xsvec.size-1]-xsvec[0])/n;
				for(i=0; i<n; i++){
					z_xs[i] = dxz*i;
					z_ys[i] = linspline.linterp(xsvec,ysvec,z_xs[i]);
					WriteLine($"{z_xs[i]} {z_ys[i]}");
				}
				
			}
			if(run=="-B"){
				vector xsvec = new vector("1,2,3,4,5"), ysvec = new vector(5);
				for(int i=0; i<5; i++) ysvec[i] = 1;
				qspline testspline = new qspline(xsvec,ysvec);
				vector b = testspline.get_b();
				vector c = testspline.get_c();
				WriteLine("For x = {1,2,3,4,5}, y = {1,1,1,1,1}");
				b.print("b = ");
				c.print("c = ");
				for(int i=0; i<5; i++) ysvec[i] = i+1;
				testspline = new qspline(xsvec, ysvec);
				b = testspline.get_b();
				c = testspline.get_c();
				WriteLine("For x = {1,2,3,4,5}, y = {1,2,3,4,5}");
				b.print("b = ");
				c.print("c = ");
				for(int i=0; i<5; i++) ysvec[i] = Pow(i+1,2);
				testspline = new qspline(xsvec, ysvec);
				b = testspline.get_b();
				c = testspline.get_c();
				WriteLine("For x = {1,2,3,4,5}, y = {1,4,9,16,25}");
				b.print("b = ");
				c.print("c = ");
				WriteLine("These are also the results I got manually :D");	

			}
		// int x = Splines.binsearch(new vector(0.5,0.7,0.9,1.1), 0.8);
		// WriteLine(x);
		}
	}
}
