using System;
using static System.Console;
using static System.Math;
class main{
	public static void Main(){
	genlist<double> listd = new genlist<double>();
	listd.add(1.0);
	listd.add(2.0);
	listd.add(3.0);
	Func<double,double> f;
	f = Sin;
	f = delegate(double x){return x*x;}; // Old notation
	f = (double x) => x*x;
	double y = f(2.0);
	WriteLine($"f(2.0) = 2*2 = {y}");
	}
}
