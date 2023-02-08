using System;
using static System.Console;
using static System.Math;

class main{
	public static void Main(){
		double i1 = Sqrt(2);
		double i2 = Pow(2,1.0/5);
		double i3 = Pow(E,PI);
		double i4 = Pow(PI,E);
		Write($"sqrt2^2 = {i1*i1} (should equal 2)\n");
		Write($"(2^(1/5))^5 = {Pow(i2,5)}\n");
		Write($"(e^pi = {i3})\n");
		Write($"(pi^e = {i4})\n");
		WriteLine("****************************************");
		Write($"gamma(1): {sfuncs.gamma(1)}\n");
		Write($"gamma(2): {sfuncs.gamma(2)}\n");
		Write($"gamma(3): {sfuncs.gamma(3)}\n");
		Write($"gamma(10): {sfuncs.gamma(10)}\n");
	}
}
