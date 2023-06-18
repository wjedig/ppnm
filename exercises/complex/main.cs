using System;
using static System.Console;
using static System.Math;

class main{
	public static void Main(){
		complex z;
		z = cmath.sqrt(-complex.One);
		WriteLine($"sqrt(-1) = {z} (expected: +/- i). {correct(z,0,-1)}");
		z = cmath.sqrt(complex.I);
		WriteLine($"sqrt(i) = {z} (expected: 1/√2+i/√2). {correct(z,1.0/Sqrt(2),1.0/Sqrt(2))}");
		z = cmath.exp(complex.I);
		WriteLine($"e^i = {z} (expected around: 0.54 + 0.84i). {correct(z,Cos(1),Sin(1))}");
		z = cmath.exp(complex.I*PI);
		WriteLine($"e^i*pi = {z} (expected: -1). {correct(z,-1,0)}");
		z = cmath.pow(complex.I,complex.I);
		WriteLine($"i^i = {z} (expected around: 0.208). {correct(z,Exp(-PI/2.0),0)}");
		z = cmath.log(complex.I);
		WriteLine($"ln(i) = {z} (expected: pi/2 * i). {correct(z,0,PI/2.0)}");
		//WriteLine($"sin(i*pi) = {cmath.log(complex.I)} (expected: pi/2 * i)");
	}
	static string correct(complex z, double re, double im){
		bool same = complex.approx(z.Re,re) && complex.approx(z.Im,im);
		return same ? "This result is approximately correct.\n" : "This result is not correct.\n";
	}
}
