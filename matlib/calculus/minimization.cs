using System;
using static System.Console;
using static System.Math;
using linalg;

namespace calculus{
public static class min{
	public static double minlambda = Pow(2,-20);
	public static int count1 = 0;
	public static int count2 = 0;
	public static double maxcount = 1e7;
	public static vector qnewton(
		Func<vector,double> f, /* objective function */
		vector start, /* starting point */
		double acc = 1e-6, /* accuracy goal, on exit |gradient| should be < acc */
		double eps = 1e-20,
		bool trace = false
	){
		count1 = 0;
		count2 = 0;
		vector x = start.copy();
		int n = x.size;
		matrix B = matrix.id(n);
		vector grad, newtonstep;
		vector newgrad = jacobi.numderiv(f,x);
		double lambda;
		do{
			grad = newgrad;
			newtonstep = -(B*grad);
			lambda = 1.0;
			count1 += 1;
			if(trace) WriteLine($"Iteration no.: {count1}");
			if(count1 > maxcount) throw new Exception("Max iterations reached! Couldn't converge");
			if(trace){x.print("current x: "); grad.print("current grad: ");}
			do{
				count2 += 1;
				if(f(x+lambda*newtonstep) < f(x)){
					vector s = lambda * newtonstep;
					x = x + s;
					newgrad = jacobi.numderiv(f, x);
					vector y = newgrad - grad;
					vector u = s - B*y;
					double uTy = u%y;
					if(Abs(uTy) < eps) break;
					matrix uuT = matrix.outer(u,u);
					matrix dB = uuT/uTy;
					B = B + dB;
					break;
				}
				if(trace) WriteLine($"Current lambda: {lambda}");
				lambda /= 2.0;
				if(lambda < minlambda){
					x = x + lambda*newtonstep;
					newgrad = jacobi.numderiv(f,x);
					B.set_unity();
					break;
				}
			}while(true);
		}while(grad.norm() > acc);
		if(trace) grad.print("Final gradient: ");
		return x;
	}//qnewton
	public static double qnewton(
		Func<double,double> f, /* objective function */
		double start, /* starting point */
		double acc = 1e-6, /* accuracy goal, on exit |gradient| should be < acc */
		double eps = 1e-6,
		bool trace = false
	){
		vector startvec = new vector(1); startvec[0] = start;
		Func<vector,double> f2 = (vector xvec) => f(xvec[0]);
		return qnewton(f:f2,start:startvec,acc:acc,eps:eps,trace:trace)[0];
	}//qnewton2
}//min
}//namespace
