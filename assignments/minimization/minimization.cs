using System;
using static System.Console;
using static System.Math;
using linalg;

public static class min{
	public static int count1 = 0;
	public static int count2 = 0;
	public static vector qnewton(
		Func<vector,double> f, /* objective function */
		vector start, /* starting point */
		double acc = 1e-6, /* accuracy goal, on exit |gradient| should be < acc */
		double eps = 1e-6,
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
			if(trace) x.print("current x: ");
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
				//if(trace) WriteLine($"Current lambda: {lambda}");
				lambda /= 2.0;
				if(lambda < 1.0/1024.0){ // changed from 1.0/1024.0
					x = x + lambda*newtonstep;
					newgrad = jacobi.numderiv(f,x);
					B.set_unity();
					break;
				}
			}while(true);
		}while(grad.norm() > acc);
		if(trace) grad.print("Final gradient: ");
		return x;
	}
}
