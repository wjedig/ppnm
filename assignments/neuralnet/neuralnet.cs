using System;
using static System.Math;
using linalg;
using genlist;
using calculus;

public class ann{
	int n; /* number of hidden neurons */
	Func<double,double> f = (double x) => x*Exp(-x*x); /* activation function */
	vector p; /* network parameters. Set up so p[0:3] are for the 1st neuron, p[3:6] for the 2nd, etc */
	int count; /* Iterations before convergence */
	bool trace;
	vector x, y;
	double acc, eps;
   	public ann(int n,bool trace=false, double acc=1e-6, double eps=1e-20){
		this.n = n;
		this.trace = trace;
		this.acc = acc;
		this.eps = eps;
		this.p = new vector(n*3); /* 3 parameters per neuron*/
		for(int i=0; i<p.size; i++) p[i] = 1.0; /* start guess: all params = 1 */
	}
	public void fit(vector x, vector y){
		this.x = x; this.y = y;
		int N = x.size;
		if(N != y.size) throw new ArgumentException("input and output list must be same size");
		//Func<vector,double> cost = (vector fitp) => {
		//	double sum = 0.0;
		//	for(int i=0; i<N; i++){
		//		sum += Pow(F(x[i],fitp)-y[i],2);
		//	}
		//	return sum;
		//};
		vector newp = min.qnewton(Cost,p,acc:this.acc,eps:this.eps,trace:this.trace);
		this.count = min.count1;
		this.p = newp;
	}
	public double predict(double x){
		return F(x,p);
	}
	double Cost(vector fitp){
		double sum = 0.0;
		for(int i=0; i<x.size; i++){
			sum += Pow(F(x[i],fitp)-y[i],2);
		}
		return sum;
	}
	double F(double xi, vector fitp){
		double sum = 0.0;
		double ai, bi, wi;
		for(int i=0; i<p.size; i+=3){
			ai = fitp[i]; bi = fitp[i+1]; wi = fitp[i+2];
			sum += f((xi-ai)/bi)*wi;
		}
		return sum;
	}
	public void set_activation(Func<double,double> f){
		this.f = f;
	}
	public vector get_params(){
		return p;
	}
	public int get_count(){
		return count;
	}
}
