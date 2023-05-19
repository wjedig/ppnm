using System;
using static System.Math;
using static System.Random;
using linalg;
using genlist;
using calculus;

public class ann{
	int n; /* number of hidden neurons */
	static Func<double,double> Gauss = (double x) => Exp(-x*x);
	Func<double,double> f = (double x) => x*Gauss(x); /* activation function, Gaussian Wavelet */
	Func<double,double> fp = (double x) => (1-2*x*x)*Gauss(x); /* f'(x) */
	Func<double,double> fpp = (double x) => (4*x*x*x - 6*x)*Gauss(x); /* f''(x) */
	Func<double,double> fa = (double x) => (-1.0/2.0)*Gauss(x); /* Antiderivative F(x) */
	bool gaussWave = true;
	//Func<double,double> f = (double x) => Exp(-x*x); /* activation function */
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
		this.p = new vector(this.n*3); /* 3 parameters per neuron*/
		var rand = new Random(1);
		for(int i=0; i<p.size; i++) p[i] = rand.NextDouble(); /* start guess: all params = 1 */
	}
	public void fit(vector x, vector y){
		this.x = x; this.y = y;
		for(int i=0; i<n; i++) p[3*i] = x[0]+(i+0.5)*(x[x.size-1]-x[0])/n;
		//x.fprint(Console.Error,"\nx=\n");
		//p.fprint(Console.Error,"\np=\n");
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
	public double predict(double x, string func = "f"){
		if(gaussWave){
			if(func=="fp") return F(x,p,this.fp);
			if(func=="fpp") return F(x,p,this.fpp);
			if(func=="fa") return F(x,p,this.fa) - F(0.0,p,this.fa);
		}
		if(func == "f") return F(x,p,this.f);
		throw new Exception("f has been updated, analytical derivatives and antiderivatives unknown");
	}
	double Cost(vector fitp){
		double sum = 0.0;
		for(int i=0; i<x.size; i++){
			sum += Pow(F(x[i],fitp,this.f)-y[i],2);
		}
		return sum;
	}
	double F(double xi, vector fitp, Func<double,double> func){
		double sum = 0.0;
		double ai, bi, wi;
		for(int i=0; i<p.size; i+=3){
			ai = fitp[i]; bi = fitp[i+1]; wi = fitp[i+2];
			sum += func((xi-ai)/bi)*wi;
		}
		return sum;
	}
	public void set_activation(Func<double,double> f){
		this.f = f;
		this.gaussWave = false;
	}
	public vector get_params(){
		return p;
	}
	public int get_count(){
		return count;
	}
}
