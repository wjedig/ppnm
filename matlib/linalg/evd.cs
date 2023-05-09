using System;
using static System.Console;
using static System.Math;

namespace linalg{

public static class jacobi{
	public static matrix numderiv(Func<vector,vector> f, vector x){ 
		vector fx = f(x);
		int n = fx.size, m = x.size;
		double dx = x.norm()*Pow(2.0,-26);
		if(dx == 0) dx = Pow(2.0,-26);
		matrix J = new matrix(n,m);
		for(int i=0; i<m; i++){
			vector newx = x.copy();
			newx[i] += dx;
			for(int j=0; j<n; j++){
				J[j,i] = (f(newx)[j] - fx[j])/dx;
			}
		}
		return J;
	}
	public static vector numderiv(Func<vector,double> f, vector x){
		double fx = f(x);
		int n = x.size;
		double dx = x.norm()*Pow(2.0,-26);
		if(dx == 0) dx = Pow(2.0,-26);
		vector J = new vector(n);
		for(int i=0; i<n; i++){
			vector newx = x.copy();
			newx[i] += dx;
			J[i] = (f(newx) - fx)/dx;
		}
		return J;
	}
	public static void timesJ(matrix A, int p, int q, double theta){
		double c = Cos(theta), s = Sin(theta);
		for(int i=0 ; i<A.size1 ; i++){
			double aip = A[i,p], aiq = A[i,q];
			A[i,p]=c*aip-s*aiq;
			A[i,q]=s*aip+c*aiq;
		}
	}
	public static void Jtimes(matrix A, int p, int q, double theta){
		double c=Cos(theta),s=Sin(theta);
		for(int j=0;j<A.size1;j++){
			double apj=A[p,j],aqj=A[q,j];
			A[p,j]= c*apj+s*aqj;
			A[q,j]=-s*apj+c*aqj;
		}
	}
	public static (vector,matrix) cyclic(matrix M){
		matrix A=M.copy();
		int n = M.size1;
		matrix V=matrix.id(n);
		vector w=new vector(n);
		/* run Jacobi rotations on A and update V */
		bool changed;
		do{
			changed=false;
			for(int p=0;p<n-1;p++){
				for(int q=p+1;q<n;q++){
					double apq=A[p,q], app=A[p,p], aqq=A[q,q];
					double theta=0.5*Atan2(2*apq,aqq-app);
					double c=Cos(theta),s=Sin(theta);
					double new_app=c*c*app-2*s*c*apq+s*s*aqq;
					double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
					if(new_app!=app || new_aqq!=aqq) // do rotation
						{
						changed=true;
						timesJ(A,p,q, theta); // A←A*J 
						Jtimes(A,p,q,-theta); // A←JT*A 
						timesJ(V,p,q, theta); // V←V*J
						}
				}
			}
		}while(changed);
		/* copy diagonal elements into w */
		for(int i=0 ; i<A.size1 ; i++){
			w[i] = A[i,i];
		}
		return (w,V);
	}

}

}
