using System;
using static System.Console;
using static System.Math;
using System.Diagnostics;
using linalg;

public static class QRGS{
	public static (matrix, matrix) decomp(matrix A){ // a: (nxm)
		matrix Q = A.copy();
		int m = A.size2;
		matrix R = new matrix(m,m);
		for(int i=0; i<m; i++){
			R[i,i] = Q[i].norm();
			Q[i]/= R[i,i];
			for(int j=i+1; j<m; j++){
				R[i,j] = Q[i].dot(Q[j]);
				Q[j] -= Q[i]*R[i,j];
			}
		}
		
		return (Q,R);
	}
	public static vector solve(matrix Q, matrix R, vector b){
		vector c = (Q.T)*b;
		for(int i=c.size-1; i>=0; i--){
			double sum = 0;
			for(int k=i+1; k<c.size; k++){
				sum += R[i,k]*c[k];
			}
			c[i] = (c[i]-sum)/R[i,i];
		}
		return c;
	}
	public static matrix inverse(matrix A, matrix Q = null, matrix R = null){
		//(int n, int m) = (A.size1, A.size2);
		//Debug.Assert(n==m,string.Format("Only square matrices are invertible. A has dimension ({0},{1})",n,m));
		//WriteLine($"n = {n}, m = {m}");

		int n = A.size1;
		matrix invA = new matrix(n,n);
		vector ei = new vector(n);
		ei[0] = 1;
		
		if((Q == null) || (R == null)){
			(Q, R) = QRGS.decomp(A);
		}
		vector xi = QRGS.solve(Q,R,ei);
		invA[0] = xi;
		for(int i=1; i<n; i++){
			ei[i-1] = 0;
			ei[i] = 1;
			xi = QRGS.solve(Q,R,ei);
			invA[i] = xi;
		}
		return invA;
	}


}

public static class Test{
	public static void test(){	
		vector a = new vector(1,2,3);
		a.dot(a);
		
		// Test QR Decomp
		var rnd = new System.Random(1);
		(int n, int m) = (30,10);
		matrix A = new matrix(n,m);
		for(int i=0; i<n; i++){
			for(int j=0; j<m; j++){
				A[i,j] = rnd.NextDouble()*100.0; // Fill matrix with numbers between 0-100
			}
		}
		(matrix Q, matrix R) = QRGS.decomp(A);
		R.print("Check that R is upper triangular, R = \n");
		matrix QTQ = (Q.T)*Q;
		//QTQ.print("\n Check that QTQ is the identity matrix: \n");
		Write($"\n Check that QTQ is the identity matrix by using Approx method: \n QTQ = 1: {QTQ.approx(matrix.id(m))}\n");
		WriteLine("Check that QR = A");
		WriteLine($"QR = A: {A.approx(Q*R)}");
		
		// Test Solve:
		
		matrix A2 = new matrix(m,m);
		vector b2 = new vector(m);
		for(int i=0; i<m; i++){
			b2[i] = rnd.NextDouble()*100.0;
			for(int j=0; j<m; j++){
				A2[i,j] = rnd.NextDouble()*100.0;
			}
		}
		
		(matrix Q2, matrix R2) = QRGS.decomp(A2);
		vector x2 = QRGS.solve(Q2, R2, b2);
		
		WriteLine("Check that Ax = b");
		WriteLine($"Ax = b : {b2.approx(A2*x2)}");

		// Test inverse:
		
		//matrix Id = matrix.id(10);
		WriteLine("************************************");
		matrix A3 = QRGS.inverse(A2);
		WriteLine("Check that A*A^-1 = I");
		matrix A4 = A2*A3;
		A4.print("A*A^-1 = ");
		WriteLine($"A*A^-1 = I: {A4.approx(matrix.id(A2.size1))}");
	}
}

class main{
	public static void Main(string[] args){
		//WriteLine("running main i guess uwu");
		foreach(string arg in args){
			//WriteLine(arg);
			string[] words=arg.Split(':');
			if(words[0]=="-size"){
				int n = int.Parse(words[1]);
				WriteLine($"n = {n}");
				matrix A = new matrix(n,n);
				var rnd = new System.Random(1);
				for(int i=0; i<n; i++){
					for(int j=0; j<n; j++){
						A[i,j] = rnd.NextDouble()*100.0; // Fill matrix with floats 0-100
					}
				}
				(matrix Q, matrix R) = QRGS.decomp(A);
			}
			if(words[0]=="-test"){
				Test.test();
			}
		}//parse args
	}//Main
}//main
