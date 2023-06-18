using System;
using static System.Console;
using static System.Math;
using linalg;

public static class utils{
	public static matrix makePosDefMatrix(int size){
		var rand = new System.Random();
		matrix A = new matrix(size,size);
		for(int i=0; i<size; i++){
			for(int j=0; j<size; j++){
				A[i,j] = rand.NextDouble(); 
			}
		}
		matrix AAT = A*(A.T);
		return AAT;
	}
	public static vector makeVector(int size){
		var rand = new System.Random(1);
		vector b = new vector(size);
		for(int i=0; i<size; i++) b[i] = rand.NextDouble();
		return b;
	}
}

public static class tasks{
	public static void showDecomp(int size){
		matrix A = utils.makePosDefMatrix(size);
		matrix L = cholesky.decomp(A);
		WriteLine("---Testing Cholesky-Banachiewicz decomposition algorithm---");
		A.print("A given by:");
		L.print("L is calculated to be:");
		matrix LLT = L*(L.T);
		WriteLine("This can be verified by checking that LL^T is equal to A");
		WriteLine($"LL^T = A: {A.approx(LLT)}");
	}	
	public static void showLinEqSolve(int size){
		WriteLine("---Testing linear equation solver using Cholesky decomposition---");
		matrix A = utils.makePosDefMatrix(size);
		vector b = utils.makeVector(size);
		A.print("A = ");
		b.print("b = ");
		vector x = cholesky.solve(A,b);
		x.print("Solution given by: ");
		WriteLine("This can be verified by checking that A*x is equal to b -");
		WriteLine($"Ax=b: {b.approx(A*x)}");
	}
	public static void showDeterminant(int size){
		WriteLine("---Testing determinant solver using Cholesky decomposition---");
		matrix A = utils.makePosDefMatrix(size);
		A.print("A = ");
		double det = cholesky.determinant(A);
		WriteLine("Determinant of A found to be:");
		WriteLine($"Det(A) = {det}");
		WriteLine("It is not easy to verify whether this result is correct, but we know it must be positive:");
		WriteLine($"Det(A) > 0: {det>0.0}");
		WriteLine("I also verified the results using an online determinant calculator for 2 test matrices.");
	}
	public static void showInverse(int size){
		WriteLine("---Testing inverse solver using Cholesky decomposition---");
		matrix A = utils.makePosDefMatrix(size);
		A.print("A = ");
		matrix invA = cholesky.inverse(A);
		invA.print("A^-1 found to be: ");
		WriteLine("This can be verified by checking AA^-1 = I");
		matrix I = matrix.id(size);
		WriteLine($"AA^-1 = I: {I.approx(A*invA)}");
	}
}


class main{
	public static int size = 10;
	public static void Main(string[] args){
		foreach(string arg in args){
			string[] words = arg.Split(":");
			if(words[0]=="-size") size = int.Parse(words[1]);
		}
		foreach(string run in args){
			if(run=="-decomp") tasks.showDecomp(size);
			if(run=="-linSolve") tasks.showLinEqSolve(size);
			if(run=="-determinant") tasks.showDeterminant(size);
			if(run=="-inverse") tasks.showInverse(size);
		}
	}
}
