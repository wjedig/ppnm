using System;
using static System.Console;
using static System.Math;
using linalg;

public static class cholesky{
	public static matrix decomp(matrix A){
		int n = A.size1; // obtain dimension of A
		if(A.size2!=n) throw new ArgumentException("Matrix must be symmetric");
		matrix L = new matrix(n,n); // same dimensions as A
		for(int i=0; i<n; i++){ // iteration over rows
			for(int j=0; j<=i; j++){ // iterate over columns up to and including diagonal elements
				double sum = 0.0;  // this term is required both for diagonal and non diagonal elements, so we calculate now
				for(int k=0; k<j; k++){
					sum += L[i,k] * L[j,k];
				}
				if(i==j) L[i,j] = Sqrt(A[i,i] - sum); 
				else L[i,j] = (1.0/L[j,j]) * (A[i,j] - sum);
			}
		}
		return L;
	}
	public static vector solve(matrix A, vector b, matrix L = null, matrix LT = null){
		int n = b.size;
		if(A.size1!=n) throw new ArgumentException("vector b must have same dimension as the rows in the matrix");
		if(L==null || LT==null){
			L = decomp(A);
			LT = L.T;
		}
		vector y = new vector(n); // I'm not solving in-place, as I find it less intuitive and explanatory
		for(int i=0; i<n; i++){ // Forward substitution solving Ly=b
			double sum = 0.0;
			for(int k=0; k<i; k++){
				sum += L[i,k]*y[k];
			}
			y[i] = (1.0/L[i,i]) * (b[i] - sum);
		}
		vector x = new vector(n);
		for(int i=n-1; i>=0; i--){
			double sum = 0.0;
			for(int k=i+1; k<n; k++){
				sum += LT[i,k]*x[k];
			}
			x[i] = (1.0/LT[i,i]) * (y[i] - sum);
		}
		return x;
	}
	public static double determinant(matrix A){
		matrix L = decomp(A);
		double det = 1.0;
		for(int i=0; i<A.size1; i++){
			det *= L[i,i]; // for triangular matrix, the determinant is the product of diagonals
		}
		return det*det; // the diagonals on L and LT are identical, thus det(L)det(LT) = (det(L))^2
	}
	public static matrix inverse(matrix A){
		matrix L = decomp(A);
		matrix LT = L.T;
		int n = A.size1;
		matrix invA = new matrix(n,n);
		vector ei = new vector(n);
		vector xi = new vector(n);
		for(int i=0; i<n; i++){
			ei[i] = 1;
			if(i!=0) ei[i-1]=0;
			xi = solve(A,ei,L,LT);
			invA[i] = xi;
		}
		return invA;
	}
}
