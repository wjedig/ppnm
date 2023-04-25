using System;
using static System.Console;
using static System.Math;
using linalg;

namespace linalg{
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
}//QRQS
}//namespace
