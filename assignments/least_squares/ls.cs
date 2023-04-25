using System;
using linalg;

namespace linalg{
public static class LS{
	public static (vector, matrix) lsfit(Func<double,double>[] fs, vector x, vector y, vector dy){
		int n = x.size;
		int m = fs.Length;
		matrix A = new matrix(n,m);
		vector b = new vector(n);
		for(int i=0 ; i<n ; i++){
			b[i] = y[i]/dy[i];
			for(int j=0 ; j<m ; j++){
				A[i,j] = fs[j](x[i])/dy[i];
			}
		}
		(matrix Q, matrix R) = QRGS.decomp(A);
		vector c = QRGS.solve(Q,R,b);
		matrix ATA = (A.T)*A;
		matrix S = QRGS.inverse(ATA);
		return (c,S);
	}
}//class
}//namespace
