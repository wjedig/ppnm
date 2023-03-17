using System;
using static System.Console;
using static System.Math;
using linalg;

public static class Utils{
	public static matrix make_sym_matrix(int n = 10, double scale = 100.0){
		matrix A = new matrix(n,n);
		var rnd = new System.Random(1);
		for(int i=0 ; i<n ; i++){
			for(int j=i ; j<n ; j++){
				double ij = rnd.NextDouble()*scale;
				A[i,j] = ij;
				if(i!=j){A[j,i] = ij;}
			}
			//A.print($"in the {i}th iteration, A = ");
		}
		return A;
	}
}

public static class Test{
	public static void testA(int n=10, double scale = 100.0, bool trace=false){
		matrix A = Utils.make_sym_matrix(n,scale);
		(vector w, matrix V) = jacobi.cyclic(A);
		matrix D = new matrix(n,n);
		for(int i=0 ; i<n ; i++){
			D[i,i] = w[i];
		}
		matrix VTAV = (V.T)*A*V;
		WriteLine($"VTAV == D : {VTAV.approx(D)}");
		if(trace){VTAV.print("VTAV = "); D.print("D = ");}
		matrix VDVT = V*D*(V.T);
		WriteLine($"VDVT == A : {VDVT.approx(A)}");
		matrix VTV = (V.T)*V;
		matrix ID = matrix.id(n);
		WriteLine($"VTV = I : {VTV.approx(ID)}");
		matrix VVT = V*(V.T);
		WriteLine($"VVT = I : {VVT.approx(ID)}");

	}
}

class main{
	public static void Main(string[] args){
		int rmax = 10;
		double dr = 0.3;
		int size = 10;
		double scale = 100.0;
		bool trace = false;
		foreach(string param in args){ // set all params
			string[] words = param.Split(":");
			if(words[0] == "-rmax"){rmax = int.Parse(words[1]);}
			if(words[0] == "-dr"){dr = double.Parse(words[1]);}
			if(words[0] == "-size"){size = int.Parse(words[1]);}
			if(words[0] == "-scale"){scale = double.Parse(words[1]);}
			if(words[0] == "-trace"){trace = true;}
		}
		foreach(string run in args){
			if(run == "-test"){Test.testA(size,scale,trace);}
			if(run == "-hydrogen"){Test.testA();}

		}
	}
}
