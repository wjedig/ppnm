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
	public static matrix build_hamiltonian(double rmax, double dr){
		int npoints = (int)(rmax/dr)-1;
		vector r = new vector(npoints);
		for(int i=0 ; i<npoints ; i++){r[i]=dr*(i+1);}
		matrix H = new matrix(npoints,npoints);
		for(int i=0 ; i<npoints-1 ; i++){
			H[i,i] = -2;
			H[i,i+1] = 1;
			H[i+1,i] = 1;
		}
		H[npoints-1,npoints-1] = -2;
		matrix.scale(H,-(0.5/dr/dr));
		for(int i=0;i<npoints;i++){H[i,i]+=-1/r[i];}
		return H;
	}
}

public static class Test{
	public static void testA(int n=10, double scale = 1.0, bool trace=false){
		WriteLine("---Generating random symmetric matrix A---");
		matrix A = Utils.make_sym_matrix(n,scale);
		A.print("A given by:");
		WriteLine("---Calculating matrices V and D---");
		(vector w, matrix V) = jacobi.cyclic(A);
		matrix D = new matrix(n,n);
		for(int i=0 ; i<n ; i++){
			D[i,i] = w[i];
		}
		V.print("V given by:");
		D.print("D given by:");
		WriteLine("---Checking the requirements from the assignment description---");
		matrix VTAV = (V.T)*A*V;
		WriteLine($"VTAV == D : {VTAV.approx(D,acc:1e-3,eps:0)}");
		if(trace){VTAV.print("VTAV = ",format:"{0,20:g15}"); D.print("D = ",format:"{0,20:g15}");}
		matrix VDVT = V*D*(V.T);
		WriteLine($"VDVT == A : {VDVT.approx(A)}");
		matrix VTV = (V.T)*V;
		matrix ID = matrix.id(n);
		WriteLine($"VTV = I : {VTV.approx(ID)}");
		matrix VVT = V*(V.T);
		WriteLine($"VVT = I : {VVT.approx(ID)}");
	}
	public static void testB_eigenvals(double rmax, double dr){
		matrix ham = Utils.build_hamiltonian(rmax, dr);
		(vector e, matrix V) = jacobi.cyclic(ham);
		e.print("Eigenvalues given by:");
		V.print("With corresponding Eigenvectors given by:");
	}
	public static void testB_plot(double rmax, double dr){
		matrix ham = Utils.build_hamiltonian(rmax, dr);
		var (e,V) = jacobi.cyclic(ham);
		for(int i=0 ; i<V.size2 ; i++){
			Write($"{i*dr} ");
			for(int j=0 ; j<V.size1 ; j++){
				Write($"{V[j][i]/Sqrt(dr)} ");
			}
			Write("\n");
		}
	}
	public static void testB_rmax(double rmax, double dr){
		matrix ham = Utils.build_hamiltonian(rmax, dr);
		(vector e, _) = jacobi.cyclic(ham);
		WriteLine($"{rmax} {e[0]}");
	}
	public static void testB_dr(double rmax, double dr){
		matrix ham = Utils.build_hamiltonian(rmax, dr);
		(vector e, _) = jacobi.cyclic(ham);
		WriteLine($"{dr} {e[0]}");
	}
}

class main{
	public static void Main(string[] args){
		double rmax = 10;
		double dr = 0.3;
		int size = 10;
		double scale = 1.0;
		bool trace = false;
		foreach(string param in args){ // set all params
			string[] words = param.Split(":");
			if(words[0] == "-rmax"){rmax = double.Parse(words[1]);}
			if(words[0] == "-dr"){dr = double.Parse(words[1]);}
			if(words[0] == "-size"){size = int.Parse(words[1]);}
			if(words[0] == "-scale"){scale = double.Parse(words[1]);}
			if(words[0] == "-trace"){trace = true;}
		}
		foreach(string run in args){
			if(run == "-test"){Test.testA(size,scale,trace);}
			if(run == "-hydrogen"){Test.testB_eigenvals(rmax, dr);}
			if(run == "-plotrmax"){Test.testB_rmax(rmax, dr);}
			if(run == "-plotdr"){Test.testB_dr(rmax, dr);}
			if(run == "-plotfuncs"){Test.testB_plot(rmax=40, dr=0.08);}


		}
	}
}
