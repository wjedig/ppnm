using System;
using static System.Console;
using static System.Math;

class main{
	public static void Main(){
		tiny();
		fpre();
	}
	public static void max_e(){	
		int i = 1;
		while(i+1>i){i++;}
		WriteLine($"my max int = {i}");
		WriteLine($"Tabulated max int: {int.MaxValue}");
	}
	public static void min_e(){
		int i = 1;
		while(i-1<i){i--;}
		WriteLine($"my min int = {i}");
		WriteLine($"Tabulated min int: {int.MinValue}");
	}
	public static void machine_epsilon(){
		double x=1; while(1+x!=1){x/=2;} x*=2;
		WriteLine($"double machine epsilon: {x}");
		WriteLine($"vs theoretical: {Pow(2,-52)}");	
		float y=1F; while((float)(1F+y) != 1F){y/=2F;} y*=2F;
		WriteLine($"float machine epsilon: {y}");	
		WriteLine($"vs theoretical: {Pow(2,-23)}");	
	}
	public static void tiny(){
		int n=(int)1e6;
		double epsilon=Pow(2,-52);
		double tiny=epsilon/2;
		double sumA=0,sumB=0;
		sumA+=1; for(int i=0;i<n;i++){sumA+=tiny;}
		for(int i=0;i<n;i++){sumB+=tiny;} sumB+=1;
		WriteLine($"sumA-1 = {sumA-1:e} should be {n*tiny:e}");
		WriteLine($"sumB-1 = {sumB-1:e} should be {n*tiny:e}");
	}
	public static void fpre(){
		double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
		double d2 = 8*0.1;

		WriteLine($"d1={d1:e15}");
		WriteLine($"d2={d2:e15}");
		WriteLine($"d1==d2 ? => {d1==d2}");

		WriteLine($"Are they approximately equal? => {approx(d1,d2)}");
	}
	public static bool approx(double a, double b, double acc=1e-9, double eps=1e-9){
		if(Abs(b-a) < acc) return true;
		else if(Abs(b-a) < Max(Abs(a),Abs(b))*eps) return true;
		else return false;
	}
}
