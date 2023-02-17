using System;
using static System.Console;
using static System.Math;

public static class main{
	public static void print(this double x,string s){ /* x.print("x="); */
		Write(s);WriteLine(x);
		}

	public static void print(this double x){ /* x.print() */
		x.print("");
	}


	public static void Main(){
		vec u = new vec(1,2,3);
		u.print("u = ");
		vec v = new vec(2,3,4);
		v.print("v = ");
		vec w = 2*v;
		w.print("2*v = ");
		double z = u.dot(v);
		z.print("u.dot(v) = ");
		z = vec.dot(u,v);
		z.print("vec.dot(u,v) = ");
		w = u.cross(v);
		w.print("u.cross(v) = ");
		w = vec.cross(u,v);
		w.print("vec.cross(u,v) = ");
		w = -v;
		w.print("-v = ");
	}
}
