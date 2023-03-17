using System;
using static System.Math;
using static System.Console;

class main{
	public static void Main(string[] args){
		double start = 0;
		double stop = 1000;
		foreach(string arg in args){
			string[] words = arg.Split(':');
			if(words[0] == "-start"){
				start = double.Parse(words[1]);
			}
			if(words[0] == "-stop"){
				stop = double.Parse(words[1]);
			}
		}
		
		double spacing = (stop-start)/1000.0;
		for(double i=start ; i<=stop ; i+=spacing){
			WriteLine($"{i} {0.4 + (Pow(i,3)/160000000)}");
		}
	}
}
