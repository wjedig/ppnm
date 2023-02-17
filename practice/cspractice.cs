using System;
using static System.Console;

class Test{
	static void Main(){
		double s = 0;
		for(double i=1;i<1000;i++){
			s += 1.0/i;
		}
		WriteLine(s);
	}
}	
