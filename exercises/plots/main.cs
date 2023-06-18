using static System.Console;

class main{
static void Main(string[] args){
	foreach(string arg in args){
		if(arg=="-gamma"){
			for(double x=-5.0 ; x<=5.0 ; x+=1.0/128){
				WriteLine($"{x} {sfuns.gamma(x)}");

			}
		}
		if(arg=="-erf"){
			for(double x=-3.5; x<=3.5; x+= 1.0/128){
				WriteLine($"{x} {sfuns.erf(x)}");
			}
		}
	}
}

}
