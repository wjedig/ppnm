Out.txt: main.exe
	mono main.exe > Out.txt
sfuncs.dll : sfuncs.cs
	mcs -target:library -out:sfuncs.dll sfuncs.cs
main.exe : main.cs sfuncs.dll
	mcs -target:exe -out:main.exe -reference:sfuncs.dll main.cs
