Out.txt: main.exe
	mono main.exe > Out.txt
main.exe : cspractice.cs
	mcs -target:exe -out:main.exe cspractice.cs
clean:
	rm --force Out.txt main.exe


