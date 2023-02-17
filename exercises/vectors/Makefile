Out.txt: main.exe
	mono $< > Out.txt

vec.dll : vec.cs
	mcs $< -target:library -out:$@

main.exe : main.cs vec.dll
	mcs -target:exe -out:$@ $(addprefix -reference:,$(filter %.dll,$^)) $(filter %.cs,$^)

clean:
	rm --force *.exe *.dll [Oo]ut* log*
