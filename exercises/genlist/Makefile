Out.txt: main.exe
	mono main.exe > Out.txt

main.exe: main.cs genlist.dll
	mcs -target:exe -out:$@ \
		$(addprefix -reference:,$(filter %.dll,$^)) \
		$(filter %.cs,$^)
		
genlist.dll : genlist.cs
	mcs -target:library -out:$@ $<
