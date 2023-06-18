DIR = ../../matlib/complex

Out.txt: main.exe
	mono $< >$@

main.exe : main.cs cmath.dll
	mcs -reference:cmath.dll main.cs

cmath.dll : $(DIR)/cmath.cs $(DIR)/complex.cs
	mcs -target:library -out:./cmath.dll $^

test:
	echo  $(DIR)/cmath.cs

.PHONEY : clean test


clean:
	$(RM) *.exe *.dll [Oo]ut* log*
