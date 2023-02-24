DIR = $(HOME)/repos/ppnm/matlib/complex

cmath.dll : $(DIR)/cmath.cs $(DIR)/complex.cs
	mcs -target:library -out:./cmath.dll $^

#main.exe : main.cs cmath.dll
#	mcs -reference:cmath.dll main.cs

test:
	echo  $(DIR)/cmath.cs

.PHONEY : clean test


clean:
	$(RM) *.exe *.dll [Oo]ut* log*
