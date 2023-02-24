#RM = rm ---force

Out.txt: main.exe Makefile
	mono main.exe > Out.txt

main.exe : main.cs # + any library.dll files
	mcs -target:exe -out:$@ $(filter %.cs,$^) $(addprefix -reference:,$(filter %.dll,$^))

.PHONEY : clean test

clean:
	$(RM) *.exe *.dll [Oo]ut* log*
test: 
	echo $(RM)


#How to make library
#
#library.dll : library.cs
#	mcs -target:library -out:$@ $(filter %.cs,$^)
