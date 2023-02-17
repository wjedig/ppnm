#RM = rm ---force

Out.txt: main.exe Makefile
	mono main.exe > Out.txt

main.exe : main.cs input.dll
	mcs -target:exe -out:$@ $(filter %.cs,$^) $(addprefix -reference:,$(filter %.dll,$^))

input.dll : input.cs
	mcs -target:library -out:$@ $(filter %.cs,$^)

numbers.txt : Makefile
	echo -"numbers:1,2,3,4,5,6,7,8,9" > $@

.PHONEY : clean test

clean:
	$(RM) *.exe *.dll [Oo]ut* log*
test: 
	echo $(RM)
