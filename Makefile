#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = $(HOME)/repos/ppnm/matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')


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
