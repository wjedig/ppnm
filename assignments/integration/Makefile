#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = ../../matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')

all: Out.A.txt Out.B.txt Erf.svg

Out.A.txt: main.exe
	mono main.exe -A

Out.B.txt: main.exe
	mono main.exe -B >$@

Erf.svg: erf.gpi a.data
	gnuplot $<

main.exe : main.cs # + any library.dll files
	mcs -target:exe -out:$@ $(CODE) $(DLLS)

.PHONEY : clean test all

clean:
	$(RM) *.exe *.dll [Oo]ut* log*
test: 
	echo $(RM)


#How to make library
#
#library.dll : library.cs
#	mcs -target:library -out:$@ $(filter %.cs,$^)
