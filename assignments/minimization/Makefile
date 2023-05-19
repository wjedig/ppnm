#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = $(HOME)/repos/ppnm/matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')
GENLIST := $(shell find $(MATLIB)/genlist -name '*.cs')
CALCULUS := $(shell find $(MATLIB)/calculus -name '*.cs')

all: Out.A.txt Out.B.txt Higgs.svg

Out.A.txt: main.exe
	mono main.exe -A > $@

Out.B.txt out.bparams.par: main.exe Makefile
	mono main.exe -B < higgs.data > $@

Higgs.svg: higgs.gpi out.bparams.par higgs.data
	gnuplot $<

main.exe : main.cs linalg.dll calculus.dll
	$(MKEXE)

linalg.dll : $(LINALG) $(GENLIST)
	$(MKLIB)

calculus.dll : linalg.dll $(CALCULUS)
	$(MKLIB)

.PHONEY : clean test all

clean:
	$(RM) *.exe *.dll [Oo]ut* log*
test: 
	echo $(RM)


#How to make library
#
#library.dll : library.cs
#	mcs -target:library -out:$@ $(filter %.cs,$^)
