#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = $(HOME)/repos/ppnm/matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')
CALCULUS := $(shell find $(MATLIB)/calculus -name '*.cs')
GENLIST := $(shell find $(MATLIB)/genlist -name '*.cs')

all: Out.A.txt ANN.svg

Out.A.txt out.afit.data out.atrain.data: main.exe Makefile
	mono $< -A  # > $@ 

ANN.svg: ann.gpi out.afit.data out.atrain.data
	gnuplot $<

main.exe : main.cs neuralnet.cs linalg.dll calculus.dll
	$(MKEXE)

linalg.dll : $(LINALG) $(GENLIST)
	$(MKLIB)

calculus.dll : linalg.dll $(CALCULUS)
	$(MKLIB)

.PHONY : all clean test

clean:
	$(RM) *.exe *.dll [Oo]ut* log*
test: 
	echo $(RM)
