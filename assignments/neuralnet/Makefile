#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = ../../matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')
CALCULUS := $(shell find $(MATLIB)/calculus -name '*.cs')
GENLIST := $(shell find $(MATLIB)/genlist -name '*.cs')

all: Out.A.txt ANN_A.svg Out.B.txt ANN_B.svg

Out.A.txt out.afit1.data out.afit2 out.atrain.data: main.exe Makefile
	mono $< -A > $@ 

ANN_A.svg: ann_A.gpi out.afit1.data out.afit2.data out.atrain.data
	gnuplot $<

Out.B.txt out.bfit1.data out.bfit2.data out.btrain1.data out.btrain2.data: main.exe Makefile
	mono $< -B -acc:1e-3 > $@

ANN_B.svg: ann_B.gpi out.bfit1.data out.btrain1.data out.bfit2.data out.btrain2.data
	gnuplot $<

main.exe : main.cs neuralnet.cs linalg.dll calculus.dll
	$(MKEXE)

linalg.dll : $(LINALG) $(GENLIST)
	$(MKLIB)

calculus.dll : linalg.dll $(CALCULUS)
	$(MKLIB)

.PHONY : all clean test

clean:
	$(RM) *.exe *.dll [Oo]ut* log* *.svg
test: 
	echo $(RM)
