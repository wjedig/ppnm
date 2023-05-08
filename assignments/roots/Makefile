#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = $(HOME)/repos/ppnm/matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')
CALCULUS := $(shell find $(MATLIB)/calculus -name '*.cs')
GENLIST := $(shell find $(MATLIB)/genlist -name '*.cs')

all: Out.TestJacobi.txt Out.A.txt Convergence.svg

Out.TestJacobi.txt: main.exe Makefile
	mono $< -test > $@

Out.A.txt: main.exe Makefile
	mono $< -A > $@ 

Convergence.svg: convergence.gpi out.b.data
	gnuplot $<

out.b.data: main.exe Makefile
	mono $< -B:plot

main.exe : main.cs linalg.dll calculus.dll
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
