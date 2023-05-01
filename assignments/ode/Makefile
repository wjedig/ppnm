#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = $(HOME)/repos/ppnm/matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')
GENLIST := $(shell find $(MATLIB)/genlist -name '*.cs')

all: ode_test.svg ode_damp.svg

ode_test.svg : ode.gpi Out.ode1.data Out.ode2.data
	gnuplot ode.gpi

ode_damp.svg : damp.gpi Out.ode3.data
	gnuplot damp.gpi

Out.ode1.data: main.exe
	mono main.exe -A1 > $@

Out.ode2.data: main.exe
	mono main.exe -A2 > $@

Out.ode3.data: main.exe
	mono main.exe -A3 > $@

main.exe : main.cs matlib.dll# + any library.dll files
	mcs -target:exe -out:$@ $(CODE) $(DLLS)

matlib.dll : $(LINALG) $(GENLIST)
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
