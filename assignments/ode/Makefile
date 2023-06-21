#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = ../../matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')
GENLIST := $(shell find $(MATLIB)/genlist -name '*.cs')

all: ode_test.svg ode_damp.svg ode_planets.svg

ode_test.svg : ode.gpi out.ode1.data out.ode2.data
	gnuplot ode.gpi

ode_damp.svg : damp.gpi out.ode3.data
	gnuplot damp.gpi

ode_planets.svg : planets.gpi BData
	gnuplot planets.gpi

out.ode1.data: main.exe
	mono main.exe -A1 > $@

out.ode2.data: main.exe
	mono main.exe -A2 > $@

out.ode3.data: main.exe
	mono main.exe -A3 > $@

BData: main.exe
	mono main.exe -B

main.exe : main.cs matlib.dll# + any library.dll files
	mcs -target:exe -out:$@ $(CODE) $(DLLS)

matlib.dll : $(LINALG) $(GENLIST)
	$(MKLIB)

.PHONEY : clean test all BData

clean:
	$(RM) *.exe *.dll [Oo]ut* log* *.svg
test: 
	echo $(RM)


#How to make library
#
#library.dll : library.cs
#	mcs -target:library -out:$@ $(filter %.cs,$^)
