#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = $(HOME)/repos/ppnm/matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')

all: Out.A.txt Out.B.txt Error.svg Quasi.svg


Out.A.txt: main.exe Makefile
	mono main.exe -A > $@

Out.B.txt: main.exe Makefile
	mono main.exe -B > $@

Error.svg: error.gpi out.a.data
	gnuplot $<

Quasi.svg: quasi.gpi out.b.data
	gnuplot $<

main.exe : main.cs matlib.dll # + any library.dll files
	$(MKEXE)

matlib.dll : $(LINALG)
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
