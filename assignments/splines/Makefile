#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = ../../matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')

all: Out.A.txt Linspline.svg Out.B.txt

Out.A.txt: main.exe Makefile
	mono main.exe -A < input.A.txt > $@

Out.plot.data: main.exe Makefile
	mono main.exe -plot < input.A.txt > $@

Out.B.txt: main.exe Makefile
	mono main.exe -B > $@

Linspline.svg: linspline.gpi Out.plot.data
	gnuplot $<

main.exe : main.cs matlib.dll
	$(MKEXE)

matlib.dll : $(LINALG)
	$(MKLIB)

.PHONY : all clean test

clean:
	$(RM) *.exe *.dll [Oo]ut* log* *.svg
test: 
	echo $(RM)


#How to make library
#
#library.dll : library.cs
#	mcs -target:library -out:$@ $(filter %.cs,$^)
