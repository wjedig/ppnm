#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = ../../matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')

all : Out.txt Decay.svg

Out.txt: main.exe
	mono main.exe -testB -input:decay.txt >$@

#out.lndecay.data : main.exe
#	mono main.exe -input:decay.txt -writeData > $@

out.decayfit.par: main.exe
	mono main.exe -input:decay.txt -writeFit > $@

Decay.svg : decay.gpi out.decayfit.par
	gnuplot $<

main.exe : main.cs matlib.dll
	$(MKEXE)

matlib.dll : ls.cs $(LINALG) 
	$(MKLIB)

.PHONEY : clean test

clean:
	$(RM) *.exe *.dll [Oo]ut* log* err start.par *.data *.svg
test: 
	echo $(RM)




