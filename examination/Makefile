#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = $(HOME)/repos/ppnm/matlib
LINALG := $(shell find $(MATLIB)/linalg -name '*.cs')
GENLIST := $(shell find $(MATLIB)/genlist -name '*.cs')
CALCULUS := $(shell find $(MATLIB)/calculus -name '*.cs')
space = echo "\n\n\n" >>$@

all: Out.txt #Out.Decomp.txt Out.LinSolve.txt Out.Determinant.txt Out.Inverse.txt

Out.txt: main.exe
	>$@
	>err
	>log
	mono main.exe -decomp 0>>err 1>>$@ 2>>log
	$(space)
	mono main.exe -linSolve 0>>err 1>>$@ 2>>log
	$(space)
	mono main.exe -determinant 0>>err 1>>$@ 2>>log
	$(space)
	mono main.exe -inverse 0>>err 1>>$@ 2>>log

main.exe : main.cs cholesky.cs linalg.dll calculus.dll
	$(MKEXE)

linalg.dll : $(LINALG) $(GENLIST)
	$(MKLIB)

calculus.dll : linalg.dll $(CALCULUS)
	$(MKLIB)

.PHONEY : clean test all

clean:
	$(RM) *.exe *.dll [Oo]ut* log* err*
test: 
	echo $(RM)
