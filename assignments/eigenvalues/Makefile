#RM = rm ---force
DLLS = $(addprefix -reference:,$(filter %.dll,$^))
CODE = $(filter %.cs,$^)
MKEXE = mcs -target:exe -out:$@ $(DLLS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(DLLS) $(CODE)
MATLIB = $(HOME)/repos/ppnm/matlib

all : Out.testA.txt dr.svg rmax.svg

Out.testA.txt: main.exe Makefile
	mono main.exe -test > $@

main.exe : main.cs matlib.dll
	$(MKEXE)

dr.svg: dr.gpi out.dr.data
	gnuplot $<

rmax.svg: rmax.gpi out.rmax.data
	gnuplot $<

eigenfuncs.svg: eigenfuncs.gpi out.eigenfuncs.data
	gnuplot $<

out.dr.data : main.exe
	>$@
	for N in $$(seq 0.05 0.05 3); do \
		mono $< -plotdr -rmax:10 -dr:$$N 1>>$@ 2>err ;\
	done

out.rmax.data : main.exe
	>$@
	for N in $$(seq 1 0.1 10); do \
		mono $< -plotrmax -rmax:$$N -dr:0.1 1>>$@ 2>err ;\
	done

out.eigenfuncs.data : main.exe
	>$@
	mono main.exe -plotfuncs >>$@

matlib.dll : evd.cs $(MATLIB)/linalg/vector.cs  $(MATLIB)/linalg/matrix.cs
	$(MKLIB)

.PHONEY : clean test

clean:
	$(RM) *.exe *.dll [Oo]ut* log* *.data err
test: 
	echo $(RM)


#How to make library
#
#library.dll : library.cs
#	mcs -target:library -out:$@ $(filter %.cs,$^)
