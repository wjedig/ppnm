
#RM = rm ---force

MATLIB = ../../matlib

all: main.exe Out.txt out.times.data Opcount.svg

Opcount.svg: opcount.gpi out.times.data
	gnuplot $<

Out.txt: main.exe Makefile
	>$@
	mono main.exe -test 1>>$@

main.exe : main.cs matlib.dll
	mcs -target:exe -out:$@ $(filter %.cs,$^) $(addprefix -reference:,$(filter %.dll,$^))

matlib.dll : $(MATLIB)/linalg/vector.cs  $(MATLIB)/linalg/matrix.cs $(MATLIB)/complex/complex.cs
	mcs -target:library -out:$@ $(filter %.cs,$^)

#fit.data : fit.exe
#	>$@
#	mono $< 1>>$@

#fit.exe : fit.cs
#	mcs -target:exe -out:$@ $<

out.times.data : main.exe
	>$@
	for N in $$(seq 100 100 800); do \
		time --format "$$N %e" --output $@ --append \
		mono $< -size:$$N 2>err ;\
	done

.PHONEY : clean test

clean:
	$(RM) *.exe *.dll [Oo]ut* log* err* *.data *.log *.svg
test: 
	echo $(RM)


#How to make library
#
#library.dll : library.cs
#	mcs -target:library -out:$@ $(filter %.cs,$^)
