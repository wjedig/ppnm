

all: Gamma.svg Erf.svg

Gamma.svg: gamma.gpi gamma.data factorials.data
	gnuplot $<

Erf.svg: erf.gpi tabulated_erf.data erf.data
	gnuplot $<

gamma.data: main.exe
	mono $< -gamma >$@

erf.data: main.exe
	mono $< -erf >$@

main.exe : main.cs sfuns.dll
	mcs -target:exe -out:$@ \
		$(filter %.cs, $^) $(addprefix -reference:,$(filter %.dll,$^))

sfuns.dll : gamma.cs erf.cs
	mcs -target:library -out:$@ $^

factorials.data : Makefile
	>$@
	echo 0 1 >> $@
	echo 1 1 >> $@
	echo 2 2 >> $@
	echo 3 6 >> $@
