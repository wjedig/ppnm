CSC = mcs
RUN = mono

all: out.txt; cat $<

out.txt: main.exe
	$(RUN) $< > $@

main.exe: main.cs matrix.dll
	$(CSC) -reference:$(word 2,$^) $<

matrix.dll : vector.cs matrix.cs
	$(CSC) -target:library -out:$@ $^

vector.dll: vector.cs
	$(CSC) -target:library -out:$@ $^

clean:
	$(RM) *.dll *.exe [Oo]ut*
