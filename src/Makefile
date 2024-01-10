all: ../inst/bin/hawkes_bin

../inst/bin/hawkes_bin: hawkes/hawkes
	mkdir -p ../exec/
	mv hawkes/hawkes ../exec/hawkes_bin

hawkes/hawkes: hawkes/Makefile
	cd hawkes \
	&& make
