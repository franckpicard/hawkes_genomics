all: ../inst/bin/hawkes_bin

../inst/bin/hawkes_bin: hawkes/hawkes
	pwd
	ls -R ..
	cp hawkes/hawkes ../inst/bin/hawkes_bin

hawkes/hawkes: hawkes/Makefile
	cd hawkes \
	&& make
