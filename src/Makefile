all: ../inst/bin/hawkes_bin

../inst/bin/hawkes_bin: hawkes/hawkes
	mkdir -p ../inst/bin/ \
	&& mv hawkes/hawkes ../inst/bin/hawkes_bin

hawkes/hawkes: Makefile_hawkes hawkes/Makefile
	cp Makefile_hawkes hawkes/Makefile \
	&& cd hawkes \
	&& make
