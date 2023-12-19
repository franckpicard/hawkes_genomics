all: hawkes_bin

hawkes_bin: hawkes/hawkes
	cp hawkes/hawkes hawkes_bin

hawkes/hawkes: hawkes/Makefile
	cd hawkes \
	&& make
