# Find subdirectories and $(MAKE) each of them.
subdirs=$(filter-out Makefile%, $(wildcard *))
all :
	for entry in $(subdirs); do echo $$entry && cd $$entry && $(MAKE); done