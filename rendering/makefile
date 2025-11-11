CC=gcc
BASEFLAGS= -Wall -lm
CFLAGS= $(BASEFLAGS) -O3 -std=gnu99
CDBGFLAGS= $(BASEFLAGS) -g
DEPS=linear_algebra.h hgmt_lor_creator.h helper_functions.h read_write.h imager.h hgmt_structs.h

hgmt_debug: imager.o ../src/hgmt_structs.o ../src/linear_algebra.o ../src/helper_functions.o ../src/read_write.o ../src/pore_geometry.o
	$(CC) -o imager $^  $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
