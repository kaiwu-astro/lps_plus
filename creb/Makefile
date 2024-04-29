export OPENGL=0

# ifeq ($(shell $(CC) -v 2>&1 | grep -c "clang"), 1)
# export OPENMPCLANG=1
# else
# export OPENMP=1
# endif

REB_DIR=/path/to/rebound
PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)

include $(REB_DIR)/src/Makefile.defs

all: librebound.so
	@echo ""
	@echo "Compiling creb.c file using "$CC
	$(CC) -I$(REB_DIR)/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) creb.c -L. -lrebound $(LIB) -mcmodel=large -O3 -o creb
	@echo ""
	@echo "creb.c file compiled successfully. Name of executable: creb . Run ./creb for usage, or see https://github.com/kaiwu-astro/lps_plus"

librebound.so:
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(REB_DIR)/src/
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
	@ln -s $(REB_DIR)/src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C $(REB_DIR)/src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
