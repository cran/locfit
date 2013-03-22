SHELL=/bin/sh
LIBNAME=locfit
LD=ld

OBJS=   S_enter.o adap.o density.o family.o frend.o kappa0.o kdtree.o lfstr.o linalg.o locfit.o nbhd.o random.o resid.o weight.o simul.o


locfit.so: $(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so

realclean:
	@rm -f Makefile *.o *.so
