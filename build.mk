VPATH = $(MAINDIR):$(SRCDIR)

all: $(PYTHONMODULE)

$(PYTHONMODULE): build_mbd.py mbd.h libmbd.a
	@$(PYTHON) $^
	@$(MPIFC) $(LINKFLAG) -o $(PYTHONMODULE) _mbd.o libmbd.a
	@echo Created $(PYTHONMODULE).

libmbd.a: $(FOBJS)
	@ar -cr $@ $^

$(FOBJS): fcompile ;

fcompile: | fcompile.py
	@python $| <config.json

clean:
	rm -f $(wildcard *.mod *.o *.a *.so *.dylib *.c _fcompile_cache.json)

fcompile.py:
	curl https://raw.githubusercontent.com/azag0/fcompile/cbcb4fd2318b22368214e29557799a3711baf49f/$@ >$@
