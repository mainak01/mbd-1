all: libmbd.a

libmbd.a: $(FOBJS)
	@$(AR) -cr $@ $^
	@rm -f $(wildcard lib*/_mbd_backend*)

$(FOBJS): fcompile ;

fcompile: | fcompile.py
	@python $| <fcompile.json

clean:
	rm -f $(wildcard *.mod *.o *.a _fcompile_cache.json)

fcompile.py:
	curl -k https://raw.githubusercontent.com/azag0/fcompile/07b7133500b4953b9f90d0e888f5908ce15ae316/$@ >$@
