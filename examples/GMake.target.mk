.SUFFIXES:

include ../GMake.defs.mk

ODIR := _build.$(suf)
$(warning Using ODIR=$(ODIR))

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/GNUmakefile \
	     SRCDIR=$(CURDIR) $(MAKECMDGOALS)

.PHONY: $(ODIR)
$(ODIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

GNUmakefile : ;
%.mk :: ;

% :: $(ODIR) ; :

$(warning Using Goal: $(MAKECMDGOALS))
$(warning Using Target: $(MAKETARGET))

.PHONY: realclean
realclean:
	rm -rf $(ODIR)
	rm -f *.exe

.PHONY: clean
clean:
	rm $(ODIR)/*.o $(ODIR)/*.mod
