name = pqgen
cmds = pqgen pq-theta pq-codon2nonsyn pq-dna2codon
docs = pqgen.1 pq-theta.1 pq-dna2codon.1
headers = pq_genetics pq_sfobj pq_sfstats
libname = lib$(name)
info = codon2aa codon2syn

HOME = $(shell echo $$HOME)/
BASE = $(HOME).local/
INSTALL = $(HOME).local/

Wgcc = -Wall -Wextra -Wpedantic

###

bin = bin/
inc = include/
lib = lib/
obj = obj/
man = doc/
share = share/
src = src/

bins = $(addprefix $(bin),$(cmds))
incs = $(addsuffix .h,$(addprefix $(inc),$(headers)))
libs = $(lib)$(libname).so
mans = $(addprefix $(man),$(docs))
objs = $(addsuffix .o,$(addprefix $(obj),$(headers)))

######

.PHONY:	all
all:	$(bins) $(libs) $(objs)

$(bin)%:	$(src)%.c $(libs)
	mkdir -p $(bin)
	gcc -I include/ -I $(BASE)include/librawk/ -L $(BASE)lib/ -L $(lib) $(Wgcc) -o $@ $(word 1,$^) -lrawk -lm -lpqgen

$(libs):	$(objs)
	mkdir -p $(lib)
	gcc -shared -o $@ $^

$(obj)%.o:	$(src)%.c
	mkdir -p $(obj)
	gcc -I $(inc) -c $(Wgcc) -fpic -o $@ $^

###

.PHONY:	clean
clean:
	-rm $(objs) $(libs)

######

ibin = $(INSTALL)$(bin)
ilib = $(INSTALL)$(lib)
iman = $(INSTALL)share/man/man1/
iinc = $(INSTALL)$(inc)$(libname)/
ishare = $(INSTALL)share/$(libname)/

IBIN = $(addprefix $(ibin),$(cmds))
IINC = $(addsuffix .h,$(addprefix $(iinc),$(headers)))
ILIB = $(addprefix $(INSTALL),$(libs))
IMAN = $(addprefix $(iman),$(docs))
ISHARE = $(addprefix $(ishare),$(info))

.PHONY:	install
install:	$(IBIN) $(IINC) $(ILIB) $(IMAN) $(ISHARE)

$(ibin)%:	$(bin)%
	mkdir -p $(ibin)
	cp $^ $@

$(iinc)%.h:	$(inc)%.h
	mkdir -p $(iinc)
	cp $^ $@

$(ilib)%.so:	$(lib)%.so
	mkdir -p $(ilib)
	cp $^ $@

$(iman)%.1:	$(man)%.1
	mkdir -p $(iman)
	cp $^ $@

$(ishare)%:	$(share)%
	mkdir -p $(ishare)
	cp $^ $@

###

.PHONY:	uninstall
uninstall:
	-rm $(IBIN) $(IINC) $(ILIB) $(IMAN) $(ISHARE)
