name = pqgen
cmds = pq-theta dna2codon
headers = pq_sfobj pq_sfstats
libname = lib$(name)

HOME = $(shell echo $$HOME)/
BASE = $(HOME).local/
INSTALL = $(HOME).local/

###

bin = bin/
inc = include/
lib = lib/
obj = obj/
src = src/

bins = $(addprefix $(bin),$(cmds))
incs = $(addsuffix .h,$(addprefix $(inc),$(headers)))
libs = $(lib)$(libname).so
objs = $(addsuffix .o,$(addprefix $(obj),$(headers)))

######

.PHONY:	all
all:	$(bins) $(libs) $(objs)

$(bin)%:	$(src)%.c $(libs)
	mkdir -p $(bin)
	gcc -I include/ -I $(BASE)include/librawk/ -L $(BASE)lib/ -L $(lib) -o $@ $(word 1,$^) -lrawk -lm -lpqgen

$(libs):	$(objs)
	mkdir -p $(lib)
	gcc -shared -o $@ $^

$(obj)%.o:	$(src)%.c
	mkdir -p $(obj)
	gcc -I $(inc) -c -Wall -fpic -o $@ $^

######

ibin = $(INSTALL)$(bin)
ilib = $(INSTALL)$(lib)
iinc = $(INSTALL)$(inc)$(libname)/

IBIN = $(addprefix $(ibin),$(cmds))
IINC = $(addsuffix .h,$(addprefix $(iinc),$(headers)))
ILIB = $(addprefix $(INSTALL),$(libs))

.PHONY:	install
install:	$(IBIN) $(IINC) $(ILIB)

$(ibin)%:	$(bin)%
	mkdir -p $(ibin)
	cp $^ $@

$(iinc)%.h:	$(inc)%.h
	mkdir -p $(iinc)
	cp $^ $@

$(ilib)%.so:	$(lib)%.so
	mkdir -p $(ilib)
	cp $^ $@

.PHONY:	clean
clean:
	-rm $(objs) $(libs)
