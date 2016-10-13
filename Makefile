name = pqgen
cmds = pqgen pq-dna2codon
docs = pqgen.1 pq-theta.1 pq-dna2codon.1
headers = pq_args pq_div pq_generics pq_genetics pq_pnds pq_sfstats pq_parse pq_theta
libname = lib$(name)
info = codon2aa codon2syn

HOME = $(shell echo $$HOME)/
BASE = $(HOME).local/
INSTALL = $(HOME).local/

CONFIG = $(HOME).config/$(name)/

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
srcs = $(addprefix $(src),$(addsuffix .c,$(headers)))

######

.PHONY:	all
all:	$(bins) $(libs) $(man)pqgen.1

$(bin)%:	$(src)%.c $(srcs)
	mkdir -p $(bin)
	gcc -I include/ -I $(BASE)include/librawk/ -L $(BASE)lib/ $(Wgcc) -o $@ $^ -lrawk -lm

$(libs):	$(objs)
	mkdir -p $(lib)
	gcc -shared -o $@ $^

$(obj)%.o:	$(src)%.c
	mkdir -p $(obj)
	gcc -I $(inc) -I $(BASE)include/librawk/ -L $(BASE)lib/ -c $(Wgcc) -fpic -o $@ $^ -lrawk

$(man)pqgen.1:	$(man)pqgen-TH $(man)pqgen-body
	cat $^ > $@

$(man)pqgen-TH:	$(bin)pqgen
	./$^ --version | cut -d' ' -f 1,3 | awk ' { print ".TH pq-genetics 1 \""strftime("%Y-%m-%d")"\" \""$$0"\" \"Population and Quantitative Genetic Tools\""} ' > $@

###

.PHONY:	clean
clean:
	-rm $(bins) $(objs) $(libs)

######

ibin = $(INSTALL)$(bin)
ilib = $(INSTALL)$(lib)
iman = $(INSTALL)share/man/man1/
iinc = $(INSTALL)$(inc)$(libname)/
ishare = $(INSTALL)share/$(libname)/
iconfig = $(CONFIG)

IBIN = $(addprefix $(ibin),$(cmds))
IINC = $(addsuffix .h,$(addprefix $(iinc),$(headers)))
ILIB = $(addprefix $(INSTALL),$(libs))
IMAN = $(addprefix $(iman),$(docs))
ISHARE = $(addprefix $(ishare),$(info))
ICONFIG = $(addprefix $(iconfig),$(info))

.PHONY:	install
install:	$(IBIN) $(IINC) $(ILIB) $(IMAN) $(ICONFIG)

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

$(iconfig)%:	$(share)%
	mkdir -p $(iconfig)
	cp $^ $@

###

.PHONY:	uninstall
uninstall:
	-rm $(IBIN) $(IINC) $(ILIB) $(IMAN) $(ICONFIG)
