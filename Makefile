HOME = $(shell echo $$HOME)/
NAME = pqgen
DOCS = pqgen.1
CMDS = pqgen pq-dna2codon pq-gen2hwe
HEADERS = pq_args pq_div pq_generics pq_genetics pq_het pq_parse pq_pnds pq_sfs pq_sfstats pq_theta pq_version
DEPINC = $(HOME).local/
PREFIX = $(HOME).local/
CONFIG = $(HOME).config/$(NAME)/

Wgcc = -Wall -Wextra -Wpedantic

bin = bin/
inc = include/
lib = lib/
man = doc/
share = share/
src = src/

bins = $(addprefix $(bin),$(CMDS))
srcs = $(addprefix $(src),$(addsuffix .c,$(HEADERS)))


######

.PHONY:	all
all:	binaries documentation


######

.PHONY:	binaries
binaries:	$(bins)

$(bin)%:	$(src)%.c $(srcs)
	mkdir -p $(bin)
	gcc -I $(inc) -I $(DEPINC)$(inc)librawk/ -L $(DEPINC)$(lib) $(Wgcc) -o $@ $^ -lrawk -lm



.PHONY:	documentation
documentation:	$(man)pqgen.1

$(man)pqgen.1:	$(man)pqgen-TH $(man)pqgen-body
	cat $^ > $@

$(man)pqgen-TH:	$(bin)pqgen
	./$^ --version | cut -d' ' -f 1,3 | awk ' { print ".TH pq-genetics 1 \""strftime("%Y-%m-%d")"\" \""$$0"\" \"Population and Quantitative Genetic Tools\""} ' > $@


######

.PHONY:	clean
clean:
	-rm $(bins) $(mans)


######

.PHONY:	install
install:	$(PREFIX)bin/pqgen $(PREFIX)share/man/man1/pqgen.1 $(CONFIG)codon2aa $(CONFIG)codon2syn

$(PREFIX)bin/%:	$(bin)%
	mkdir -p $(prefix)bin/
	cp $^ $@

$(PREFIX)share/man/man1/%:	$(man)%
	mkdir -p $(PREFIX)share/man/man1/
	cp $^ $@

$(CONFIG)%:	share/%
	mkdir -p $(CONFIG)
	cp $^ $@


######

.PHONY:	uninstall
uninstall:
	-rm $(PREFIX)bin/pqgen $(PREFIX)share/man/man1/pqgen.1 $(CONFIG)codon2aa $(CONFIG)codon2syn

