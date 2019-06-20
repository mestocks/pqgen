ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -DDATADIR=\"$(pkgdatadir)\"

bin_PROGRAMS = bin/pqgen bin/dna2codon bin/dna2div bin/gen2hwe bin/codon2pnds

bin_pqgen_SOURCES = \
src/main.c\
src/cmd_het.c src/cmd_sfs.c src/cmd_theta.c\
src/pq_div.c src/pq_pnds.c\
src/pq_args.c src/pq_generics.c src/pq_genetics.c src/pq_htable.c\
src/pq_parse.c src/pq_sfstats.c src/pq_version.c\
include/pq_args.h include/pq_div.h include/pq_generics.h include/pq_genetics.h include/pq_htable.h\
include/pq_help.h include/cmd_het.h include/pq_limits.h include/pq_parse.h\
include/pq_pnds.h include/cmd_sfs.h include/pq_sfstats.h include/cmd_theta.h include/pq_version.h

bin_codon2pnds_SOURCES = \
src/codon2pnds.c\
src/pq_div.c src/cmd_het.c src/pq_pnds.c src/cmd_sfs.c src/cmd_theta.c\
src/pq_args.c src/pq_generics.c src/pq_genetics.c src/pq_htable.c\
src/pq_parse.c src/pq_sfstats.c src/pq_version.c\
include/pq_args.h include/pq_div.h include/pq_generics.h include/pq_genetics.h include/pq_htable.h\
include/pq_help.h include/cmd_het.h include/pq_limits.h include/pq_parse.h\
include/pq_pnds.h include/cmd_sfs.h include/pq_sfstats.h include/cmd_theta.h include/pq_version.h

bin_dna2codon_SOURCES = \
src/dna2codon.c src/pq_parse.c include/pq_parse.h

bin_dna2div_SOURCES = \
src/dna2div.c\
src/pq_div.c src/cmd_het.c src/pq_pnds.c src/cmd_sfs.c src/cmd_theta.c\
src/pq_args.c src/pq_generics.c src/pq_genetics.c src/pq_htable.c\
src/pq_parse.c src/pq_sfstats.c src/pq_version.c\
include/pq_args.h include/pq_div.h include/pq_generics.h include/pq_genetics.h\
include/pq_help.h include/cmd_het.h include/pq_limits.h include/pq_parse.h\
include/pq_pnds.h include/cmd_sfs.h include/pq_sfstats.h include/cmd_theta.h include/pq_version.h

bin_gen2hwe_SOURCES = \
src/gen2hwe.c src/pq_parse.c include/pq_parse.h

bin_pqgen_CFLAGS = -I $(srcdir)/include
bin_dna2codon_CFLAGS = -I $(srcdir)/include
bin_dna2div_CFLAGS = -I $(srcdir)/include
bin_codon2pnds_CFLAGS = -I $(srcdir)/include
bin_gen2hwe_CFLAGS = -I $(srcdir)/include

dist_pkgdata_DATA = share/codon2aa share/codon2syn

man_MANS = doc/pqgen.1 doc/pqgen-het.1 doc/pqgen-sfs.1 doc/pqgen-theta.1
dist_man_MANS = doc/pqgen.1 doc/pqgen-het.1 doc/pqgen-sfs.1 doc/pqgen-theta.1
