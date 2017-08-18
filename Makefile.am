ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -DDATADIR=\"$(pkgdatadir)\"

bin_PROGRAMS = bin/pqgen bin/dna2codon bin/gen2hwe bin/codon2pnds

bin_pqgen_SOURCES = \
src/pqgen.c\
src/pq_div.c src/pq_het.c src/pq_pnds.c src/pq_sfs.c src/pq_theta.c\
src/pq_args.c src/pq_generics.c src/pq_genetics.c\
src/pq_parse.c src/pq_sfstats.c src/pq_version.c\
include/pq_args.h include/pq_div.h include/pq_generics.h include/pq_genetics.h\
include/pq_help.h include/pq_het.h include/pq_limits.h include/pq_parse.h\
include/pq_pnds.h include/pq_sfs.h include/pq_sfstats.h include/pq_theta.h include/pq_version.h

bin_codon2pnds_SOURCES = \
src/codon2pnds.c\
src/pq_div.c src/pq_het.c src/pq_pnds.c src/pq_sfs.c src/pq_theta.c\
src/pq_args.c src/pq_generics.c src/pq_genetics.c\
src/pq_parse.c src/pq_sfstats.c src/pq_version.c\
include/pq_args.h include/pq_div.h include/pq_generics.h include/pq_genetics.h\
include/pq_help.h include/pq_het.h include/pq_limits.h include/pq_parse.h\
include/pq_pnds.h include/pq_sfs.h include/pq_sfstats.h include/pq_theta.h include/pq_version.h

bin_dna2codon_SOURCES = \
src/pq-dna2codon.c

bin_gen2hwe_SOURCES = \
src/pq-gen2hwe.c

bin_pqgen_CFLAGS = -I $(srcdir)/include -I $(includedir)/librawk
bin_dna2codon_CFLAGS = -I $(srcdir)/include -I $(includedir)/librawk
bin_codon2pnds_CFLAGS = -I $(srcdir)/include -I $(includedir)/librawk
bin_gen2hwe_CFLAGS = -I $(srcdir)/include -I $(includedir)/librawk

dist_pkgdata_DATA = share/codon2aa share/codon2syn
