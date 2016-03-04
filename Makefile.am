AUTOMAKE_OPTIONS = foreign
CFLAGS += -ggdb -Wall -march=core2
CXXFLAGS += -ggdb -Wall -march=core2
LDFLAGS += -lpthread

bin_PROGRAMS = rfmix simulate
rfmix_SOURCES = cmdline-utils.c inputline.cpp genetic-map.cpp md5rng.cpp hash-table.cpp mm.cpp rfmix.cpp load-input.cpp random-forest.cpp crf.cpp output.cpp gensamples.cpp s-sample.cpp

simulate_SOURCES = cmdline-utils.c inputline.cpp genetic-map.cpp hash-table.cpp vcf.cpp simulate.cpp s-sample.cpp s-subpop.cpp

