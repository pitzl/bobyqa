
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)

student: student.cpp Makefile
	g++ -O2 -Wall $(ROOTCFLAGS) student.cpp  -o student \
	$(ROOTLIBS) -L.. -lbobyqa
	@echo 'done: student'
	@echo 'export LD_LIBRARY_PATH=sLD_LIBRARY_PATH:..'

rosenb: rosenb.cpp Makefile
	g++ -O2 -Wall $(ROOTCFLAGS) rosenb.cpp  -o rosenb \
	$(ROOTLIBS) -L.. -lbobyqa
	@echo 'done: rosenb'
	@echo 'export LD_LIBRARY_PATH=sLD_LIBRARY_PATH:..'
