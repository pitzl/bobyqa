
lib: bobyqa.cc bobyqa.h Makefile
	g++ -O2 -Wall -Wextra bobyqa.cc -shared -fPIC -o libbobyqa.so
