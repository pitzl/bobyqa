
lib: bobyqa.cc bobyqa.h
	g++ -O2 -Wall bobyqa.cc -fPIC -shared -o libbobyqa.so
