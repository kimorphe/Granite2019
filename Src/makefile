all: load gdelay

load: load.o fft.o wave1d.o
	g++ -o load load.o fft.o wave1d.o
fft.o: fft.cpp 
	g++ -c fft.cpp
load.o: load.cpp
	g++ -c load.cpp
wave1d.o: wave1d.cpp
	g++ -c wave1d.cpp

gdelay: gdelay.o fft.o wave1d.o
	g++ -o gdelay gdelay.o fft.o wave1d.o
gdelay.o: gdelay.cpp
	g++ -c gdelay.cpp
