all: load gdelay bundle kvecs

load: load.o fft.o wave1d.o array3d.o
	g++ -o load load.o fft.o wave1d.o array3d.o
gdelay: gdelay.o fft.o wave1d.o
	g++ -o gdelay gdelay.o fft.o wave1d.o
bundle: bundle.o load.o fft.o wave1d.o array3d.o
	g++ -o bundle bundle.o array3d.o fft.o wave1d.o
kvecs: kvecs.o fft.o wave1d.o array3d.o
	g++ -o kvecs kvecs.o array3d.o fft.o wave1d.o


fft.o: fft.cpp 
	g++ -c fft.cpp
load.o: load.cpp
	g++ -c load.cpp
wave1d.o: wave1d.cpp
	g++ -c wave1d.cpp
gdelay.o: gdelay.cpp
	g++ -c gdelay.cpp
bundle.o: bundle.cpp
	g++ -c bundle.cpp
array3d.o: array3d.cpp
	g++ -c array3d.cpp
kvecs.o: kvecs.cpp
	g++ -c kvecs.cpp
