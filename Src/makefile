all: load gdelay bundle kvecs kflow rnd phix heap

load: load.o fft.o wave1d.o array3d.o
	g++ -o load load.o fft.o wave1d.o array3d.o
gdelay: gdelay.o fft.o wave1d.o
	g++ -o gdelay gdelay.o fft.o wave1d.o
bundle: bundle.o load.o fft.o wave1d.o array3d.o
	g++ -o bundle bundle.o array3d.o fft.o wave1d.o
kvecs: kvecs.o fft.o wave1d.o array3d.o fslice.o
	g++ -o kvecs kvecs.o array3d.o fft.o wave1d.o fslice.o
phix: phix.o fft.o wave1d.o array3d.o fslice.o
	g++ -o phix phix.o array3d.o fft.o wave1d.o fslice.o
#kflow: kflow.o array3d.o wave1d.o fft.o
#	g++ -o kflow kflow.o array3d.o wave1d.o fft.o
kflow: kflow.o  fft.o array3d.o wave1d.o 
	g++ -o kflow kflow.o  fft.o array3d.o wave1d.o 
rnd: rnd.o fslice.o
	g++ -o rnd rnd.o fslice.o
heap: heap.o
	g++ -o heap heap.o

heap.o: heap.cpp
	g++ -c heap.cpp
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
phix.o: phix.cpp
	g++ -c phix.cpp
kflow.o: kflow.cpp
	g++ -c kflow.cpp
fslice.o: fslice.cpp
	g++ -c fslice.cpp
rnd.o: rnd.cpp
	g++ -c rnd.cpp
