all: main.o print.o stamping.o functions.o
	CC main.o print.o stamping.o functions.o -O1 -ffast-math 

main.o: main.cpp
	CC -c main.cpp

print.o: print.cpp
	CC -c print.cpp

stamping.o: stamping.cpp
	CC -c stamping.cpp

functions.o: functions.cpp
	CC -c functions.cpp 
	
clean:
	rm -rf *.o frequency
