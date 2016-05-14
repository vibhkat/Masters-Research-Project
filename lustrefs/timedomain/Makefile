all: main.o stamping.o functions.o
	CC main.o stamping.o functions.o -O1 -ffast-math -o result

main.o: main.cpp
	CC -c main.cpp

stamping.o: stamping.cpp
	CC -c stamping.cpp

functions.o: functions.cpp
	CC -c functions.cpp

clean:
	rm -rf *.o result

    
