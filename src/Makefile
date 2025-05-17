CC = g++
# CFLAGS = -std=c++11 -O3 -I /usr/local/include/boost/
CFLAGS = -std=c++11 -g

OBJS = main.o datastructure.o

PJ4: $(OBJS)
	$(CC) $(OBJS) -o ../bin/PJ4

main.o: main.cpp datastructure.hpp
	$(CC) -c main.cpp $(CFLAGS)

datastructure.o: datastructure.cpp datastructure.hpp
	$(CC) -c datastructure.cpp $(CFLAGS)

clean:
	rm -f *.o ../bin/PJ4
