#
# Makefile for HW2
#

CC=g++
# if you want to use debugger, add -g to CFLAGS and LDFLAGS
CFLAGS=
LDFLAGS=-std=c++11 -O3 -lm
SOURCES=./src/main.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=./bin/hw1
INCLUDES=./src/cell.h ./src/net.h

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o:  %.c  ${INCLUDES}
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)
