CC=mpic++
CFLAGS=
LDFLAGS=
SOURCES=src/*
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=bin/N-TR0PY

make: src/*
	mpic++ -Wall src/* -o bin/N-TR0PY
