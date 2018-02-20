# Makefile for the Montecarlo simulation exercise
# Jimmy Aguilar Mena
# 20/02/2018

all: dispersion.x

dispersion.x: main.c histogram.o
	gcc $^ -o $@ -lm

%.o: %.c
	gcc -c $< -o $@

.PHONY: clean test

clean:
	rm -r *.o *.x *.out

test: dispersion.x
	./$< 16 5 2 1000000
