all: project1

project1: project1.c
	gcc -o project1 project1.c -lm

clean:
	rm *~ project1
