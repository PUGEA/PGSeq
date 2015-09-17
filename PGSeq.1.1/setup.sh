#!/bin/bush
gcc -c -lm -fPIC -I/usr/include/python2.7 donlp2.c newx.c pgseq.c user_eval.c wrap.c
gcc -shared -fPIC -o example.so donlp2.o newx.o pgseq.o user_eval.o wrap.o

