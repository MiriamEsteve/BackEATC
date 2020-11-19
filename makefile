programa: main.o eat.o vlist.o vtree.o distrib.o data.o
	gcc -o programa main.o eat.o vlist.o vtree.o distrib.o data.o -lm
	rm -f *.o

main.o: main.c eat.h ../00lib/vtree.h ../00lib/distrib.h ../00lib/data.h
	gcc -c main.c -lm

eat.o: eat.c eat.h ../00lib/data.h
	gcc -c eat.c -lm

vlist.o: ../00lib/vlist.c ../00lib/vlist.h
	gcc -c ../00lib/vlist.c -lm

vtree.o: ../00lib/vtree.c ../00lib/vtree.h
	gcc -c ../00lib/vtree.c -lm

distrib.o: ../00lib/distrib.c ../00lib/distrib.h
	gcc -c ../00lib/distrib.c -lm

data.o: ../00lib/data.c ../00lib/data.h
	gcc -c ../00lib/data.c -lm

clean:
	rm -f programa *.o

