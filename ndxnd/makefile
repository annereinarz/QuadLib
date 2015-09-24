all: cubequad simplexquad

cubequad: src/*.c src/CubexCube/*.c
	gcc -std=c99 -Wall -lm -g -o $@ -I src -I src/CubexCube $+

simplexquad: src/*.c src/SimpxSimp/*.c
	gcc -std=c99 -Wall -lm -g -o $@ -I src -I src/SimpxSimp $+

