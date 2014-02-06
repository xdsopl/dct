CFLAGS = -std=c99 -W -Wall -O3 -D_GNU_SOURCE=1 -g
LDFLAGS = -lm $(shell pkg-config fftw3f --libs)

all: dct

test: dct
	./dct input.ppm output.ppm

dct: dct.c

clean:
	rm -f dct

