/*
dct - playing with dct and lossy image compression
Written in 2014 by <Ahmet Inan> <xdsopl@googlemail.com>
To the extent possible under law, the author(s) have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

float srgb(float v)
{
	float K0 = 0.03928f;
	float a = 0.055f;
	float phi = 12.92f;
	float gamma = 2.4f;
	return v <= K0 / phi ? v * phi : (1.0f + a) * powf(v, 1.0f / gamma) - a;
}

float linear(float v)
{
	float K0 = 0.03928f;
	float a = 0.055f;
	float phi = 12.92f;
	float gamma = 2.4f;
	return v <= K0 ? v / phi : powf((v + a) / (1.0f + a), gamma);
}

float fclampf(float x, float a, float b)
{
	return fminf(fmaxf(x, a), b);
}

void ycbcr2rgb(float *io)
{
	float WR = 0.2126f;
	float WB = 0.0722f;
	float WG = 1.0f - WR - WB;
	float UMAX = 0.5f;
	float VMAX = 0.5f;
	float y = io[0], u = io[1], v = io[2];
	io[0] = fclampf(y + (1.0f - WR) / VMAX * v, 0.0f, 1.0f);
	io[1] = fclampf(y - WB * (1.0f - WB) / (UMAX * WG) * u - WR * (1.0f - WR) / (VMAX * WG) * v, 0.0f, 1.0f);
	io[2] = fclampf(y + (1.0f - WB) / UMAX * u, 0.0f, 1.0f);
}

void rgb2ycbcr(float *io)
{
	float WR = 0.2126f;
	float WB = 0.0722f;
	float WG = 1.0f - WR - WB;
	float UMAX = 0.5f;
	float VMAX = 0.5;
	float r = io[0], g = io[1], b = io[2];
	io[0] = fclampf(WR * r + WG * g + WB * b, 0.0f, 1.0f);
	io[1] = fclampf(UMAX / (1.0f - WB) * (b - io[0]), -UMAX, UMAX);
	io[2] = fclampf(VMAX / (1.0f - WR) * (r - io[0]), -VMAX, VMAX);
}

struct image {
	float *buffer;
	int width, height, total;
	char *name;
};

void delete_image(struct image *image)
{
	free(image->buffer);
	free(image);
}

struct image *new_image(char *name, int width, int height)
{
	struct image *image = malloc(sizeof(struct image));
	image->height = height;
	image->width = width;
	image->total = width * height;
	image->name = name;
	image->buffer = malloc(3 * sizeof(float) * width * height);
	return image;
}

void ycbcr_image(struct image *image)
{
	for (int i = 0; i < image->total; i++)
		rgb2ycbcr(image->buffer + 3 * i);
}

void rgb_image(struct image *image)
{
	for (int i = 0; i < image->total; i++)
		ycbcr2rgb(image->buffer + 3 * i);
}

struct image *read_ppm(char *name)
{
	FILE *file = fopen(name, "r");
	if (!file) {
		fprintf(stderr, "could not open \"%s\" file to read.\n", name);
		return 0;
	}
	if ('P' != fgetc(file) || '6' != fgetc(file)) {
		fprintf(stderr, "file \"%s\" not P6 image.\n", name);
		fclose(file);
		return 0;
	}
	int integer[3];
	struct image *image = 0;
	int c = fgetc(file);
	if (EOF == c)
		goto eof;
	for (int i = 0; i < 3; i++) {
		while ('#' == (c = fgetc(file)))
			while ('\n' != (c = fgetc(file)))
				if (EOF == c)
					goto eof;
		while ((c < '0') || ('9' < c))
			if (EOF == (c = fgetc(file)))
				goto eof;
		char str[16];
		for (int n = 0; n < 16; n++) {
			if (('0' <= c) && (c <= '9') && n < 15) {
				str[n] = c;
				if (EOF == (c = fgetc(file)))
					goto eof;
			} else {
				str[n] = 0;
				break;
			}
		}
		integer[i] = atoi(str);
	}
	if (!(integer[0] && integer[1] && integer[2])) {
		fprintf(stderr, "could not read image file \"%s\".\n", name);
		fclose(file);
		return 0;
	}
	if (integer[2] != 255) {
		fprintf(stderr, "cant read \"%s\", only 8 bit per channel SRGB supported at the moment.\n", name);
		fclose(file);
		return 0;
	}
	image = new_image(name, integer[0], integer[1]);
	for (int i = 0; i < 3 * image->total; i++) {
		int v = fgetc(file);
		if (EOF == v)
			goto eof;
		image->buffer[i] = linear(v / 255.0f);
	}
	fclose(file);
	return image;
eof:
	fprintf(stderr, "EOF while reading from \"%s\".\n", name);
	fclose(file);
	delete_image(image);
	return 0;
}

int write_ppm(struct image *image)
{
	FILE *file = fopen(image->name, "w");
	if (!file) {
		fprintf(stderr, "could not open \"%s\" file to write.\n", image->name);
		return 0;
	}
	if (!fprintf(file, "P6 %d %d 255\n", image->width, image->height)) {
		fprintf(stderr, "could not write to file \"%s\".\n", image->name);
		fclose(file);
		return 0;
	}
	for (int i = 0; i < 3 * image->total; i++) {
		if (EOF == fputc(255.0f * srgb(image->buffer[i]), file))
			goto eof;
	}
	fclose(file);
	return 1;
eof:
	fprintf(stderr, "EOF while writing to \"%s\".\n", image->name);
	fclose(file);
	return 0;
}

float quantization(int i, int j, int N, float min, float max)
{
	(void)N;
	return min + (max - min) * powf(2.0f, - (i + j));
//	return min + (max - min) * (2 * N - (i + j)) / (2.0f * N);
//	return min + (max - min) * (N * N - (i * j)) / (float)(N * N);
}

void blah(fftwf_plan DCTII, fftwf_plan DCTIII, float *td, float *fd, float *io, int N, float min, float max)
{
	for (int i = 0; i < N * N; i++)
		td[i] = io[3 * i];

	fftwf_execute(DCTII);

	for (int i = 0; i < N * N; i++)
		fd[i] /= 2.0f * N;

	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++) {
			float q = quantization(i, j, N, min, max);
			fd[N * j + i] = roundf(q * fd[N * j + i]) / q;
		}
	}

	fftwf_execute(DCTIII);

	for (int i = 0; i < N * N; i++)
		td[i] /= 2.0f * N;

	for (int i = 0; i < N * N; i++)
		io[3 * i] = td[i];
}

void probe(float *t, int N, int i, int j, int pi, int pj)
{
	if (i != pi || j != pj)
		return;
	printf("\n");
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++)
			printf(" % 12.7f ", t[N * j + i]);
		printf("\n");
	}
}

void doit(struct image *output, struct image *input)
{
	int N = 8;
	float td[N * N], fd[N * N];
	fftwf_plan DCTII = fftwf_plan_r2r_2d(N, N, td, fd, FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
	fftwf_plan DCTIII = fftwf_plan_r2r_2d(N, N, fd, td, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
	int w = output->width;
	int h = output->height;
	int tw = w / N;
	int th = h / N;
	float *ob = output->buffer;
	float *ib = input->buffer;
	ycbcr_image(input);
	for (int tj = 0; tj < th; tj++) {
		for (int ti = 0; ti < tw; ti++) {
			float tile[3 * N * N];
			for (int j = 0; j < N; j++)
				memcpy(tile + 3 * N * j, ib + 3 * (N * (w * tj + ti) + w * j), 3 * sizeof(float) * N);

			int pi = tw / 2 - 1;
			int pj = th / 2 - 1;
			blah(DCTII, DCTIII, td, fd, tile+0, N, 16, 64);
			probe(fd, N, ti, tj, pi, pj);
			blah(DCTII, DCTIII, td, fd, tile+1, N, 4, 16);
			probe(fd, N, ti, tj, pi, pj);
			blah(DCTII, DCTIII, td, fd, tile+2, N, 4, 16);
			probe(fd, N, ti, tj, pi, pj);

			for (int j = 0; j < N; j++)
				memcpy(ob + 3 * (N * (w * tj + ti) + w * j), tile + 3 * N * j, 3 * sizeof(float) * N);
		}
	}
	rgb_image(output);
}

int main(int argc, char **argv)
{
	if (argc != 3) {
		fprintf(stderr, "usage: %s input.ppm output.ppm\n", argv[0]);
		return 1;
	}
	struct image *input = read_ppm(argv[1]);
	if (!input)
		return 1;
	struct image *output = new_image(argv[2], input->width, input->height);

	doit(output, input);

	if (!write_ppm(output))
		return 1;
	return 0;
}

