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

void yuv2rgb(float *rgb, float *yuv)
{
	float WR = 0.2126f;
	float WB = 0.0722f;
	float WG = 1.0f - WR - WB;
	float UMAX = 0.436f;
	float VMAX = 0.615f;
	rgb[0] = fclampf(yuv[0] + (1.0f - WR) / VMAX * yuv[2], 0.0f, 1.0f);
	rgb[1] = fclampf(yuv[0] - WB * (1.0f - WB) / (UMAX * WG) * yuv[1] - WR * (1.0f - WR) / (VMAX * WG) * yuv[2], 0.0f, 1.0f);
	rgb[2] = fclampf(yuv[0] + (1.0f - WB) / UMAX * yuv[1], 0.0f, 1.0f);
}

void rgb2yuv(float *yuv, float *rgb)
{
	float WR = 0.2126f;
	float WB = 0.0722f;
	float WG = 1.0f - WR - WB;
	float UMAX = 0.436f;
	float VMAX = 0.615f;
	yuv[0] = fclampf(WR * rgb[0] + WG * rgb[1] + WB * rgb[2], 0.0f, 1.0f);
	yuv[1] = fclampf(UMAX / (1.0f - WB) * (rgb[2] - yuv[0]), -UMAX, UMAX);
	yuv[2] = fclampf(VMAX / (1.0f - WR) * (rgb[0] - yuv[0]), -VMAX, VMAX);
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

void blah(fftwf_plan DCTII, fftwf_plan DCTIII, float *td, float *fd, float *io, int N, float f)
{
	for (int i = 0; i < N * N; i++)
		td[i] = io[3 * i];

	fftwf_execute(DCTII);

	for (int i = 0; i < N * N; i++)
		fd[i] = roundf(f / (2.0f * N) * fd[i]) / f;

	fftwf_execute(DCTIII);

	for (int i = 0; i < N * N; i++)
		io[3 * i] = td[i] / (2.0f * N);
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
	for (int tj = 0; tj < th; tj++) {
		for (int ti = 0; ti < tw; ti++) {
			float yuv[3 * N * N];
			for (int j = 0; j < N; j++) {
				for (int i = 0; i < N; i++) {
					int idx = w * N * tj + w * j + N * ti + i;
					rgb2yuv(yuv + 3 * (N * j + i), ib + 3 * idx);
				}
			}

			blah(DCTII, DCTIII, td, fd, yuv+0, N, 64);
			blah(DCTII, DCTIII, td, fd, yuv+1, N, 16);
			blah(DCTII, DCTIII, td, fd, yuv+2, N, 16);

			for (int j = 0; j < N; j++) {
				for (int i = 0; i < N; i++) {
					int idx = w * N * tj + w * j + N * ti + i;
					yuv2rgb(ob + 3 * idx, yuv + 3 * (N * j + i));
				}
			}
		}
	}
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

