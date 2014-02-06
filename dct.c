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

struct rgb {
	float r, g, b;
};

struct image {
	struct rgb *buffer;
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
	image->buffer = malloc(sizeof(struct rgb) * width * height);
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
	for (int i = 0; i < image->total; i++) {
		int r = fgetc(file);
		int g = fgetc(file);
		int b = fgetc(file);
		if (EOF == r || EOF == g || EOF == b)
			goto eof;
		image->buffer[i] = (struct rgb){ linear(r / 255.0f), linear(g / 255.0f), linear(b / 255.0f) };
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
	for (int i = 0; i < image->total; i++) {
		if (EOF == fputc(255.0f * srgb(image->buffer[i].r), file))
			goto eof;
		if (EOF == fputc(255.0f * srgb(image->buffer[i].g), file))
			goto eof;
		if (EOF == fputc(255.0f * srgb(image->buffer[i].b), file))
			goto eof;
	}
	fclose(file);
	return 1;
eof:
	fprintf(stderr, "EOF while writing to \"%s\".\n", image->name);
	fclose(file);
	return 0;
}

struct yuv {
	float y, u, v;
};

float fclampf(float x, float a, float b)
{
	return fminf(fmaxf(x, a), b);
}

struct rgb yuv2rgb(struct yuv c)
{
	float WR = 0.2126f;
	float WB = 0.0722f;
	float WG = 1.0f - WR - WB;
	float UMAX = 0.436f;
	float VMAX = 0.615f;
	return (struct rgb) {
		fclampf(c.y + (1.0f - WR) / VMAX * c.v, 0.0f, 1.0f),
		fclampf(c.y - WB * (1.0f - WB) / (UMAX * WG) * c.u - WR * (1.0f - WR) / (VMAX * WG) * c.v, 0.0f, 1.0f),
		fclampf(c.y + (1.0f - WB) / UMAX * c.u, 0.0f, 1.0f)
	};
}

struct yuv rgb2yuv(struct rgb c)
{
	float WR = 0.2126f;
	float WB = 0.0722f;
	float WG = 1.0f - WR - WB;
	float UMAX = 0.436f;
	float VMAX = 0.615f;
	float y = WR * c.r + WG * c.g + WB * c.b;
	return (struct yuv) {
		fclampf(y, 0.0f, 1.0f),
		fclampf(UMAX / (1.0f - WB) * (c.b - y), -UMAX, UMAX),
		fclampf(VMAX / (1.0f - WR) * (c.r - y), -VMAX, VMAX),
	};
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
	struct rgb *ob = output->buffer;
	struct rgb *ib = input->buffer;
	for (int tj = 0; tj < th; tj++) {
		for (int ti = 0; ti < tw; ti++) {
			struct yuv yuv[N * N];
			for (int j = 0; j < N; j++) {
				for (int i = 0; i < N; i++) {
					int idx = w * N * tj + w * j + N * ti + i;
					yuv[N * j + i] = rgb2yuv(ib[idx]);
				}
			}
			// Y
			for (int i = 0; i < N * N; i++)
				td[i] = yuv[i].y;

			fftwf_execute(DCTII);

			float yf = 64;
			for (int i = 0; i < N * N; i++)
				fd[i] = roundf(yf / (2.0f * N) * fd[i]) / yf;

			fftwf_execute(DCTIII);

			for (int i = 0; i < N * N; i++)
				yuv[i].y = td[i] / (2.0f * N);

			// U
			for (int i = 0; i < N * N; i++)
				td[i] = yuv[i].u;

			fftwf_execute(DCTII);

			float uf = 16;
			for (int i = 0; i < N * N; i++)
				fd[i] = roundf(uf / (2.0f * N) * fd[i]) / uf;

			fftwf_execute(DCTIII);

			for (int i = 0; i < N * N; i++)
				yuv[i].u = td[i] / (2.0f * N);

			// V
			for (int i = 0; i < N * N; i++)
				td[i] = yuv[i].v;

			fftwf_execute(DCTII);

			float vf = 16;
			for (int i = 0; i < N * N; i++)
				fd[i] = roundf(vf / (2.0f * N) * fd[i]) / vf;

			fftwf_execute(DCTIII);

			for (int i = 0; i < N * N; i++)
				yuv[i].v = td[i] / (2.0f * N);

			for (int j = 0; j < N; j++) {
				for (int i = 0; i < N; i++) {
					int idx = w * N * tj + w * j + N * ti + i;
					ob[idx] = yuv2rgb(yuv[N * j + i]);
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

