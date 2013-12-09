/*----------------------------------------------------------------------------
 
 RHF - Ray Histogram Fusion
 
 Copyright (c) 2013, A. Buades <toni.buades@uib.es>,
 M. Delbracio <mdelbra@gmail.com>, 
 J-M. Morel <morel@cmla.ens-cachan.fr>,
 P. Muse <muse@fing.edu.uy>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU Affero General Public License for more details.
 
 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 
 ----------------------------------------------------------------------------*/

#ifndef _IO_PNG_H
#define _IO_PNG_H

#ifdef __cplusplus
extern "C" {
#endif

#define IO_PNG_VERSION "0.20110825"

#include <stddef.h>

/* io_png.c */
char *io_png_info(void);
unsigned char *io_png_read_u8(const char *fname, size_t *nxp, size_t *nyp, size_t *ncp);
unsigned char *io_png_read_u8_rgb(const char *fname, size_t *nxp, size_t *nyp);
unsigned char *io_png_read_u8_gray(const char *fname, size_t *nxp, size_t *nyp);
float *io_png_read_f32(const char *fname, size_t *nxp, size_t *nyp, size_t *ncp);
float *io_png_read_f32_rgb(const char *fname, size_t *nxp, size_t *nyp);
float *io_png_read_f32_gray(const char *fname, size_t *nxp, size_t *nyp);
int io_png_write_u8(const char *fname, const unsigned char *data, size_t nx, size_t ny, size_t nc);
int io_png_write_f32(const char *fname, const float *data, size_t nx, size_t ny, size_t nc);

#ifdef __cplusplus
}
#endif

#endif /* !_IO_PNG_H */
