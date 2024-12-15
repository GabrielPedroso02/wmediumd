/*
 *	wmediumd, wireless medium simulator for mac80211_hwsim kernel module
 *	Copyright (c) 2011 cozybit Inc.
 *
 *	Author: Javier Lopez    <jlopex@cozybit.com>
 *		Javier Cardona  <javier@cozybit.com>
 *
 *	This program is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU General Public License
 *	as published by the Free Software Foundation; either version 2
 *	of the License, or (at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *	02110-1301, USA.
 */

#include <sys/timerfd.h>
#include <libconfig.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>

#include "wmediumd.h"

#define MAP_AREA = 1500
int **vegetation_matrix = NULL;

void printNumberOfCalls(const char *filename, int x1, int y1, int x2, int y2, int call_count, int matrix_value) {
    FILE *file = fopen(filename, "w");

    if (file == NULL) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

	// debug

    fprintf(file, "Call_count: %d\n", call_count);
    fprintf(file, "X1: %d\n", x1);
    fprintf(file, "y1: %d\n", y1);
    fprintf(file, "X2: %d\n", x2);
    fprintf(file, "y2: %d\n", y2);
    fprintf(file, "Matrix value at that pos: %d\n", matrix_value);

    fclose(file);
}

void printWeissbergerOutput(const char *filename, int log_distance_path_loss, double PL, double freq, struct station *src, struct station *dst, int vegetation_depth) {
    FILE *file = fopen(filename, "w");

    if (file == NULL) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

	// debug

	fprintf(file, "Vegetation depth: %d\n", vegetation_depth);
    fprintf(file, "LogDistance: %d\n", log_distance_path_loss);
    fprintf(file, "Weissberger: %.1f\n", PL);
    fprintf(file, "freq: %.1f\n", freq);
    fprintf(file, "Sta source x=%.1f y=%.1f z=%.1f\n", src->x, src->y, src->z);
    fprintf(file, "Sta dest x=%.1f y=%.1f z=%.1f\n", dst->x, dst->y, dst->z);
    fprintf(file, "Sta src txpower:%d\n", src->tx_power);
    fprintf(file, "Sta dest txpower:%d\n", dst->tx_power);

    fclose(file);
}

void printVegetationDepth(const char *filename, int depth, double PL) {
    FILE *file = fopen(filename, "w");

    if (file == NULL) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

	// debug

    fprintf(file, "depth: %d\n", depth);
    fprintf(file, "PL: %.1f\n", PL);

    fclose(file);
}

void printErrorMessage(const char *filename, char *message, int matrix_size) {
    FILE *file = fopen(filename, "w");

    if (file == NULL) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

	// debug

    fprintf(file, "message: %s\n", message);
    fprintf(file, "matrix_size: %d\n", matrix_size);

    fclose(file);
}

void allocate_mem_vegetation_matrix(int matrix_size) {

    vegetation_matrix = (int **)malloc(matrix_size * sizeof(int *));
    if (vegetation_matrix == NULL) {
		// printErrorMessage("allocate_mem_error.txt", "Could not allocate memory for vegetation matrix rows", 1);
        fprintf(stderr, "Error: Could not allocate memory for vegetation matrix rows\n");
        exit(1);
    }

    for (int i = 0; i < matrix_size; i++) {
        vegetation_matrix[i] = (int *)malloc(matrix_size * sizeof(int));
        if (vegetation_matrix[i] == NULL) {
			// printErrorMessage("allocate_mem_error.txt", "Could not allocate memory for vegetation matrix columns", 2);
            fprintf(stderr, "Error: Could not allocate memory for vegetation matrix columns\n");
            exit(1);
        }
    }
}


void load_vegetation_matrix(const char *filename, int matrix_size) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
		// printErrorMessage("eror_opening_file.txt", filename, 0);
        fprintf(stderr, "Error: Could not open the matrix file\n");
        return;
    }

	allocate_mem_vegetation_matrix(matrix_size);

    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            fscanf(file, "%d", &vegetation_matrix[i][j]);
        }
    }

    fclose(file);
}

int calculate_vegetation_depth(int px1, int py1, int px2, int py2, int matrix_size) {
	//Bresenham's Line Algorithm
    int vdepth = 0;

	int x1 = px1; // to dodge C's shenanigans
	int y1 = py1;
	int x2 = px2;
	int y2 = py2;

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1; 
    int err = dx - dy;
	
	// static int call_count = 0;
	// call_count++;
	// char filename[30] = "a.txt";

	// if (call_count <= 9) // Ensure it's a single digit (valid for 1-9)
    // 	filename[0] = '0' + call_count;  // Convert int to char
	// sprintf(filename, "%d.txt", call_count);


    while (x1 != x2 || y1 != y2) {
        // Sum vegetation depth from the matrix at the current cell
        if (x1 >= 0 && x1 < matrix_size && y1 >= 0 && y1 < matrix_size) {
            vdepth += vegetation_matrix[y1][x1];
			// if (call_count == 1)
				// printNumberOfCalls("output_firstCalcDepth.txt", x1, y1, x2, y2, 1, vdepth);
			// printNumberOfCalls(filename, x1, y1, x2, y2, call_count, vegetation_matrix[y1][x1]);
        }

        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x1 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y1 += sy;
        }
    }

    // Add the last point (x2, y2)
    if (px2 >= 0 && px2 < matrix_size && py2 >= 0 && py2 < matrix_size) {
        vdepth += vegetation_matrix[py2][px2];
		// printNumberOfCalls("output_firstCalcDepthCC.txt", px1, py1, px2, py2, 1, vdepth);
    }

	// printNumberOfCalls("output_firstCalcDepthResult.txt", px1, py1, px2, py2, 1, vdepth);

    return vdepth * (MAP_AREA/matrix_size); // multiply by the size of the block
}

static void string_to_mac_address(const char *str, u8 *addr)
{
	int a[ETH_ALEN];

	sscanf(str, "%x:%x:%x:%x:%x:%x",
	       &a[0], &a[1], &a[2], &a[3], &a[4], &a[5]);

	addr[0] = (u8) a[0];
	addr[1] = (u8) a[1];
	addr[2] = (u8) a[2];
	addr[3] = (u8) a[3];
	addr[4] = (u8) a[4];
	addr[5] = (u8) a[5];
}

static int get_link_snr_default(struct wmediumd *ctx, struct station *sender,
				 struct station *receiver)
{
	return SNR_DEFAULT;
}

static int get_link_snr_from_snr_matrix(struct wmediumd *ctx,
					struct station *sender,
					struct station *receiver)
{
	return ctx->snr_matrix[sender->index * ctx->num_stas + receiver->index];
}

static double _get_error_prob_from_snr(struct wmediumd *ctx, double snr,
					   unsigned int rate_idx, u32 freq,
					   int frame_len,
				       struct station *src, struct station *dst)
{
	return get_error_prob_from_snr(snr, rate_idx, freq, frame_len);
}

static double get_error_prob_from_matrix(struct wmediumd *ctx, double snr,
					 unsigned int rate_idx, u32 freq,
					 int frame_len, struct station *src,
					 struct station *dst)
{
	if (dst == NULL) // dst is multicast. returned value will not be used.
		return 0.0;

	return ctx->error_prob_matrix[ctx->num_stas * src->index + dst->index];
}

int use_fixed_random_value(struct wmediumd *ctx)
{
	return ctx->error_prob_matrix != NULL || ctx->station_err_matrix != NULL;
}

#define FREQ_1CH (2.412e9)		// [Hz]
#define SPEED_LIGHT (2.99792458e8)	// [meter/sec]
/*
 * Calculate path loss based on a free-space path loss
 *
 * This function returns path loss [dBm].
 */
static int calc_path_loss_free_space(void *model_param,
			  struct station *dst, struct station *src)
{
	struct free_space_model_param *param;
	double PL, d, denominator, numerator, lambda;
	double f = src->freq * pow(10,6);

	if (f < 0.1)
		f = FREQ_1CH;

	param = model_param;

	d = sqrt((src->x - dst->x) * (src->x - dst->x) +
			 (src->y - dst->y) * (src->y - dst->y) +
			 (src->z - dst->z) * (src->z - dst->z));

	/*
	 * Calculate PL0 with Free-space path loss in decibels
	 *
	 * 20 * log10 * (4 * M_PI * d * f / c)
	 *   d: distance [meter]
	 *   f: frequency [Hz]
	 *   c: speed of light in a vacuum [meter/second]
	 *
	 * https://en.wikipedia.org/wiki/Free-space_path_loss
	 */
	lambda = SPEED_LIGHT / f;
	denominator = pow(lambda, 2);
	numerator = pow((4.0 * M_PI * d), 2) * param->sL;
	PL = 10.0 * log10(numerator / denominator);
	return PL;
}
/*
 * Calculate path loss based on a log distance model
 *
 * This function returns path loss [dBm].
 */
static int calc_path_loss_log_distance(void *model_param,
			  struct station *dst, struct station *src)
{
	struct log_distance_model_param *param;
	double PL, PL0, d;
	double f = src->freq * pow(10,6);

	if (f < 0.1)
		f = FREQ_1CH;

	param = model_param;

	d = sqrt((src->x - dst->x) * (src->x - dst->x) +
		 (src->y - dst->y) * (src->y - dst->y) +
		 (src->z - dst->z) * (src->z - dst->z));

	/*
	 * Calculate PL0 with Free-space path loss in decibels
	 *
	 * 20 * log10 * (4 * M_PI * d * f / c)
	 *   d: distance [meter]
	 *   f: frequency [Hz]
	 *   c: speed of light in a vacuum [meter/second]
	 *
	 * https://en.wikipedia.org/wiki/Free-space_path_loss
	 */
	PL0 = 20.0 * log10(4.0 * M_PI * 1.0 * f / SPEED_LIGHT);

	/*
	 * Calculate signal strength with Log-distance path loss model
	 * https://en.wikipedia.org/wiki/Log-distance_path_loss_model
	 */
	PL = PL0 + 10.0 * param->path_loss_exponent * log10(d) + param->Xg;
	return PL;
}
static int calc_path_loss_weissberger(void *model_param, 
				struct station *dst, struct station *src)
{
	struct weissberger_model_param *param;
	int log_distance_path_loss;
	double PL = 0;
	double f = src->freq / 1000;

	param = model_param;

	log_distance_path_loss = calc_path_loss_log_distance(&param->logd_param, dst, src);

	int x1 = (int)((src->x) / 10.0);
	int y1 = 149 - (int)(src->y / 10.0);
	int x2 = (int)((dst->x) / 10.0);
	int y2 = 149 -(int)(dst->y / 10.0);

	// static int call_count = 0;
	// call_count++;

	// if (call_count == 1)
	// 	printNumberOfCalls("outputTest.txt", x1, y1, x2, y2, call_count, vegetation_matrix[1][15]);

	// printNumberOfCalls("output4.txt", x1, y1, x2, y2, call_count, -1);

	int vegetation_depth = calculate_vegetation_depth(x1, y1, x2, y2, param->matrix_size);

	if (0 < vegetation_depth && vegetation_depth <= 14) {
		PL = 0.45 * vegetation_depth * pow(f, 0.284);
	}
	else if (14 < vegetation_depth && vegetation_depth <= 400) {
		PL = 1.33 * pow(vegetation_depth, 0.588) * pow(f, 0.284);
		// printVegetationDepth("output5.txt", vegetation_depth, PL);
	}
	else if (vegetation_depth > 400) { //assume signal is definitely loss after 400m of vegetation
		PL = 10000;
	}

	// printWeissbergerOutput("output_weissberger.txt", log_distance_path_loss, PL, f, src, dst, vegetation_depth);

	// if (call_count % 3 == 0) { //it generates A LOT of logs, so only get a third of it
	// 	FILE *log_file = fopen("vegetation_report.log", "a");
	// if (log_file == NULL) {
	// 	perror("Error opening log file");
	// 	}
	// 	else {
	// // Write the vegetation depth and PL to the log file
	// if (PL < 200000)
	// 	fprintf(log_file, "From: %d (%.1f,%.1f), To %d (%.1f,%.1f), Vegetation Depth: %d, Weissberger: %.2f, Total Path Loss: %d\n",
    //         src->index, src->x, src->y, dst->index, dst->x, dst->y, vegetation_depth, PL, (int)(PL + log_distance_path_loss));

	// // Close the log file
	// fclose(log_file);
	// 	}
	// }

	return PL + log_distance_path_loss;
}
/*
 * Calculate path loss based on a itu model
 *
 * This function returns path loss [dBm].
 */
static int calc_path_loss_itu(void *model_param,
			  struct station *dst, struct station *src)
{
	struct itu_model_param *param;
	double PL, d;
	double f = src->freq;
	int N=28, pL;

	if (f < 0.1)
		f = FREQ_1CH;

	param = model_param;
	pL = param->pL;

	d = sqrt((src->x - dst->x) * (src->x - dst->x) +
			 (src->y - dst->y) * (src->y - dst->y) +
			 (src->z - dst->z) * (src->z - dst->z));

	if (d>16)
		N=38;
	if (pL!=0)
		N=pL;
	/*
	 * Calculate signal strength with ITU path loss model
	 * Power Loss Coefficient Based on the Paper
     * Site-Specific Validation of ITU Indoor Path Loss Model at 2.4 GHz
     * from Theofilos Chrysikos, Giannis Georgopoulos and Stavros Kotsopoulos
     * LF: floor penetration loss factor
     * nFLOORS: number of floors
	 */

	PL = 20.0 * log10(f) + N * log10(d) + param->lF * param->nFLOORS - 28;
	return PL;
}
/*
 * Calculate path loss based on a log-normal shadowing model
 *
 * This function returns path loss [dBm].
 */
static int calc_path_loss_log_normal_shadowing(void *model_param,
			  struct station *dst, struct station *src)
{
	struct log_normal_shadowing_model_param *param;
	double PL, PL0, d;
	double f = src->freq * pow(10,6);
	double gRandom = src->gRandom;

	if (f < 0.1)
		f = FREQ_1CH;

	param = model_param;

	d = sqrt((src->x - dst->x) * (src->x - dst->x) +
		 (src->y - dst->y) * (src->y - dst->y) +
		 (src->z - dst->z) * (src->z - dst->z));

	/*
	 * Calculate PL0 with Free-space path loss in decibels
	 *
	 * 20 * log10 * (4 * M_PI * d * f / c)
	 *   d: distance [meter]
	 *   f: frequency [Hz]
	 *   c: speed of light in a vacuum [meter/second]
	 *
	 * https://en.wikipedia.org/wiki/Free-space_path_loss
	 */
	PL0 = 20.0 * log10(4.0 * M_PI * 1.0 * f / SPEED_LIGHT);

	/*
	 * Calculate signal strength with Log-distance path loss model + gRandom (Gaussian random variable)
	 * https://en.wikipedia.org/wiki/Log-distance_path_loss_model
	 */
	PL = PL0 + 10.0 * param->path_loss_exponent * log10(d) - gRandom;
	return PL;
}
/*
 * Calculate path loss based on a two ray ground model
 *
 * This function returns path loss [dBm].
 */
static int calc_path_loss_two_ray_ground(void *model_param,
			  struct station *dst, struct station *src)
{
	//struct two_ray_ground_model_param *param;
        double PL, d;
        double f = 20 * 1000000; //frequency in Hz
        double lambda = SPEED_LIGHT / f;
        int ht = 1;
        int hr = 1;
        double dCross = (4 * M_PI * ht * hr) / (lambda / 1000);

        //param = model_param;
        d = sqrt((src->x - dst->x) * (src->x - dst->x) +
                         (src->y - dst->y) * (src->y - dst->y) +
                         (src->z - dst->z) * (src->z - dst->z));

        if (d < dCross){
            PL = calc_path_loss_free_space(model_param, dst, src);
            return PL;
        }
        else{
            double numerator = src->tx_power * src->gain * dst->gain * pow(ht, 2) * pow(hr, 2);
            double denominator = pow(d, 4);
            PL = (numerator / denominator);
            return PL;
        }
}

/* Existing link is from from -> to; copy to other dir */
static void mirror_link(struct wmediumd *ctx, int from, int to)
{
	ctx->snr_matrix[ctx->num_stas * to + from] =
		ctx->snr_matrix[ctx->num_stas * from + to];

	if (ctx->error_prob_matrix) {
		ctx->error_prob_matrix[ctx->num_stas * to + from] =
			ctx->error_prob_matrix[ctx->num_stas * from + to];
    }
}

static void recalc_path_loss(struct wmediumd *ctx)
{
	int start, end, path_loss, gains, txpower, signal;

	for (start = 0; start < ctx->num_stas; start++) {
		for (end = 0; end < ctx->num_stas; end++) {
			if (start == end)
				continue;
			txpower = ctx->sta_array[start]->tx_power;
			if (ctx->sta_array[end]->isap == 1)
				txpower = ctx->sta_array[end]->tx_power;

			path_loss = ctx->calc_path_loss(ctx->path_loss_param,
				ctx->sta_array[end], ctx->sta_array[start]);
			gains = txpower + ctx->sta_array[start]->gain + ctx->sta_array[end]->gain;
			signal = gains - path_loss - ctx->noise_threshold;
            ctx->snr_matrix[ctx->num_stas * start + end] = signal;
            ctx->snr_matrix[ctx->num_stas * end + start] = signal;
	}
    }
}

static void move_stations_to_direction(struct wmediumd *ctx)
{
	struct station *station;
	struct timespec now;

	static int call_count = 0;
	call_count++;

	clock_gettime(CLOCK_MONOTONIC, &now);
	if (!timespec_before(&ctx->next_move, &now))
		return;

	list_for_each_entry(station, &ctx->stations, list) {
		station->x += station->dir_x;
		station->y += station->dir_y;
	}
	recalc_path_loss(ctx);

	clock_gettime(CLOCK_MONOTONIC, &ctx->next_move);
	ctx->next_move.tv_sec += MOVE_INTERVAL;
}

static void move_stations_donothing(struct wmediumd *ctx)
{
}

static int parse_path_loss(struct wmediumd *ctx, config_t *cf)
{
	struct station *station;
	const config_setting_t *positions, *position;
	const config_setting_t *directions, *direction;
	const config_setting_t *tx_powers, *model;
	const config_setting_t *isnodeaps;
	const char *path_loss_model_name;

	positions = config_lookup(cf, "model.positions");
	if (!positions) {
		w_flogf(ctx, LOG_ERR, stderr,
			"No positions found in model\n");
		return -EINVAL;
	}
	if (config_setting_length(positions) != ctx->num_stas) {
		w_flogf(ctx, LOG_ERR, stderr,
			"Specify %d positions\n", ctx->num_stas);
		return -EINVAL;
	}

	directions = config_lookup(cf, "model.directions");
	if (directions) {
		if (config_setting_length(directions) != ctx->num_stas) {
			w_flogf(ctx, LOG_ERR, stderr,
				"Specify %d directions\n", ctx->num_stas);
			return -EINVAL;
		}
		ctx->move_stations = move_stations_to_direction;
	}

	tx_powers = config_lookup(cf, "model.tx_powers");
	if (!tx_powers) {
		w_flogf(ctx, LOG_ERR, stderr,
			"No tx_powers found in model\n");
		return -EINVAL;
	}
	if (config_setting_length(tx_powers) != ctx->num_stas) {
		w_flogf(ctx, LOG_ERR, stderr,
			"Specify %d tx_powers\n", ctx->num_stas);
		return -EINVAL;
	}

	isnodeaps = config_lookup(cf, "model.isnodeaps");

	model = config_lookup(cf, "model");
	if (config_setting_lookup_string(model, "model_name",
		&path_loss_model_name) != CONFIG_TRUE) {
		w_flogf(ctx, LOG_ERR, stderr, "Specify model_name\n");
		return -EINVAL;
	}
	if (strncmp(path_loss_model_name, "log_distance",
		    sizeof("log_distance")) == 0) {
		struct log_distance_model_param *param;
		ctx->calc_path_loss = calc_path_loss_log_distance;
		param = malloc(sizeof(*param));
		if (!param) {
			w_flogf(ctx, LOG_ERR, stderr,
				"Out of memory(path_loss_param)\n");
			return -EINVAL;
		}

		if (config_setting_lookup_float(model, "path_loss_exp",
			&param->path_loss_exponent) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr,
				"path_loss_exponent not found\n");
			return -EINVAL;
		}

		if (config_setting_lookup_float(model, "xg",
			&param->Xg) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr, "xg not found\n");
			return -EINVAL;
		}
		ctx->path_loss_param = param;
	}
	else if (strncmp(path_loss_model_name, "free_space",
			sizeof("free_space")) == 0) {
		struct free_space_model_param *param;
		ctx->calc_path_loss = calc_path_loss_free_space;
		param = malloc(sizeof(*param));
		if (!param) {
			w_flogf(ctx, LOG_ERR, stderr,
				"Out of memory(path_loss_param)\n");
			return -EINVAL;
		}

		if (config_setting_lookup_int(model, "sL",
			&param->sL) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr,
				"system loss not found\n");
			return -EINVAL;
		}
		ctx->path_loss_param = param;
	}
	else if (strncmp(path_loss_model_name, "log_normal_shadowing",
			sizeof("log_normal_shadowing")) == 0) {
		struct log_normal_shadowing_model_param *param;
		ctx->calc_path_loss = calc_path_loss_log_normal_shadowing;
		param = malloc(sizeof(*param));
		if (!param) {
			w_flogf(ctx, LOG_ERR, stderr,
				"Out of memory(path_loss_param)\n");
			return -EINVAL;
		}

		if (config_setting_lookup_int(model, "sL",
			&param->sL) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr,
				"system loss not found\n");
			return -EINVAL;
		}

		if (config_setting_lookup_float(model, "path_loss_exp",
			&param->path_loss_exponent) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr,
				"path_loss_exponent not found\n");
			return -EINVAL;
		}
		ctx->path_loss_param = param;
	}
	else if (strncmp(path_loss_model_name, "two_ray_ground",
				sizeof("two_ray_ground")) == 0) {
		struct two_ray_ground_model_param *param;
		ctx->calc_path_loss = calc_path_loss_two_ray_ground;
		param = malloc(sizeof(*param));
		if (!param) {
			w_flogf(ctx, LOG_ERR, stderr,
				"Out of memory(path_loss_param)\n");
			return -EINVAL;
		}

		if (config_setting_lookup_int(model, "sL",
			&param->sL) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr,
				"system loss not found\n");
			return -EINVAL;
		}
		ctx->path_loss_param = param;
	}
	else if (strncmp(path_loss_model_name, "itu",
			sizeof("itu")) == 0) {
		struct itu_model_param *param;
		ctx->calc_path_loss = calc_path_loss_itu;
		param = malloc(sizeof(*param));
		if (!param) {
			w_flogf(ctx, LOG_ERR, stderr,
				"Out of memory(path_loss_param)\n");
			return -EINVAL;
		}

		if (config_setting_lookup_int(model, "nFLOORS",
			&param->nFLOORS) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr,
				"nFLOORS not found\n");
			return -EINVAL;
		}

		if (config_setting_lookup_int(model, "lF",
			&param->lF) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr, "LF not found\n");
			return -EINVAL;
		}

		if (config_setting_lookup_int(model, "pL",
			&param->pL) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr, "PL not found\n");
			return -EINVAL;
		}
		ctx->path_loss_param = param;
	}
	else if(strncmp(path_loss_model_name, "weissberger",
			sizeof("weissberger")) == 0) {
		struct weissberger_model_param *param;
		ctx->calc_path_loss = calc_path_loss_weissberger;
		param = malloc(sizeof(*param));
		if (!param) {
			w_flogf(ctx, LOG_ERR, stderr,
				"Out of memory(path_loss_param)\n");
			return -EINVAL;
		}

		if (config_setting_lookup_float(model, "path_loss_exp",
			&param->logd_param.path_loss_exponent) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr,
				"path_loss_exponent not found\n");
			return -EINVAL;
		}

		if (config_setting_lookup_float(model, "xg",
			&param->logd_param.Xg) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr, "xg not found\n");
			return -EINVAL;
		}

		if (config_setting_lookup_string(model, "vegetation_matrix_file",
			&param->vegetation_matrix_file) != CONFIG_TRUE) {
			// printErrorMessage("config_lookup_error.txt", "vegetation_matrix_file not found", 1);
			w_flogf(ctx, LOG_ERR, stderr, "vegetation_matrix not found\n");
			return -EINVAL;
		}

		if (config_setting_lookup_int(model, "matrix_size",
			&param->matrix_size) != CONFIG_TRUE) {
			w_flogf(ctx, LOG_ERR, stderr, "matrix_size not found\n");
			return -EINVAL;
		} 

		// printErrorMessage("before_load.txt", param->vegetation_matrix_file, param->matrix_size);

		load_vegetation_matrix(param->vegetation_matrix_file, param->matrix_size);

		ctx->path_loss_param = param;
	}
	else {
		w_flogf(ctx, LOG_ERR, stderr, "No path loss model found\n");
		return -EINVAL;
	}

	list_for_each_entry(station, &ctx->stations, list) {
		position = config_setting_get_elem(positions, station->index);
		if (config_setting_length(position) != 3) {
			w_flogf(ctx, LOG_ERR, stderr,
				"Invalid position: expected (double,double,double)\n");
			return -EINVAL;
		}
		station->x = config_setting_get_float_elem(position, 0);
		station->y = config_setting_get_float_elem(position, 1);
		station->z = config_setting_get_float_elem(position, 2);

		if (directions) {
			direction = config_setting_get_elem(directions,
				station->index);
			if (config_setting_length(direction) != 2) {
				w_flogf(ctx, LOG_ERR, stderr,
					"Invalid direction: expected (double,double)\n");
				return -EINVAL;
			}
			station->dir_x = config_setting_get_float_elem(
				direction, 0);
			station->dir_y = config_setting_get_float_elem(
				direction, 1);
		}

		station->tx_power = config_setting_get_float_elem(
			tx_powers, station->index);
		if (isnodeaps) {
			station->isap = config_setting_get_int_elem(
				isnodeaps, station->index);
		}
	}

	recalc_path_loss(ctx);

	return 0;
}

static double pseudo_normal_distribution(void)
{
	int i;
	double normal = -6.0;

	for (i = 0; i < 12; i++)
		normal += drand48();

	return normal;
}

static int _get_fading_signal(struct wmediumd *ctx)
{
	return ctx->fading_coefficient * pseudo_normal_distribution();
}

static int get_no_fading_signal(struct wmediumd *ctx)
{
	return 0;
}

/*
 *	Loads a config file into memory
 */
int load_config(struct wmediumd *ctx, const char *file, const char *per_file, bool full_dynamic)
{
	config_t cfg, *cf;
	const config_setting_t *ids, *links, *model_type;
	const config_setting_t *error_probs = NULL, *error_prob;
	const config_setting_t *enable_interference;
	const config_setting_t *fading_coefficient, *noise_threshold, *default_prob;
	const config_setting_t *mediums, *medium_data,*interface_data, *medium_detection;
	int count_ids, count_mediums, count_interfaces, station_id, i, j;
	int start, end, snr;
	struct station *station;
	const char *model_type_str;
	float default_prob_value = 0.0;
	bool *link_map = NULL;

	if (full_dynamic) {
		ctx->sta_array = malloc(0);
		ctx->num_stas = 0;
		ctx->intf = NULL;
		ctx->get_fading_signal = get_no_fading_signal;
		ctx->fading_coefficient = 0;
		ctx->move_stations = move_stations_donothing;
		ctx->snr_matrix = malloc(0);
		ctx->per_matrix = NULL;
		ctx->per_matrix_row_num = 0;
		ctx->error_prob_matrix = NULL;
		ctx->get_link_snr = get_link_snr_default;
		ctx->station_err_matrix = malloc(0);
		return 0;
	}
	ctx->station_err_matrix = NULL;

	/*initialize the config file*/
	cf = &cfg;
	config_init(cf);

	/*read the file*/
	if (!config_read_file(cf, file)) {
		w_logf(ctx, LOG_ERR, "Error loading file %s at line:%d, reason: %s\n",
				file,
				config_error_line(cf),
				config_error_text(cf));
		config_destroy(cf);
		return -EIO;
	}

	ids = config_lookup(cf, "ifaces.ids");
	if (!ids) {
		w_logf(ctx, LOG_ERR, "ids not found in config file\n");
		return -EIO;
	}
	count_ids = config_setting_length(ids);

	w_logf(ctx, LOG_NOTICE, "#_if = %d\n", count_ids);

	/* Fill the mac_addr */
	ctx->sta_array = malloc(sizeof(struct station *) * count_ids);
	if (!ctx->sta_array) {
		w_flogf(ctx, LOG_ERR, stderr, "Out of memory(sta_array)!\n");
		return -ENOMEM;
	}
	for (i = 0; i < count_ids; i++) {
		u8 addr[ETH_ALEN];
		const char *str =  config_setting_get_string_elem(ids, i);
		string_to_mac_address(str, addr);

		station = malloc(sizeof(*station));
		if (!station) {
			w_flogf(ctx, LOG_ERR, stderr, "Out of memory!\n");
			return -ENOMEM;
		}
		station->index = i;
		memcpy(station->addr, addr, ETH_ALEN);
		memcpy(station->hwaddr, addr, ETH_ALEN);
		station->tx_power = SNR_DEFAULT;
		station->gain = GAIN_DEFAULT;
		//station->height = HEIGHT_DEFAULT;
		station->gRandom = GAUSS_RANDOM_DEFAULT;
		station->isap = AP_DEFAULT;
		station->medium_id = MEDIUM_ID_DEFAULT;
		station_init_queues(station);
		list_add_tail(&station->list, &ctx->stations);
		ctx->sta_array[i] = station;

		w_logf(ctx, LOG_NOTICE, "Added station %d: " MAC_FMT "\n", i, MAC_ARGS(addr));
	}
	ctx->num_stas = count_ids;

	enable_interference = config_lookup(cf, "ifaces.enable_interference");
	if (enable_interference &&
	    config_setting_get_bool(enable_interference)) {
		ctx->intf = calloc(ctx->num_stas * ctx->num_stas,
				   sizeof(struct intf_info));
		if (!ctx->intf) {
			w_flogf(ctx, LOG_ERR, stderr, "Out of memory(intf)\n");
			return -ENOMEM;
		}
		for (i = 0; i < ctx->num_stas; i++)
			for (j = 0; j < ctx->num_stas; j++)
				ctx->intf[i * ctx->num_stas + j].signal = -200;
	} else {
		ctx->intf = NULL;
	}

    noise_threshold =
		config_lookup(cf, "model.noise_threshold");
	if (noise_threshold &&
	    config_setting_get_int(noise_threshold) < 0) {
		ctx->get_fading_signal = _get_fading_signal;
		ctx->noise_threshold =
			config_setting_get_int(noise_threshold);
	} else {
		ctx->noise_threshold = NOISE_LEVEL;
	}


	fading_coefficient =
		config_lookup(cf, "model.fading_coefficient");
	if (fading_coefficient &&
	    config_setting_get_int(fading_coefficient) > 0) {
		ctx->get_fading_signal = _get_fading_signal;
		ctx->fading_coefficient =
			config_setting_get_int(fading_coefficient);
	} else {
		ctx->get_fading_signal = get_no_fading_signal;
		ctx->fading_coefficient = 0;
	}

	ctx->move_stations = move_stations_donothing;

	/* create link quality matrix */
	ctx->snr_matrix = calloc(sizeof(int), count_ids * count_ids);
	if (!ctx->snr_matrix) {
		w_flogf(ctx, LOG_ERR, stderr, "Out of memory!\n");
		return -ENOMEM;
	}
	/* set default snrs */
	for (i = 0; i < count_ids * count_ids; i++)
		ctx->snr_matrix[i] = SNR_DEFAULT;

	links = config_lookup(cf, "ifaces.links");
	if (!links) {
		model_type = config_lookup(cf, "model.type");
		if (model_type) {
			model_type_str = config_setting_get_string(model_type);
			if (memcmp("snr", model_type_str, strlen("snr")) == 0) {
				links = config_lookup(cf, "model.links");
			} else if (memcmp("prob", model_type_str,
				strlen("prob")) == 0) {
				error_probs = config_lookup(cf, "model.links");
			} else if (memcmp("path_loss", model_type_str,
				strlen("path_loss")) == 0) {
				/* calculate signal from positions */
				if (parse_path_loss(ctx, cf))
					goto fail;
			}
		}
	}
    mediums = config_lookup(cf, "ifaces.medium_array");
    if (mediums) {
        count_mediums = config_setting_length(mediums);
        for (i = 0; i < count_mediums; i++) {
            medium_data = config_setting_get_elem(mediums, i);
            count_interfaces = config_setting_length(medium_data);
            for (j = 0; j < count_interfaces; j++) {
                interface_data = config_setting_get_elem(medium_data, j);
                station_id = config_setting_get_int(interface_data);
                ctx->sta_array[station_id]->medium_id = i+1;
            }
        }
    }
    medium_detection = config_lookup(cf, "ifaces.enable_medium_detection");
    if (medium_detection) {
        ctx->enable_medium_detection =config_setting_get_bool(enable_interference);
    }else{
        ctx->enable_medium_detection = ENABLE_MEDIUM_DETECTION;
    }

	if (per_file && error_probs) {
		w_flogf(ctx, LOG_ERR, stderr,
			"per_file and error_probs could not be used at the same time\n");
		goto fail;
	}

	ctx->get_link_snr = get_link_snr_from_snr_matrix;
	ctx->get_error_prob = _get_error_prob_from_snr;

	ctx->per_matrix = NULL;
	ctx->per_matrix_row_num = 0;
	if (per_file && read_per_file(ctx, per_file))
		goto fail;

	ctx->error_prob_matrix = NULL;
	if (error_probs) {
		ctx->error_prob_matrix = calloc(sizeof(double),
						count_ids * count_ids);
		if (!ctx->error_prob_matrix) {
			w_flogf(ctx, LOG_ERR, stderr,
				"Out of memory(error_prob_matrix)\n");
			goto fail;
		}

		ctx->get_link_snr = get_link_snr_default;
		ctx->get_error_prob = get_error_prob_from_matrix;

		default_prob = config_lookup(cf, "model.default_prob");
		if (default_prob) {
			default_prob_value = config_setting_get_float(
				default_prob);
			if (default_prob_value < 0.0 ||
			    default_prob_value > 1.0) {
					w_flogf(ctx, LOG_ERR, stderr,
						"model.default_prob should be in [0.0, 1.0]\n");
					goto fail;
			}
		}
	}

	link_map = calloc(sizeof(bool), count_ids * count_ids);
	if (!link_map) {
		w_flogf(ctx, LOG_ERR, stderr, "Out of memory\n");
		goto fail;
	}

	/* read snr values */
	for (i = 0; links && i < config_setting_length(links); i++) {
		config_setting_t *link;

		link = config_setting_get_elem(links, i);
		if (config_setting_length(link) != 3) {
			w_flogf(ctx, LOG_ERR, stderr, "Invalid link: expected (int,int,int)\n");
			goto fail;
		}
		start = config_setting_get_int_elem(link, 0);
		end = config_setting_get_int_elem(link, 1);
		snr = config_setting_get_int_elem(link, 2);

		if (start < 0 || start >= ctx->num_stas ||
		    end < 0 || end >= ctx->num_stas) {
				w_flogf(ctx, LOG_ERR, stderr, "Invalid link [%d,%d,%d]: index out of range\n",
						start, end, snr);
				goto fail;
		}
		ctx->snr_matrix[ctx->num_stas * start + end] = snr;
		link_map[ctx->num_stas * start + end] = true;
	}

	/* initialize with default_prob */
	for (start = 0; error_probs && start < ctx->num_stas; start++)
		for (end = start + 1; end < ctx->num_stas; end++) {
			ctx->error_prob_matrix[ctx->num_stas *
				start + end] =
			ctx->error_prob_matrix[ctx->num_stas *
				end + start] = default_prob_value;
		}

	/* read error probabilities */
	for (i = 0; error_probs &&
	     i < config_setting_length(error_probs); i++) {
			float error_prob_value;

			error_prob = config_setting_get_elem(error_probs, i);
			if (config_setting_length(error_prob) != 3) {
				w_flogf(ctx, LOG_ERR, stderr, "Invalid error probability: expected (int,int,float)\n");
				goto fail;
			}

			start = config_setting_get_int_elem(error_prob, 0);
			end = config_setting_get_int_elem(error_prob, 1);
			error_prob_value = config_setting_get_float_elem(error_prob, 2);

			if (start < 0 || start >= ctx->num_stas ||
				end < 0 || end >= ctx->num_stas ||
				error_prob_value < 0.0 || error_prob_value > 1.0) {
					w_flogf(ctx, LOG_ERR, stderr, "Invalid error probability [%d,%d,%f]\n",
						start, end, error_prob_value);
					goto fail;
			}

			ctx->error_prob_matrix[ctx->num_stas * start + end] =
				error_prob_value;
			link_map[ctx->num_stas * start + end] = true;
	}

	/*
	 * If any links are specified in only one direction, mirror them,
	 * making them symmetric.  If specified in both directions they
	 * can be asymmetric.
	 */
	for (i = 0; i < ctx->num_stas; i++) {
		for (j = 0; j < ctx->num_stas; j++) {
			if (i == j)
				continue;

			if (link_map[ctx->num_stas * i + j] &&
			    !link_map[ctx->num_stas * j + i]) {
				mirror_link(ctx, i, j);
			}
		}
	}

    free(link_map);
	config_destroy(cf);
	return 0;

fail:
	free(ctx->snr_matrix);
	free(ctx->error_prob_matrix);
	config_destroy(cf);
	return -EINVAL;
}
