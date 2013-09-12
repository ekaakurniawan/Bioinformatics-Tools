// Copyright (C) 2013 by Eka A. Kurniawan
// eka.a.kurniawan(ta)gmail(tod)com
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the
// Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

// -----------------------------------------------------------------------------
// Sequential Hierarchical Agglomerative Clustering Using Single-Linkage Method
// -----------------------------------------------------------------------------
// Normal compilation:
//    gcc -lm main.c
// Optimized compilation:
//    gcc -lm -O3 main.c
// To run:
//    ./a.out dataFile1.data 90000 90000.txt

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX_FLOAT 1000
#define MIN_CLUSTERS 9
#define E_DIST(a, b) sqrtf(pow(a.x - b.x, 2) + pow(a.y - b.y, 2))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

typedef struct dFloat {
    float x;
    float y;
} point2D;

int main (int argc, char *argv[]) {
    FILE *fin_p, *fout_p;
    int ttl_points;
    
    point2D *points;
    float **dist_m;
    float *min_row_dist;
    int *min_row_loc;
    int *clusters;

    float dist, min_dist;
    int min_loc, min_loc_i, min_loc_j;
    int ci, cj;
    
    int loop, i, j;
    time_t t_start, t_tic;
    double t_diff;

    t_start = time(NULL);

    // check and acquire parameters
    if (argc != 4) {
        fprintf(stdout, "Use: ./a.out input_file total_points output_file\n");
        fprintf(stdout, "Example: ./a.out dataFile1.data 90000 out.txt\n");
        exit(EXIT_FAILURE);
    }
    fin_p = fopen(argv[1], "r");
    fout_p = fopen(argv[3], "w");
    ttl_points = atoi(argv[2]);

    t_tic = time(NULL);
    // allocate memory for points
    points = (point2D *)malloc(ttl_points * sizeof(point2D));
    // get points from input file
    for (i = 0; i < ttl_points; i++) {
        fscanf(fin_p, "%f %f", &points[i].x, &points[i].y);
    }
    t_diff = difftime(time(NULL), t_tic);
    fprintf(stdout, "TIME - Allocate Points                     : %lf seconds\n", t_diff);
    fprintf(fout_p, "TIME - Allocate Points                     : %lf seconds\n", t_diff);

    t_tic = time(NULL);
    // allocate memory for distance matrix
    dist_m = (float **)malloc(ttl_points * sizeof(float *));
    for (i = 0; i < ttl_points; i++) {
        dist_m[i] = (float *)malloc(ttl_points * sizeof(float));
    }
    // calculate distances and store them at distance matrix
    for (i = 0; i < ttl_points; i++) {
        for (j = i + 1; j < ttl_points; j++) {
            dist_m[i][j] = E_DIST(points[i], points[j]);
        }
    }
    t_diff = difftime(time(NULL), t_tic);
    fprintf(stdout, "TIME - Allocate And Calculate Distance     : %lf seconds\n", t_diff);
    fprintf(fout_p, "TIME - Allocate And Calculate Distance     : %lf seconds\n", t_diff);

    // allocate memory for minimum row
    min_row_dist = (float *)malloc(ttl_points * sizeof(float));
    min_row_loc = (int *)malloc(ttl_points * sizeof(int));
    
    // allocate memory for cluster allocation
    clusters = (int *)malloc(ttl_points * sizeof(int));
    for (i = 0; i < ttl_points; i++) {
        clusters[i] = i;
    }

    t_tic = time(NULL);
    // get minimum distance for each row
    for (i = 0; i < ttl_points; i++) {
        min_dist = MAX_FLOAT;
        for (j = i + 1; j < ttl_points; j++) {
            dist = dist_m[i][j];
            if (dist < min_dist) {
                min_dist = dist;
                min_loc = j;
            }
        }
        min_row_dist[i] = min_dist;
        min_row_loc[i] = min_loc;
    }

    t_diff = difftime(time(NULL), t_tic);
    fprintf(stdout, "TIME - Get Minimum Distance Per Row        : %lf seconds\n", t_diff);
    fprintf(fout_p, "TIME - Get Minimum Distance Per Row        : %lf seconds\n", t_diff);

    t_tic = time(NULL);
    for (loop = 0; loop < ttl_points - MIN_CLUSTERS; loop++) {
        if (!(loop % 1000)) fprintf(stdout, "loop: %d\n", loop);

        // get global minimum distance from minimum row distance
        min_dist = MAX_FLOAT;
        for (i = 0; i < ttl_points; i++) {
            dist = min_row_dist[i];
            if (dist < min_dist) {
                min_dist = dist;
                min_loc = i;
            }
        }

        // update distance matrix
        min_loc_i = min_loc;
        min_loc_j = min_row_loc[min_loc_i];
        for (i = 0; i < min_loc_i; i++) {
            dist_m[i][min_loc_i] = MIN(dist_m[i][min_loc_i], dist_m[i][min_loc_j]);
        }
        for (j = min_loc_i + 1; j < min_loc_j; j++) {
            dist_m[min_loc_i][j] = MIN(dist_m[min_loc_i][j], dist_m[j][min_loc_j]);
        }
        for (j = min_loc_j + 1; j < ttl_points; j++) {
            dist_m[min_loc_i][j] = MIN(dist_m[min_loc_i][j], dist_m[min_loc_j][j]);
        }
        
        // remove one of combined cluster
        for (i = 0; i < min_loc_j; i++) {
            dist_m[i][min_loc_j] = MAX_FLOAT;
        }
        for (j = min_loc_j + 1; j < ttl_points; j++) {
            dist_m[min_loc_j][j] = MAX_FLOAT;
        }

        // update minimum distance and location for new row distance
        for (i = 0; i < min_loc_i; i++) {
            if (min_row_loc[i] == min_loc_j) {
                min_row_loc[i] = min_loc_i;
            }
        }

        for (i = min_loc_i; i < min_loc_j; i++) {
            if (min_row_loc[i] == min_loc_j) {
                min_dist = MAX_FLOAT;
                for (j = i + 1; j < ttl_points; j++) {
                    dist = dist_m[i][j];
                    if (dist < min_dist) {
                        min_dist = dist;
                        min_loc = j;
                    }
                }
                min_row_dist[i] = min_dist;
                min_row_loc[i] = min_loc;
            }
        }

        min_row_dist[min_loc_j] = MAX_FLOAT;

        // update clusters
        ci = clusters[min_loc_i];
        cj = clusters[min_loc_j];
        for (i = 0; i < ttl_points; i++) {
            if (clusters[i] == cj) {
                clusters[i] = ci;
            }
        }
    }
    t_diff = difftime(time(NULL), t_tic);
    fprintf(stdout, "TIME - Clustering                          : %lf seconds\n", t_diff);
    fprintf(fout_p, "TIME - Clustering                          : %lf seconds\n", t_diff);

    t_diff = difftime(time(NULL), t_start);
    fprintf(stdout, "TIME - Total                               : %lf seconds\n", t_diff);
    fprintf(fout_p, "TIME - Total                               : %lf seconds\n", t_diff);

    // store final cluster at output file
    fprintf(fout_p, "\nClusters :\n");
    for (i = 0; i < ttl_points; i++) {
        fprintf(fout_p, " %d", clusters[i]);
    }

    fprintf(stdout, "Done!\n");
    fclose(fin_p);
    fclose(fout_p);
    free(points);
    for (i = 0; i < ttl_points; i++) {
        free(dist_m[i]);
    }
    free(dist_m);
    free(min_row_dist);
    free(min_row_loc);
    free(clusters);
    return 0;
}
