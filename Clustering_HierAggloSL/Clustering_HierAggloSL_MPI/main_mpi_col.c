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
// Parallel Hierarchical Agglomerative Clustering Using Single-Linkage Method
// Implemented Using MPI
// -----------------------------------------------------------------------------
// Normal compilation:
//    mpicc -lm main.c
// Optimized compilation:
//    mpicc -lm -O3 main.c
// To run:
//    lamboot
//    mpirun -np 8 a.out dataFile1.data 90000 90000.txt

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX_FLOAT 1000
#define MIN_CLUSTERS 9
#define E_DIST(a, b) sqrtf(pow(a.x - b.x, 2) + pow(a.y - b.y, 2))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define BLOCK_START(pid, np, n) ((pid * n) / np)
#define BLOCK_END(pid, np, n) ((((pid + 1) * n) / np) - 1)
// get processor ID from element ID
#define EID_to_PID(eid, np, n) ((((eid + 1) * np) - 1) / n)

typedef struct dFloat {
    float x;
    float y;
} point2D;

typedef struct dInt {
    int i;
    int j;
} location2D;

typedef struct floatInt{
    float distance;
    int location;
} distLoc;

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

    int pid, np, master;
    int block_start, block_end, block_size;
    distLoc minLocSrc, minLocDst;
    location2D min_loc_bcast;
    int i_pid, j_pid;

    // ---------------------------------------------------------- Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    // Set last processor as the master processor
    master = np - 1;

    // ---------------------------------------------------------- Starting Point
    MPI_Barrier(MPI_COMM_WORLD);
    t_start = time(NULL);

    // ------------------------------------------------------- Distribute Points
    MPI_Barrier(MPI_COMM_WORLD);
    t_tic = time(NULL);
    
    // Calculate block low and high value and the size
    ttl_points = atoi(argv[2]);
    block_start = BLOCK_START(pid, np, ttl_points);
    block_end = BLOCK_END(pid, np, ttl_points);
    block_size = (block_end - block_start) + 1;

    // allocate memory for points
    points = (point2D *)malloc(ttl_points * sizeof(point2D));
    if (pid == master) {
        fin_p = fopen(argv[1], "r");
        fout_p = fopen(argv[3], "w");
        for (i = 0; i < ttl_points; i++) {
            fscanf(fin_p, "%f %f", &points[i].x, &points[i].y);
        }
    }
    MPI_Bcast(points, ttl_points * 2, MPI_FLOAT, master, MPI_COMM_WORLD);

    t_diff = difftime(time(NULL), t_tic);
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stdout, "%2d: TIME - Allocate Points                    : %lf seconds\n", pid, t_diff);
    if (pid == master) {
        fprintf(fout_p, "%2d: TIME - Allocate Points                    : %lf seconds\n", pid, t_diff);
    }

    // ----------------------------------------- Allocate And Calculate Distance
    MPI_Barrier(MPI_COMM_WORLD);
    t_tic = time(NULL);

    // allocate memory for distance matrix
    dist_m = (float **)malloc(ttl_points * sizeof(float *));
    for (i = 0; i < ttl_points; i++) {
        dist_m[i] = (float *)malloc(block_size * sizeof(float));
    }
    // calculate distances and store them at distance matrix
    for (i = 0; i < ttl_points; i++) {
        for (j = 0; j < block_size; j++) {
            dist_m[i][j] = E_DIST(points[i], points[j + block_start]);
        }
    }
    for (j = 0; j < block_size; j++) {
        dist_m[j + block_start][j] = MAX_FLOAT;
    }

    t_diff = difftime(time(NULL), t_tic);
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stdout, "%2d: TIME - Allocate And Calculate Distance    : %lf seconds\n", pid, t_diff);
    if (pid == master) {
        fprintf(fout_p, "%2d: TIME - Allocate And Calculate Distance    : %lf seconds\n", pid, t_diff);
    }

    // ------------------------------------------------------- Memory Allocation
    // allocate memory for minimum row
    min_row_dist = (float *)malloc(block_size * sizeof(float));
    min_row_loc = (int *)malloc(block_size * sizeof(int));

    // allocate memory for whole clusters
    clusters = (int *)malloc(ttl_points * sizeof(int));
    for (j = 0; j < ttl_points; j++) {
        clusters[j] = j;
    }

    // ------------------------------------------------ Get Row Minimum Distance
    MPI_Barrier(MPI_COMM_WORLD);
    t_tic = time(NULL);
    
    // get minimum distance for each row
    for (j = 0; j < block_size; j++) {
        min_dist = MAX_FLOAT;
        for (i = 0; i < ttl_points; i++) {
            dist = dist_m[i][j];
            if (dist < min_dist) {
                min_dist = dist;
                min_loc = i;
            }
        }
        min_row_dist[j] = min_dist;
        min_row_loc[j] = min_loc;
    }
    t_diff = difftime(time(NULL), t_tic);
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stdout, "%2d: TIME - Get Minimum Distance Per Row       : %lf seconds\n", pid, t_diff);
    if (pid == master) {
        fprintf(fout_p, "%2d: TIME - Get Minimum Distance Per Row       : %lf seconds\n", pid, t_diff);
    }

    // ---------------------------------------------------- Clustering Main Loop
    MPI_Barrier(MPI_COMM_WORLD);
    t_tic = time(NULL);
    for (loop = 0; loop < ttl_points - MIN_CLUSTERS; loop++) {
        if (pid == master) {
            if (!(loop % 1000)) fprintf(stdout, "loop: %d\n", loop);
        }

        // ----------------------------------------- Get Global Minimum Distance
        // get global minimum distance from minimum row distance
        min_dist = MAX_FLOAT;
        for (j = 0; j < block_size; j++) {
            dist = min_row_dist[j];
            if (dist < min_dist) {
                min_dist = dist;
                min_loc = j;
            }
        }
        // get global minimum using MPI reduce function
        minLocSrc.distance = min_dist;
        minLocSrc.location = pid;
        MPI_Reduce(&minLocSrc, &minLocDst, 1, MPI_FLOAT_INT, MPI_MINLOC, \
                   master, MPI_COMM_WORLD);
        // master broadcasts the processor that has minimum distance
        MPI_Bcast(&minLocDst.location, 1, MPI_INT, master, MPI_COMM_WORLD);
        // processor with minimum distance broadcasts the location
        min_loc_bcast.i = min_row_loc[min_loc];
        min_loc_bcast.j = block_start + min_loc;
        MPI_Bcast(&min_loc_bcast, 2, MPI_INT, minLocDst.location, \
                  MPI_COMM_WORLD);
        // each processor now has the minimum location to update distance matrix
        min_loc_i = min_loc_bcast.i;
        min_loc_j = min_loc_bcast.j;
        // get processors reponsible for min_loc_i and min_loc_j
        i_pid = EID_to_PID(min_loc_i, np, ttl_points);
        j_pid = EID_to_PID(min_loc_j, np, ttl_points);

        // ------------------------------------------------ Duplication Checking
        if (clusters[min_loc_i] != min_loc_i) {
            if (pid == i_pid) {
                min_row_dist[min_loc_i - block_start] = MAX_FLOAT;
            }
            min_loc_i = clusters[min_loc_i];
            i_pid = EID_to_PID(min_loc_i, np, ttl_points);
        }

        // ---------------------------------------------- Update Distance Matrix
        // combine minimum values of new cluster
        for (j = 0; j < block_size; j++) {
            if ((dist_m[min_loc_i][j] < MAX_FLOAT) && \
                (dist_m[min_loc_j][j] < MAX_FLOAT)) {
                dist_m[min_loc_i][j] = MIN(dist_m[min_loc_i][j], \
                                           dist_m[min_loc_j][j]);
            } else {
                dist_m[min_loc_i][j] = MAX_FLOAT;
            }
        }

        // remove one of column combined cluster
        for (j = 0; j < block_size; j++) {
            dist_m[min_loc_j][j] = MAX_FLOAT;
        }

        if (pid == j_pid) {
            j = min_loc_j - block_start;
            // remove current minimum distance
            dist_m[min_loc_i][j] = MAX_FLOAT;
            // remove the minimum row distance as well
            min_row_dist[j] = MAX_FLOAT;
        }

        // ----------------------------------------------------- Update Clusters
        ci = clusters[min_loc_i];
        cj = clusters[min_loc_j];
        for (i = 0; i < ttl_points; i++) {
            if (clusters[i] == cj) {
                clusters[i] = ci;
            }
        }

        // ------------------------------- Update Row Minimum Value and Location
        // update row minimum value and location of new combined cluster
        min_dist = MAX_FLOAT;
        for (j = 0; j < block_size; j++) {
            dist = dist_m[min_loc_i][j];
            if (dist < min_dist) {
                min_dist = dist;
                min_loc = j;
            }
        }
        minLocSrc.distance = min_dist;
        minLocSrc.location = block_start + min_loc;
        MPI_Reduce(&minLocSrc, &minLocDst, 1, MPI_FLOAT_INT, MPI_MINLOC, \
                   i_pid, MPI_COMM_WORLD);
        if (pid == i_pid) {
            j = min_loc_i - block_start;
            min_row_dist[j] = minLocDst.distance;
            min_row_loc[j] = minLocDst.location;
        }

        // update location
        for (j = 0; j < block_size; j++) {
            if (min_row_loc[j] == min_loc_j) {
                min_row_loc[j] = min_loc_i;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    t_diff = difftime(time(NULL), t_tic);
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stdout, "%2d: TIME - Clustering                         : %lf seconds\n", pid, t_diff);
    if (pid == master) {
        fprintf(fout_p, "%2d: TIME - Clustering                         : %lf seconds\n", pid, t_diff);
    }

    t_diff = difftime(time(NULL), t_start);
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stdout, "%2d: TIME - Total                              : %lf seconds\n", pid, t_diff);
    if (pid == master) {
        fprintf(fout_p, "%2d: TIME - Total                              : %lf seconds\n", pid, t_diff);
    }

    // store final cluster at output file
    if (pid == master) {
        fprintf(fout_p, "\nClusters :\n");
        for (i = 0; i < ttl_points; i++) {
            fprintf(fout_p, "%d ", clusters[i]);
        }
        fprintf(fout_p, "\n");
    }

    if (pid == master) {
        fprintf(stdout, "%2d: Done!\n", pid);
        fclose(fin_p);
        fclose(fout_p);
    }
    free(points);
    for (i = 0; i < ttl_points; i++) {
        free(dist_m[i]);
    }
    free(dist_m);
    free(min_row_dist);
    free(min_row_loc);
    free(clusters);
    fflush(stdout);
    MPI_Finalize();
    return 0;
}
