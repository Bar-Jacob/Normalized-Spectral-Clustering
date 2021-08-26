#ifndef KMEANS_H_
#define KMEANS_H_
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

typedef struct Cluster
{
    int num_of_points;
    double *sum_of_points;
    double *centroid;
} Cluster;

void kmeans(double **clusters_spk, int k,
            int dimension, int num_of_points);
void shared_kmeans_algorithm(Cluster *clusters, double **clusters_spk, int k,
                             int dimension, int num_of_points);
void kmeans_pp(int k, int dimension_p, int num_of_points_p,
               PyObject *centroids_locations, PyObject *data_points_p);
double Euclidian_Distance(double *vector1, double *vector2, int dimension);
void finding_cluster(double *vector, Cluster *clusters, int k, int dimension);
int update_mean(Cluster *clusters, int same_average, int k, int dimension);
void update_sum_of_elements_in_cluster(double *vector, int loc, Cluster *clusters, int dimension);
void free_clusters(Cluster *clusters, int k, int num_of_points);
void print_clusters(Cluster *clusters, int dimension, int k);

#endif
