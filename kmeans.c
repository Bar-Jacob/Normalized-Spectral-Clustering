#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

typedef struct Cluster
{
    int num_of_points;
    double* sum_of_points;
    double* centroid;
}Cluster;

void kmeans(double** clusters_spk, int k, 
                        int dimension, int num_of_points);
double Euclidian_Distance(double* vector1, double* vector2, int dimension);
void finding_cluster(double* vector, Cluster* clusters, int k, int dimension);
int update_mean(Cluster* clusters, int same_average, int k, int dimension);
void update_sum_of_elements_in_cluster(double* vector, int loc, Cluster* clusters, int dimension);
void free_clusters(Cluster* clusters, int k, int num_of_points);

void kmeans(double** clusters_spk, int k, 
                        int dimension, int num_of_points)
{
    Cluster* clusters;
    int same_average = 0;
    int cnt = 0;
    int i = 0;
    int j = 0;
    /*
    Initializing k clusters
    */
    clusters = (Cluster*)calloc(k,sizeof(struct Cluster));
    for(i = 0; i < k; i++){
        clusters[i].centroid = (double*)calloc(dimension, sizeof(double));
        assert(clusters[i].centroid != NULL && "An Error Has Occured");

        memcpy(clusters[i].centroid, clusters_spk[i], sizeof(double)*dimension); /*will be equal to the i'th vector
        */
        clusters[i].num_of_points = 0;
        clusters[i].sum_of_points = (double*)calloc(dimension, sizeof(double));
        assert(clusters[i].sum_of_points != NULL && "An Error Has Occured");
    }

    cnt = 0;
    while ((cnt < 300) && (!same_average))
    {
        same_average = 1;
        for(i = 0; i < num_of_points; i++){
            finding_cluster(clusters_spk[i], clusters, k, dimension);
        }
        
        same_average = update_mean(clusters, same_average, k, dimension);
        if(same_average == 1){
            break;
        }

        for(i = 0; i < k; i++){
            clusters[i].num_of_points = 0;
            for(j = 0; j < dimension; j++){
                clusters[i].sum_of_points[j] = 0;
            }
        }
        cnt++;
    }

    for(i = 0; i < k; i++){
        for(j = 0; j < dimension; j++){
            if(j == dimension-1){
                printf("%.4f", clusters[i].centroid[j]);
            }else{
                printf("%.4f,", clusters[i].centroid[j]);
            }
        }
        printf("\n");
    }
    free_clusters(clusters, k, num_of_points);
}

void finding_cluster(double* vector, Cluster* clusters, int k, int dimension){
    double min_distance = -1.0;
    int num_of_cluster = -1;
    double distance;
    int i = 0;

    for(i = 0; i < k; i++){
        distance = Euclidian_Distance(vector, clusters[i].centroid, dimension);
        if((distance < min_distance) || (min_distance < 0)){
            min_distance = distance;
            num_of_cluster = i;
        }
    }
    clusters[num_of_cluster].num_of_points++;
    update_sum_of_elements_in_cluster(vector, num_of_cluster, clusters, dimension);
}

void update_sum_of_elements_in_cluster(double* vector, int loc, Cluster* clusters, int dimension){
    int i = 0;
    for(i = 0; i < dimension; i++){
        clusters[loc].sum_of_points[i] += vector[i];
    }
}


int update_mean(Cluster* clusters, int same_average, int k, int dimension){
    int i = 0;
    int j = 0;
    for(i = 0; i < k; i++){
        for(j = 0; j < dimension; j++){
            if((clusters[i].sum_of_points[j]/clusters[i].num_of_points)!=
                clusters[i].centroid[j]){
                    same_average = 0;
                    clusters[i].centroid[j] = clusters[i].sum_of_points[j]/clusters[i].num_of_points;
                }
        }
    }
    return same_average;
}

double Euclidian_Distance(double* vector1, double* centroid, int dimension){
    double sum = 0.0;
    int xi = 0;
    for(xi = 0; xi < dimension; xi++){
        sum += (vector1[xi]-centroid[xi])*(vector1[xi]-centroid[xi]);
    } 
    return sum;
}

void free_clusters(Cluster* clusters, int k, int num_of_points){
    int j = 0;
    
    for(j = 0; j < k; j++){
        free(clusters[j].centroid);
        free(clusters[j].sum_of_points);
    }
    free(clusters);
}