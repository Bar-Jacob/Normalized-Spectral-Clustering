#include "kmeans.h"

void kmeans_pp(int k, int dimension, int num_of_points, 
                int* centroids_locations, double** data_points)
{
    Cluster* clusters;
    int i = 0;

    /*
    converting k centroids from python to k centroids in C
    */
    clusters = (Cluster *)calloc(k, sizeof(struct Cluster));
    assert(clusters != NULL && "An Error Has Occured");
    for (i = 0; i < k; i++)
    {
        clusters[i].centroid = (double *)calloc(dimension, sizeof(double));
        assert(clusters[i].centroid != NULL && "An Error Has Occured");

        memcpy(clusters[i].centroid, data_points[centroids_locations[i]], 
        sizeof(double) * dimension); /*will be equal to the i'th vector*/

        clusters[i].num_of_points = 0;
        clusters[i].sum_of_points = (double *)calloc(dimension, sizeof(double));
        assert(clusters[i].sum_of_points != NULL && "An Error Has Occured");
    }
    print_locations(centroids_locations, k);
    shared_kmeans_algorithm(clusters, data_points, k, dimension, num_of_points);
    free(centroids_locations);
}

void kmeans_func(double** clusters_spk, int k, 
                        int dimension, int num_of_points)
{
    Cluster* clusters;
    int i = 0; 

    /*
    Initializing k clusters
    */
    clusters = (Cluster*)calloc(k,sizeof(struct Cluster));
    assert(clusters != NULL && "An Error Has Occured");
    for(i = 0; i < k; i++){
        clusters[i].centroid = (double*)calloc(dimension, sizeof(double));
        assert(clusters[i].centroid != NULL && "An Error Has Occured");

        memcpy(clusters[i].centroid, clusters_spk[i], 
        sizeof(double)*dimension); /*will be equal to the i'th vector*/

        clusters[i].num_of_points = 0;
        clusters[i].sum_of_points = (double*)calloc(dimension, sizeof(double));
        assert(clusters[i].sum_of_points != NULL && "An Error Has Occured");
    }
    shared_kmeans_algorithm(clusters, clusters_spk, k, dimension, num_of_points);
}

void shared_kmeans_algorithm(Cluster* clusters, double** clusters_spk, int k,
                            int dimension, int num_of_points){
    int same_average = 0;
    int cnt = 0;
    int i = 0;
    int j = 0;
    cnt = 0;

    while ((cnt < 300) && (!same_average))
    {
        /*if there's a change in the centroids it will change to 0,
        else it will stay 1*/
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
    print_clusters(clusters, dimension, k);
    free_clusters(clusters, k);
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

void free_clusters(Cluster* clusters, int k){
    int j = 0;
    
    for(j = 0; j < k; j++){
        free(clusters[j].centroid);
        free(clusters[j].sum_of_points);
    }
    free(clusters);
}

void print_clusters(Cluster* clusters, int dimension, int k){
    int i = 0;
    int j = 0;
    for(i = 0; i < k; i++){
        for(j = 0; j < dimension; j++){
            if(j == dimension-1){
                if (clusters[i].centroid[j] < 0 && clusters[i].centroid[j] > -0.00005)
                {
                    clusters[i].centroid[j] = clusters[i].centroid[j] * -1;
                }
                printf("%.4f", clusters[i].centroid[j]);
            }else{
                if (clusters[i].centroid[j] < 0 && clusters[i].centroid[j] > -0.00005)
                {
                    clusters[i].centroid[j] = clusters[i].centroid[j] * -1;
                }
                printf("%.4f,", clusters[i].centroid[j]);
            }
        }
        if(i != k - 1){
            printf("\n");
        }
    }
}

/*prints the locations of the chosen centroids*/
void print_locations(int* array, int num_of_cells)
{
    int i = 0;
    for (i = 0; i < num_of_cells; i++)
    {
        if (i == num_of_cells - 1)
        {
            printf("%d", array[i]);
        }
        else
        {
            printf("%d,", array[i]);
        }
    }
    printf("\n");
}