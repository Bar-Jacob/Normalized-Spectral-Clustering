import spkmeansmodule
import pandas as pd
import numpy as np
import sys

def goal_decode(goal):
    if(goal == 's'):
        return 0
    if(goal == 'w'):
        return 1
    if(goal == 'd'):
        return 2
    if(goal == 'l'):
        return 3
    if(goal == 'j'):
        return 4

#recieving input from cmd
k = int(sys.argv[1])
goal = sys.argv[2][0]
goal = goal_decode(goal)

input_file = sys.argv[3]

#reading the given files
f1 = open(input_file, 'r')
len_f1 = len(f1.readline().split(","))
f1.close

#convert the files into dataframes + give the titles to the colums 
points = pd.read_csv(input_file, names=["c" + str(i) for i in range(len_f1)])
points = pd.DataFrame(points)

points_numpy = np.array(points) #convert to numpy array

if(k != 0 and k >= len(points)):
    print("An Error Has Occured")
    assert()

# setting an argument for C
dimension = len(points_numpy[0])
data_points_size = len(points)

# creating a 1D list of the data points to give as an argument to C
data_points_p = []
for vector in points_numpy:
    for i in range(dimension):
        data_points_p.append(vector[i])

#using kmeanssp module by calling the fit function
#goals legend: spk -> 0, wam -> 1, ddg -> 2, lnorm -> 3, jacobi -> 4
data_points_p = spkmeansmodule.fit(k, goal, dimension, data_points_size, data_points_p)

distances = [-1.0 for i in range(len(data_points_p))]
probabilities = [0.0 for i in range(len(data_points_p))]
indexes = [0 for i in range(len(data_points_p))]

data_points_numpy = np.array(data_points_numpy) #convert to numpy array

centroids = np.array(
    [[0.0 for i in range(len(data_points_numpy[0]))] for i in range(k)])
# creating a list of the locations of the centorinds in the merged_points_numpy list
centroids_locations = [0 for i in range(k)]

#calculating probabilities after calculating distances. 
#the closer you are to one of the centroids, the lower your probability
def probabilities_calc():
    sum = 0.0
    for distance in distances:
        sum += distance
    for i in range(len(distances)):
        probabilities[i] = distances[i]/sum

#calculating the minimum distance of each point from all current clusters 
#by comparing the minimum distance from the last iteration to the distance from the new centroid
def min_distances(Z):
    for index in range(len(data_points_numpy)):
        curr_distance = pow(np.linalg.norm(
            (data_points_numpy[index]-centroids[Z-1])), 2)
        if(curr_distance < distances[index] or distances[index] == -1.0):
            distances[index] = curr_distance

#initializing the first k centroids
def KMeansPP():
    cnt = 1
    np.random.seed(0)

    random_index = (int)(np.random.choice(indexes))
    centroids_locations[0] = random_index
    centroids[0] = np.ndarray.copy(data_points_numpy[random_index])
    Z = 1
    while(Z != k):
        min_distances(Z)
        probabilities_calc()
        random_index = (int)(np.random.choice(
            indexes, p=probabilities))
        centroids_locations[cnt] = random_index
        centroids[Z] = np.ndarray.copy(data_points_numpy[random_index])
        Z += 1
        cnt += 1

KMeansPP()
# setting an argument for C
dimension = len(centroids[0])
data_points_size = len(data_points_numpy)

# creating a 1D list of the data points to give as an argument to C
data_points_p = []
for vector in data_points_numpy:
    for i in range(dimension):
        data_points_p.append(vector[i])

#using kmeanssp module by calling the fit function
np.array(spkmeansmodule.fit_pp(
    k, dimension, data_points_size, centroids_locations, data_points_p))

