import spkmeansmodule
import pandas as pd
import numpy as np
import sys

#recieving input from cmd
k = int(sys.argv[1])
goal = sys.argv[2][0]
input_file = sys.argv[3]

#reading the given files
f1 = open(input_file, 'r')
len_f1 = len(f1.readline().split(","))
f1.close

#convert the files into dataframes + give the titles to the colums 
points = pd.read_csv(input_file, names=["c" + str(i) for i in range(len_f1)])
points = pd.DataFrame(points)

points_numpy = points.to_numpy() #convert to numpy array

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
print("before fit")
spkmeansmodule.fit(k, goal, dimension, data_points_size, data_points_p)