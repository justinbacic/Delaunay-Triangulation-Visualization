import Triangulation as Tri
import time
import random
from numpy import random
import numpy as np
import matplotlib.pyplot as plt

#Generates points in a spiral pattern, this serves as one of the worst case scenarios for the triangulation
def caseSpiral(num_points,num_runs, randomized=False):
    theta = np.linspace(0, 10 * np.pi, num_points)  # Angle increases to create a spiral
    radius = np.linspace(0, 5, num_points)  # Radius increases outward

    # Convert polar coordinates to Cartesian
    x_values = radius * np.cos(theta)
    y_values = radius * np.sin(theta)

    # Store points in a list as (x, y) tuples in spiral order
    points = [(x_values[i], y_values[i]) for i in range(num_points)]
    average_rumtime = 0
    average_flips = 0
    average_rebuckets = 0
    for i in range(num_runs):
        if(randomized):
            random.shuffle(points)
        triangulation = Tri.Triangulation(points)
        start_time = time.time()
        triangulation.incremental_delaunay()
        average_rumtime += time.time()-start_time
        print("Runtime = ",str(time.time()-start_time))
        #print("Inserted ",triangulation.numInserted,"points")
        average_flips += triangulation.edgeFlips
        average_rebuckets += triangulation.rebuckets
    average_rumtime = average_rumtime / number_of_runs
    average_flips = average_flips / number_of_runs
    average_rebuckets = average_rebuckets / number_of_runs
    print("Average Runtime is ",average_rumtime, ", Average edge flips is ",average_flips,", Average Rebuckets is ", average_rebuckets,sep="")

#Generates completely random points on each run, to serve as a benchmark for the other cases 
def caseBenchmark(num_points,num_runs, randomized = False):
    average_rumtime = 0
    average_flips = 0
    average_rebuckets = 0
    for i in range(num_runs):
        points = []
        for _ in range(num_points):
                point = (round(random.uniform(0, 100),2), (round(random.uniform(0, 100),2)))
                if point not in points:
                    points.append(point)
        triangulation = Tri.Triangulation(points)
        start_time = time.time()
        triangulation.incremental_delaunay()
        average_rumtime += time.time()-start_time
        print("Runtime = ",str(time.time()-start_time))
        #print("Inserted ",triangulation.numInserted,"points")
        average_flips += triangulation.edgeFlips
        average_rebuckets += triangulation.rebuckets
    average_rumtime = average_rumtime / number_of_runs
    average_flips = average_flips / number_of_runs
    average_rebuckets = average_rebuckets / number_of_runs
    print("Average Runtime is ",average_rumtime, ", Average edge flips is ",average_flips,", Average Rebuckets is ", average_rebuckets,sep="")

#Generates points on a line
def caseLinePoints( num_points, num_runs,randomized=False):
    m = 2    # Slope of the line
    b = 1    # Y-intercept
    points = []
    x_start=0
    x_end = num_points * 10 #Prevent collisions of points
    """Generates a specified number of points along a line y = mx + b"""
    x_values = np.linspace(x_start, x_end, num_points)  # Evenly spaced x values
    y_values = m * x_values + b  # Compute corresponding y values
    for i in range(num_points):
        points.append((x_values[i],y_values[i]))
    average_rumtime = 0
    average_flips = 0
    average_rebuckets = 0
    for i in range(num_runs):
        if(randomized):
            random.shuffle(points)
        triangulation = Tri.Triangulation(points)
        start_time = time.time()
        triangulation.incremental_delaunay()
        average_rumtime += time.time()-start_time
        print("Runtime = ",str(time.time()-start_time))
        average_flips += triangulation.edgeFlips
        average_rebuckets += triangulation.rebuckets
    average_rumtime = average_rumtime / number_of_runs
    average_flips = average_flips / number_of_runs
    average_rebuckets = average_rebuckets / number_of_runs
    print("Average Runtime is ",average_rumtime, ", Average edge flips is ",average_flips,", Average Rebuckets is ", average_rebuckets,sep="")

#Generates points on the outside of a circle to triangulate
def caseCircle(num_points, num_runs, randomized=False):
    theta = np.linspace(0, 2 * np.pi, num_points, endpoint=False)  # Evenly spaced angles
    radius = num_points /2  # Fixed radius for a circle

    # Convert polar coordinates to Cartesian
    x_values = radius * np.cos(theta)
    y_values = radius * np.sin(theta)

    # Store points in a list as (x, y) tuples
    points = [(x_values[i], y_values[i]) for i in range(num_points)]
    average_runtime = 0
    average_flips = 0
    average_rebuckets = 0

    for i in range(num_runs):
        if randomized:
            random.shuffle(points)
        triangulation = Tri.Triangulation(points)
        start_time = time.time()
        triangulation.incremental_delaunay()
        average_runtime += time.time() - start_time
        print("Runtime = ", str(time.time() - start_time))
        average_flips += triangulation.edgeFlips
        average_rebuckets += triangulation.rebuckets
    average_runtime /= num_runs
    average_flips /= num_runs
    average_rebuckets /= num_runs

    print("Average Runtime is ", average_runtime, ", Average edge flips is ", average_flips, ", Average Rebuckets is ", average_rebuckets, sep="")

num_points = 2500
number_of_runs = 100
#Instances of different test cases
caseCircle(num_points,number_of_runs, randomized=True)
#caseLinePoints(num_points, number_of_runs)
#caseSpiral(num_points,number_of_runs,True)
#caseBenchmark(num_points,number_of_runs)



