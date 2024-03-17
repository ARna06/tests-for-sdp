library("volesti")
library("rgl")
polytope <- gen_cube(100, "H")

samples <- sample_points(polytope, n = 1000, random_walk = list(
    "walk" = "CDHR", "burn-in" = 1000, "walk_length" = 5
))
#here we use coordinate direction hit and run random walk scheme,
#to sample the points from the hypercube.

#we generate a plot of the sampled points using only the coordinates.
plot3d(samples, xlab = "x", ylab = "y", zlab = "z") 
rglwidget()