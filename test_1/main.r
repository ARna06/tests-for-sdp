library('volesti')
library('rgl')

walk = "RDHR"
b = c(1,1,1,1,1)
A = matrix(c(-1,0,0,0,-1,0,0,0,-1,1,0,0,0,1,0,0,0,1), nrow = 6, ncol = 3, byrow = T)
P = Hpolytope(A=A, b=b)
points = sample_points(P, 100)
plot3d(x = points[1,], y = points[2,], z = points[3,],
       xlab = "X", ylab = "Y", zlab = "Z", 
       type = 'p',
       title3d(main = sprintf("Sampling a random cube with walk %s", walk)))
rglwidget()
