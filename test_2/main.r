library('MASS')
library('volesti')
library('geigen')

randomize_direction <- function(number, dimension) {
  samples = mvrnorm(n = number, mu = rep(0, dimension), Sigma = diag(dimension))
  direction = apply(samples, 1, function(x){
    x/sqrt(sum(x^2))
  })
  return(direction)
}

constraint = function(vector, spec) {
  cons = spec@matrices[[1]]
  for (x in 1:length(vector)) {
    cons = cons + vector[x]*spec@matrices[[x+1]]    
  }
  return(cons)
}

sampling_from_boundary = function(number, S) {
  dimen = length(S@matrices)-1
  direction = randomize_direction(number = (number+1), dimension = dimen)
  boundary_points = c()
  
  for (i in 1:number) {
    random_direction = direction[,i]
    B = constraint(random_direction, spec = S) - S@matrices[[1]]
    A = constraint(rep(0, dimen), spec = S)
    generalized_eigens = geigen::geigen(A, -B, symmetric = FALSE)$values
    positive_eigens = generalized_eigens[generalized_eigens>0]
    negative_eigens = generalized_eigens[generalized_eigens<0]
    
    if (length(positive_eigens) != 0) {
      max_lambda = min(positive_eigens)
    } else {
      max_lambda = Inf
    }
    if (length(negative_eigens) != 0) {
      min_lambda = max(negative_eigens)
    } else {
      min_lambda = -Inf
    }
    boundary_points = c(boundary_points, max_lambda*random_direction,
                        min_lambda*random_direction)
  }
  df = t(matrix(boundary_points, nrow = dimen, ncol = 2*number, byrow = TRUE))
  sampled_points = df[1:number,]
  return(sampled_points)
}


#Example usage:
A0 = matrix(c(-1,0,0,0,-2,1,0,1,-2), nrow=3, ncol=3, byrow = TRUE)
A1 = matrix(c(-1,0,0,0,0,1,0,1,0), nrow=3, ncol=3, byrow = TRUE)
A2 = matrix(c(0,0,-1,0,0,0,-1,0,0), nrow=3, ncol=3, byrow = TRUE)
A3 = matrix(c(-1,-2,0,-2,0,0,0,0,1), nrow=3, ncol=3, byrow = TRUE)
lmi = list(A0, A1, A2, A3)
S = Spectrahedron(matrices = lmi)
samples = sampling_from_boundary(number = 100, S = S)