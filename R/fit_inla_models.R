EXAMPLE 1

### Define coordinates
 coords<-as.matrix(data[,1:2])

 ### Boundary
 bnd<-inla.nonconvex.hull(coords, convex=-0.1)

 ### Mesh
 cmv.mesh<-inla.mesh.2d(boundary=bnd, offset=c(.06,.1), max.edge=c(.05,.2))


 ### Create SPDE object
 cmv.spde<-inla.spde2.matern(mesh=cmv.mesh, alpha=2)

 ### List of indexes for the spatial effect
 s.index<-inla.spde.make.index("spatial", n.spde = cmv.spde$n.spde)

 ### Projector matrix
 A.est<-inla.spde.make.A(mesh=cmv.mesh, loc=coords, index=rep(1:nrow(coords)))
 dim(A.est)
[1] 3504  800

 cmv.stack.est<-inla.stack(
+   data=list(y=data$result),
+   tag="est",
+   A=list(A.est),
+   effects=list(c(s.index, list(intercept=1))))

 formula<-  y~ -1 + intercept + f(spatial, model=cmv.spde)

output1<- inla(formula, data = list(y=data$result), family="binomial", Ntrials=1)
Error in eval(expr, envir, enclos) : object 'intercept' not found

# CODE

# required packages
library(INLA)
library(tidyverse)

# load data
Data <- read.csv("data/Koala_Data_Final.csv")
#convert to tibble
Data <- as_tibble(Data)

# extract coordinates and create mesh
Coords <- Data %>% select(X_Coord, Y_Coord)
Mesh <- inla.mesh.2d(Coords, max.edge = c(5000, 10000), offset = c(10000, 20000))

#plot mesh
plot(Mesh, asp = 1)
points(Coords, col = 'red')

# create A matrix to map points onto mesh nodes
A <- inla.spde.make.A(mesh = Mesh, loc = as.matrix(Coords))

# create SPDE object
Spde <- inla.spde2.matern(mesh = Mesh, alpha = 1.5)

# create a data stack
Y <- Data %>% select(Pellets) # response
X <- Data %>% select(DBH) # covariates
Stk.e <- inla.stack(tag = 'est', ## tag id
            data = list(y = Y), ## response
            A=list(1, A), ## two projection matrices
            effects=list(## two elements:
            data.frame(b0 = 1, x = X), ## covariate
            s = 1:Spde$n.spde)) ## RF index

# fit model
Formula <- Pellets ~ 0 + b0 + DBH + ## fixed part
      f(s, model = Spde) ## RF term
Fit <- inla(Formula, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
      control.predictor = list(compute = TRUE, A = inla.stack.A(Stk.e)))

round(Fit$summary.fixed, 4)


#prepare data for model fit




#fit INLA model
