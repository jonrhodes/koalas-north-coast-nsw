# CODE FOR FITTING SPATIAL KOALA DISTRIBUTION MODELS

# required packages
if(!require(ggregplot)) devtools::install_github("gfalbery/ggregplot")
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggregplot)

# load data, standardise variables, and convert variables to factors where necessary
Data <- read.csv("data/Koala_Data_Final.csv") %>% as_tibble() %>% mutate_at(vars(Preference), list(as.factor)) %>%
        mutate_at(vars(c(DBH, Nitrate, Ammonium, Phosphorous, Lscape_PERC, Elevation)), list(scale))

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

# create named index of vectors
s <- inla.spde.make.index('s', n.spde = Spde$n.spde)

# create a data stack
Y <- Data %>% select(Pellets) # response
Intercept <- rep(1, nrow(Data)) # intercept
X0 <- model.matrix(~ -1 + Preference + DBH + Nitrate + Ammonium + Phosphorous + Lscape_PERC + Elevation,
                data = Data %>% select(Preference, DBH, Nitrate, Ammonium, Phosphorous, Lscape_PERC, Elevation)) # fixed effects model matrix
X <- as.data.frame(X0[,-which(colnames(X0)%in%c("Preference6"))])
ID <- Data$ID
Stk.e <- inla.stack(tag = 'est', ## tag id
            data = list(y = Y), ## response
            A = list(1, 1, 1, A), ## two projection matrices
            effects = list(## three elements:
            Intercept = Intercept,
            X = X, ## fixed effect covariates
            ID = ID, ## random effect for site
            s = s)) ## random field effect

# fit model
#no random or spatial effects
Formula1 <- Pellets ~ -1 + Intercept + Preference1 + Preference2 + Preference3 + Preference4 + Preference5 + DBH + Nitrate + Ammonium + Phosphorous + Lscape_PERC + Elevation
# random effects only
Formula2 <- Pellets ~ -1 + Intercept + Preference1 + Preference2 + Preference3 + Preference4 + Preference5 + DBH + Nitrate + Ammonium + Phosphorous + Lscape_PERC + Elevation +
            f(ID, model = 'iid')
#spatial effects only
Formula3 <- Pellets ~ -1 + Intercept + Preference1 + Preference2 + Preference3 + Preference4 + Preference5 + DBH + Nitrate + Ammonium + Phosphorous + Lscape_PERC + Elevation +
            f(s, model = Spde)
# random and spatial effects
Formula4 <- Pellets ~ -1 + Intercept + Preference1 + Preference2 + Preference3 + Preference4 + Preference5 + DBH + Nitrate + Ammonium + Phosphorous + Lscape_PERC + Elevation +
            f(ID, model = 'iid') + f(s, model = Spde)

Fit1<- inla(Formula1, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
Fit2<- inla(Formula2, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
Fit3<- inla(Formula3, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
Fit4<- inla(Formula4, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

#compare models by DIC and plot
INLADICFig(list(Fit1,Fit2,Fit3,Fit4), ModelNames = c("Base", "Random-effect", "Spatial", "Random-effect + Spatial"))

#plot coefficients
Efxplot(list(Fit1, Fit3), ModelNames = c("Base model", "Spatial model"))
