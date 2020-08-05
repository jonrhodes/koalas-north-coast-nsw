# CODE FOR FITTING SPATIAL KOALA DISTRIBUTION MODELS

# required packages
if(!require(ggregplot)) devtools::install_github("gfalbery/ggregplot")
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggregplot)
library(INLA)

# load data, standardise variables, and convert variables to factors where necessary
Data <- read.csv("data/Koala_Data_Final.csv") %>% as_tibble() %>% mutate_at(vars(Preference), list(as.factor))

# load landscape variables and standardise variables
LSData <- read.csv("data/LandscapeBuffers_Final.csv") %>% as_tibble()

# summarise the amount of each preference type
TreeSum <- Data %>% group_by(ID, Preference) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% select(-n) %>% spread(Preference, prop) %>% mutate_at(vars(-group_cols()), list(~ifelse(is.na(.), 0, .))) %>% rename(PrefPerc1 = `1`, PrefPerc2 = `2`, PrefPerc3 = `3`, PrefPerc4 = `4`, PrefPerc5 = `5`, PrefPerc6 = `6`)

# join tables to create final summarised data
JoinData <- Data %>% left_join(LSData, by = "ID") %>% left_join(TreeSum, by = "ID")

# group preference classes
JoinData <- JoinData %>% mutate(PrefGrp = ifelse(as.character(Preference) == "1", "A", "B")) %>% mutate(PrefGrp = as.factor(PrefGrp)) # A = preference 1, B = preferences 2,3,4,5,6
JoinData <- JoinData %>% mutate(PrefPercA = PrefPerc1, PrefPercB = PrefPerc2 + PrefPerc3 + PrefPerc4 + PrefPerc5 + PrefPerc6) # PrefPercA = preference 1, PrefPercB = preferences 2,3,4,5,6

#function for scaling data
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

# standardise variables
JoinData <- JoinData %>% mutate_at(vars(c(DBH, Nitrate, Ammonium, Phosphorous, Elevation, HSuit_1km, Suit_1km, Marg_1km, NotSuit_1km, Cleared_1km, HSuit_2_5km, Suit_2_5km, Marg_2_5km, NotSuit_2_5km, Cleared_2_5km, HSuit_5km, Suit_5km, Marg_5km, NotSuit_5km, Cleared_5km, PrefPerc1, PrefPerc2, PrefPerc3, PrefPerc4, PrefPerc5, PrefPerc6, PrefPercA, PrefPercB)), list(scale2))

# extract coordinates and create mesh
Coords <- JoinData %>% select(X_Coord, Y_Coord)
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
Y <- JoinData %>% select(Pellets) # response
Intercept <- rep(1, nrow(Data)) # intercept
X0 <- model.matrix(~ -1 + PrefGrp + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + PrefPercA, data = JoinData %>% select(PrefGrp, DBH, Nitrate, Ammonium, Phosphorous, Elevation, HSuit_1km, PrefPercA)) # fixed effects model matrix
X <- as.data.frame(X0[,-which(colnames(X0)%in%c("PrefGrpB"))])
ID <- JoinData$ID
Stk.e <- inla.stack(tag = 'est', ## tag id
            data = list(y = Y), ## response
            A = list(1, 1, 1, A), ## two projection matrices
            effects = list(## three elements:
            Intercept = Intercept,
            X = X, ## fixed effect covariates
            ID = ID, ## random effect for site
            s = s)) ## random field effect

# fit models all variables (1km only) to test for importance of random and spatial effects
# no random or spatial effects
Formula1 <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + PrefPercA
# random effects only
Formula2 <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + PrefPercA +
            f(ID, model = 'iid')
#spatial effects only
Formula3 <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + PrefPercA +
            f(s, model = Spde)
# random and spatial effects
Formula4 <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + PrefPercA +
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
Efxplot(list(Fit1, Fit2, Fit3, Fit4), ModelNames = c("Base model", "Random-effect model", "Spatial model", "Random-effect + spatial model"))

# Test Hypotheses - main effects tree scale

# Null (only elevation)
Formula1 <- Pellets ~ -1 + Intercept + Elevation + f(s, model = Spde)

# Only tree species
Formula2 <- Pellets ~ -1 + Intercept + Elevation + PrefGrpA + f(s, model = Spde)

# Only DBH
Formula3 <- Pellets ~ -1 + Intercept + Elevation + DBH + f(s, model = Spde)

# tree species + DBH
Formula4 <- Pellets ~ -1 + Intercept + Elevation + PrefGrpA + DBH + f(s, model = Spde)

Fit1 <- inla(Formula1, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
Fit2 <- inla(Formula2, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
Fit3 <- inla(Formula3, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
Fit4 <- inla(Formula4, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

#compare models by DIC and plot
INLADICFig(list(Fit1,Fit2,Fit3,Fit4)) # MAIN EFFECTS OF BOTH DBH AND TREE PREF CATEGORY SUPPORTED

# now fit best tree scale model with interaction

# tree species + DBH + tree species:DBH
Formula5 <- Pellets ~ -1 + Intercept + Elevation + PrefGrpA + DBH + PrefGrpA:DBH + f(s, model = Spde)

Fit5 <- inla(Formula5, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

#compare models by DIC and plot
INLADICFig(list(Fit1,Fit2,Fit3,Fit4,Fit5)) # INTERACTION IS IMPORTANT

# tree species + DBH + tree species:DBH + soil
Formula6 <- Pellets ~ -1 + Intercept + Elevation + PrefGrpA + DBH + PrefGrpA:DBH + Nitrate + Ammonium + Phosphorous + f(s, model = Spde)

Fit6 <- inla(Formula6, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

#compare models by DIC and plot
INLADICFig(list(Fit1,Fit2,Fit3,Fit4,Fit5,Fit6)) # DIC NOT LOWER, BUT COMPARABLE

# tree species + DBH + tree species:DBH + soil + tree species:soil
Formula7 <- Pellets ~ -1 + Intercept + Elevation + PrefGrpA + DBH + PrefGrpA:DBH + Nitrate + Ammonium + Phosphorous + PrefGrpA:Nitrate + PrefGrpA:Ammonium + PrefGrpA:Phosphorous + f(s, model = Spde)

Fit7 <- inla(Formula7, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

#compare models by DIC and plot
INLADICFig(list(Fit1,Fit2,Fit3,Fit4,Fit5,Fit6,Fit7)) # DIC NOT LOWER, BUT COMPARABLE

# tree species + DBH + tree species:DBH + proportion preference 1
Formula8 <- Pellets ~ -1 + Intercept + Elevation + PrefGrpA + DBH + PrefGrpA:DBH + PrefPercA + f(s, model = Spde)

Fit8 <- inla(Formula8, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

#compare models by DIC and plot
INLADICFig(list(Fit1,Fit2,Fit3,Fit4,Fit5,Fit6,Fit7,Fit8)) # DIC SLIGHTLY LOWER, BUT COMPARABLE

# tree species + DBH + tree species:DBH + proportion preference 1 + tree species:proportion preference 1
Formula9 <- Pellets ~ -1 + Intercept + Elevation + PrefGrpA + DBH + PrefGrpA:DBH + PrefPercA + PrefGrpA:PrefPercA + f(s, model = Spde)

Fit9 <- inla(Formula9, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

#compare models by DIC and plot
INLADICFig(list(Fit1,Fit2,Fit3,Fit4,Fit5,Fit6,Fit7,Fit8,Fit9)) # STRONG EVIDENCE FOR AN INTERACTION EFFECT

# tree species + DBH + tree species:DBH + proportion preference 1 + tree species:proportion preference 1 + habitat surrounding
Formula10 <- Pellets ~ -1 + Intercept + Elevation + PrefGrpA + DBH + PrefGrpA:DBH + PrefPercA + PrefGrpA:PrefPercA + HSuit_1km + f(s, model = Spde)

Fit10 <- inla(Formula10, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

#compare models by DIC and plot
INLADICFig(list(Fit1,Fit2,Fit3,Fit4,Fit5,Fit6,Fit7,Fit8,Fit9,Fit10)) # DIC NOT LOWER, BUT COMPARABLE

# tree species + DBH + tree species:DBH + proportion preference 1 + tree species:proportion preference 1 + habitat surrounding + tree species: habitat surrounding
Formula11 <- Pellets ~ -1 + Intercept + Elevation + PrefGrpA + DBH + PrefGrpA:DBH + PrefPercA + PrefGrpA:PrefPercA + HSuit_1km + PrefGrpA:HSuit_1km + f(s, model = Spde)

Fit11 <- inla(Formula11, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

#compare models by DIC and plot
INLADICFig(list(Fit1,Fit2,Fit3,Fit4,Fit5,Fit6,Fit7,Fit8,Fit9,Fit10,Fit11)) # DIC NOT LOWER - INTERACTION IMPORTANT

Models <- list(Fit1,Fit2,Fit3,Fit4,Fit5,Fit6,Fit7,Fit8,Fit9,Fit10,Fit11)
saveRDS(Models, "results/models.rds")

#plot coefficients
Efxplot(list(Models[[11]]), ModelNames = c("Best model"))

#gproj <- inla.mesh.projector(Mesh,  dims = c(300, 300))
#g.mean <- inla.mesh.project(gproj, Models[[7]]$summary.random$mean)
