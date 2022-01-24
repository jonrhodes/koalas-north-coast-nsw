# CODE FOR FITTING SPATIAL KOALA DISTRIBUTION MODELS

# load and install required packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(RColorBrewer)) install.packages("RColorBrewer")
if(!require(INLA)) install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
if(!require(ggregplot)) devtools::install_github("gfalbery/ggregplot")
if(!require(MCMCglmm)) install.packages("MCMCglmm")
if(!require(jpeg)) install.packages("jpeg")

# load data, standardise variables, and convert variables to factors where necessary
Data <- read.csv("input/Koala_Data_Final.csv") %>% as_tibble() %>% mutate_at(vars(Preference), list(as.factor))

# load landscape variables and standardise variables
LSData <- read.csv("input/LandscapeBuffers_Final.csv") %>% as_tibble()

# summarise the amount of each preference type
TreeSum <- Data %>% group_by(ID, Preference) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% dplyr::select(-n) %>% spread(Preference, prop) %>% mutate_at(vars(-group_cols()), list(~ifelse(is.na(.), 0, .))) %>% rename(PrefPerc1 = `1`, PrefPerc2 = `2`, PrefPerc3 = `3`, PrefPerc4 = `4`, PrefPerc5 = `5`, PrefPerc6 = `6`)

# join tables to create final summarised data
JoinData <- Data %>% left_join(LSData, by = "ID") %>% left_join(TreeSum, by = "ID")

# group preference classes
JoinData <- JoinData %>% mutate(PrefGrp = ifelse(as.character(Preference) == "1", "A", "B")) %>% mutate(PrefGrp = as.factor(PrefGrp)) # A = preference 1, B = preferences 2,3,4,5,6
JoinData <- JoinData %>% mutate(PrefPercA = PrefPerc1, PrefPercB = PrefPerc2 + PrefPerc3 + PrefPerc4 + PrefPerc5 + PrefPerc6) # PrefPercA = preference 1, PrefPercB = preferences 2,3,4,5,6

#function for scaling data
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

# standardise variables
JoinData <- JoinData %>% mutate_at(vars(c(DBH, Nitrate, Ammonium, Phosphorous, Elevation, HSuit_1km, Suit_1km, Marg_1km, NotSuit_1km, Cleared_1km, HSuit_2_5km, Suit_2_5km, Marg_2_5km, NotSuit_2_5km, Cleared_2_5km, HSuit_5km, Suit_5km, Marg_5km, NotSuit_5km, Cleared_5km, PrefPerc1, PrefPerc2, PrefPerc3, PrefPerc4, PrefPerc5, PrefPerc6, PrefPercA, PrefPercB)), list(scale2))

# get correlation matrix for covariates and save
CorMatrix <- cor(JoinData %>% dplyr::select(DBH, Nitrate, Ammonium, Phosphorous, Elevation, HSuit_1km, HSuit_2_5km, HSuit_5km, PrefPercA))
write.csv(CorMatrix, file="output/corr_matrix.csv")
# inspect cross plots
plot(JoinData %>% dplyr::select(DBH, Nitrate, Ammonium, Phosphorous, Elevation, HSuit_1km, HSuit_2_5km, HSuit_5km, PrefPercA))

# extract coordinates
Coords <- JoinData %>% dplyr::select(X_Coord, Y_Coord)

# create mesh
max.edge.length <- c(1500, 3000)
Mesh <- inla.mesh.2d(loc = Coords, max.edge =  max.edge.length, cutoff = 200)

# plot mesh and inspect
plot(Mesh, asp = 1)
points(Coords, col = 'red')

# create A matrix to map points onto mesh nodes
A <- inla.spde.make.A(mesh = Mesh, loc = as.matrix(Coords))

# create SPDE object
Spde <- inla.spde2.matern(mesh = Mesh, alpha = 1.5)

# create named index of vectors
s <- inla.spde.make.index('s', n.spde = Spde$n.spde)

# create a data stack
Y <- JoinData %>% dplyr::select(Pellets) # response
Intercept <- rep(1, nrow(Data)) # intercept
X0 <- model.matrix(~ -1 + PrefGrp + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + HSuit_2_5km + HSuit_5km + PrefPercA, data = JoinData %>%
            dplyr::select(PrefGrp, DBH, Nitrate, Ammonium, Phosphorous, Elevation, HSuit_1km, HSuit_2_5km, HSuit_5km, PrefPercA)) # fixed effects model matrix
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

# fit models with all variable main and interaction effects to test for importance of the random and spatial effects
# and scale of the landscape variables at the 1km, 2.5 km and 5km scale

# no random or spatial effects - 1km
FormulaNone1km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_1km:PrefGrpA
# random effects only - 1km
FormulaRand1km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_1km:PrefGrpA +
            f(ID, model = 'iid')
# spatial effects only - 1km
FormulaSpat1km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_1km:PrefGrpA +
            f(s, model = Spde)
# random and spatial effects - 1km
FormulaRandSpat1km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_1km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_1km:PrefGrpA  +
            f(ID, model = 'iid') + f(s, model = Spde)
# fit models
FitNone1km <- inla(FormulaNone1km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
FitRand1km <- inla(FormulaRand1km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
FitSpat1km <- inla(FormulaSpat1km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
FitRandSpat1km <- inla(FormulaRandSpat1km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

# no random or spatial effects - 2.5km
FormulaNone2_5km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_2_5km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_2_5km:PrefGrpA
# random effects only - 2.5km
FormulaRand2_5km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_2_5km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_2_5km:PrefGrpA +
            f(ID, model = 'iid')
# spatial effects only - 2.5km
FormulaSpat2_5km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_2_5km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_2_5km:PrefGrpA +
            f(s, model = Spde)
# random and spatial effects - 2.5km
FormulaRandSpat2_5km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_2_5km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_2_5km:PrefGrpA +
            f(ID, model = 'iid') + f(s, model = Spde)
# fit models
FitNone2_5km <- inla(FormulaNone2_5km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
FitRand2_5km <- inla(FormulaRand2_5km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
FitSpat2_5km <- inla(FormulaSpat2_5km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
FitRandSpat2_5km <- inla(FormulaRandSpat2_5km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

# no random or spatial effects - 5km
FormulaNone5km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_5km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_5km:PrefGrpA
# random effects only - 5km
FormulaRand5km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_5km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_5km:PrefGrpA +
            f(ID, model = 'iid')
# spatial effects only - 5km
FormulaSpat5km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_5km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_5km:PrefGrpA +
            f(s, model = Spde)
# random and spatial effects - 5km
FormulaRandSpat5km <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_5km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_5km:PrefGrpA +
            f(ID, model = 'iid') + f(s, model = Spde)

# fit models
FitNone5km <- inla(FormulaNone5km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
FitRand5km <- inla(FormulaRand5km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
FitSpat5km <- inla(FormulaSpat5km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
FitRandSpat5km <- inla(FormulaRandSpat5km, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))

# save models
ModelsRandSpat <- list(FitNone1km, FitRand1km, FitSpat1km, FitRandSpat1km, FitNone2_5km, FitRand2_5km, FitSpat2_5km, FitRandSpat2_5km, FitNone5km, FitRand5km, FitSpat5km, FitRandSpat5km)
saveRDS(ModelsRandSpat, file = "output/models_rand_spat.rds")

# compare models by DIC and plot
INLADICFig(ModelsRandSpat)

# choose best model - NOTE NEED TO CHANGE DEPENDING ON WHAT THE BEST MODEL IS
BestModel <- ModelsRandSpat[[12]]

# check the range parameter for the best model
rf <- inla.spde.result(inla = BestModel, name = "s", spde = Spde, do.transf = TRUE)
plot(rf$marginals.range[[1]], type = "l", xlab = "Range (m)", ylab = "Density", xlim = c(0, 20000))
abline(v = max.edge.length[2], col = 'red')

# Test Hypotheses - main effects only
FinalFormulaMain <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_5km + PrefPercA +
            f(ID, model = 'iid') + f(s, model = Spde) # ADJUST FORMULA DEPENFDING ON WHAT THE BEST MODEL IS
FinalFitMain <- inla(FinalFormulaMain, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
saveRDS(FinalFitMain, "output/final_model_main.rds")

# Test Hypotheses - main and interaction effects
# random and spatial effects - 5km
FinalFormulaMainInt <- Pellets ~ -1 + Intercept + PrefGrpA + DBH + Nitrate + Ammonium + Phosphorous + Elevation + HSuit_5km + PrefPercA +
            DBH:PrefGrpA + Nitrate:PrefGrpA + Ammonium:PrefGrpA + Phosphorous:PrefGrpA + PrefPercA:PrefGrpA + HSuit_5km:PrefGrpA +
            f(ID, model = 'iid') + f(s, model = Spde) # ADJUST FORMULA DEPENFDING ON WHAT THE BEST MODEL IS
FinalFitMainInt <- inla(FinalFormulaMainInt, family = "binomial", Ntrials = 1, data = inla.stack.data(Stk.e),
          control.compute = list(dic = TRUE),
          control.predictor = list(A = inla.stack.A(Stk.e)))
saveRDS(FinalFitMainInt, "output/final_model_interactions.rds")

# get coefficients and plot graphs
# main effects model
RNamesMain <- c("Pref 1", "DBH", "Nitrate", "Ammonium", "Phosphorous", "Elevation", "Percent Pref 1", "Hab 5km Buffer")
CoefsMain <- summary(FinalFitMain)$fixed[c(2:7, 9, 8), ] %>% as_tibble() %>%
              mutate(variable = factor(RNamesMain, ordered = TRUE, levels = c("Pref 1", "DBH", "Nitrate", "Ammonium", "Phosphorous", "Elevation", "Percent Pref 1", "Hab 5km Buffer")),
              lower = `0.025quant`, upper = `0.975quant`) %>% dplyr::select(variable, mean, lower, upper)
PlotMain <- CoefsMain %>% ggplot(aes(x = variable, y = mean)) + theme_bw() + geom_col(width = 0.5, fill = "grey") + geom_errorbar(aes(ymin = lower,  ymax = upper), width = .2, position=position_dodge(.9), color = "black") +
              labs(y = "Coefficient Estimates (95% Credible Intervals)", x = "") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("output/fig_main.jpg", PlotMain)
# interaction effects model
RNamesInt <- c("Pref 1 x DBH", "Pref 1 x Nitrate", "Pref 1 x Ammonium", "Pref1 x Phosphorous", "Pref 1 x Percent Pref 1", "Pref 1 x Hab 5km Buffer")
CoefsInt <- summary(FinalFitMainInt)$fixed[10:15, ] %>% as_tibble() %>%
              mutate(variable = factor(RNamesInt, ordered = TRUE, levels = c("Pref 1 x DBH", "Pref 1 x Nitrate", "Pref 1 x Ammonium", "Pref1 x Phosphorous", "Pref 1 x Percent Pref 1", "Pref 1 x Hab 5km Buffer")),
              lower = `0.025quant`, upper = `0.975quant`) %>% dplyr::select(variable, mean, lower, upper)
PlotInt <- CoefsInt %>% ggplot(aes(x = variable, y = mean))  + theme_bw() + geom_col(width = 0.5, fill = "grey") + geom_errorbar(aes(ymin = lower,  ymax = upper), width = .2, position=position_dodge(.9), color = "black") +
              labs(y = "Coefficient Estimates (95% Credible Intervals)", x = "") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("output/fig_int.jpg", PlotInt)
