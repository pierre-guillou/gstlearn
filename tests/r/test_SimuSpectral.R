knitr::opts_chunk$set(echo = TRUE)
rm(list = ls())
library(gstlearn)

# global parameters
flag.verbose = FALSE
opers = EStatOption_fromKeys(c("NUM", "MINI", "MAXI", "MEAN", "STDV"))

prefix = "Simu"
var_name = paste(prefix, "V1", "S1", sep = ".")

idx_type = 3
type_cov <- c("GAUSSIAN", "EXPONENTIAL", "MATERN")[idx_type]
nu <- c(10, 1/2, 3/4)[idx_type]
ranges = c(5, 10, 15)
angles = c(30., 0, 0)

# initialization of the simulator
nb <- 1000 # number of spectral components
ns <- 20   # number of simulations
seed <- 13112 # Seed for the first simulation

# Define the title (used for plots)
title = paste0(type_cov, " model")
if (type_cov == "MATERN") {title = paste0(title, " (nu = ", nu, ")")}

ndim = 1
# definition of the Space and creation of the grid
err = defineDefaultSpace(ESpaceType_RN(), ndim)
nx   = c(10000)
dx   = c(1.)
grid = DbGrid_create(nx = nx, dx = dx)

# specification of the model using Model
mod  = Model_createFromParam(type = ECov_fromKey(type_cov), 
                             ranges = ranges[1:ndim], param = nu, flagRange=FALSE)
# using the *SimuSpectral* interface
sim = SimuSpectralRN(mod$getCov())
err = sim$simulate(ns = nb, seed = 123, verbose = TRUE)
# 
gamma = sim$getGamma()$toTL()
omega = sim$getOmega()$toTL()
phi   = sim$getPhi()
coor  = grid$getColumnsAsMatrix(names = c("x1"), useSel = FALSE)$toTL()
val   = as.numeric(t(gamma) %*% cos(omega %*% t(coor) + phi))

# # using the interface
err = sim$compute(dbout = grid, verbose = TRUE, namconv = NamingConvention(prefix))
# 
stopifnot (all(abs(val - grid[paste(prefix, "V1", "simu", sep = ".")]) < 1e-6))

# using the simuSpectral interface
err  = simuSpectral(NULL, dbout = grid, cova = mod$getCov(), nbsimu = ns, ns = nb, seed=seed, cov0 = NULL, namconv = NamingConvention(prefix))

# statistics
knitr::kable(dbStatisticsMono(db = grid, names = paste(prefix, "*", sep = "."), opers = opers)$toTL(),
             digits = 3, caption = title)
