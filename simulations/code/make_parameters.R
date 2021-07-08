# a grid of simulation parameters
simulations <- c("discrete", "discrete2", "discrete3", 
                 "continuous", "robustness")
ns <- c(250, 500, 1000, 2000)
seed <- 1:1000
parm <- expand.grid(simulation = simulations, 
                    n = ns, 
                    seed = seed)

save(parm, file = here::here("output", "parm.RData"))