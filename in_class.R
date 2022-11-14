library(pika)

x <- readRDS(here::here("lvMultiSim.Rds"))


# automatically runs likelihood maximization etc
compNB <- sad(x, model = "tnegb")

plot(compNB, ptype = "rad", log = 'x')

library(roleR)


p <- untbParams(1000, 100000, 200, 0.01, 0.1, 'oceanic_island',
                50000, 50000)
neutMod <- roleModel(p)
# neutMod <- roleModel(p)
# Error in qfun(seq(1, 1/S, length = S) - 1/(2 * S)) :
#     could not find function "qfun"

neutMod <- iterModel(neutMod)
neutModFin <- getFinalState(neutMod)
y <- getSumStats(neutModFin, list(abund = rawAbundance))
y <- y$abund$abund
y <- y[y > 0]
