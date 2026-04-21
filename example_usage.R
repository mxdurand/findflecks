source("findzeroes.R")
source("findflecks.R")
source("plotflecks.R")


# Taken on 2021.07.17 in the oat fields of Helsinki university
df <- read.csv("data/example_oat.csv", row.names = 1)

df1 <- df[df$height == "low",]
df2 <- df[df$height == "mid",]
df3 <- df[df$height == "top",]

dfX <- df1
Z <- findZeros(time = dfX$Time, var = dfX$PAR_q)
dfS <- findFlecks(time = dfX$Time, var = dfX$PAR_q, zeroes = Z, minTime = 0, minAmp = 5, minPdiff = 0.05, shadeflecks = F, asmMethod = "max", verbose = F)
plotTSfleckEz(time = dfX$Time, var = dfX$PAR_q, zeroes = Z, fleck_data = dfS)
