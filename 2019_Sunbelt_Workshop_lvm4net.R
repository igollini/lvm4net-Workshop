# install.packages(c("lvm4net", "statnet", "latentnet", "RSiena", "sna", "manet", "sand", "colorspace"))

Y <- RSiena::s501
alcohol <- RSiena::s50a[,1]
smoking <- RSiena::s50s[,1]

library(lvm4net)
D <- 2
modLSM <- lsm(Y, D, nstart = 3)

plot(modLSM, Y, drawCB = TRUE)

gof_s501 <- goflsm(modLSM, Y, doplot = TRUE)

modLSM2 <- lsm(Y, D, sigma = 2.5, xi = -1, psi2 = 1, nstart = 3)
gof_s501_2 <- goflsm(modLSM2, Y, doplot = TRUE)

plot(modLSM, Y, drawCB = TRUE,
  colPl = rev(heat.colors(5))[alcohol], main = "Alcohol Consumption")
legend("topright", legend = c("1-Low", 2:4, "5-High"), 
  col = rev(heat.colors(5)), pch = 20, pt.cex = 2)

plot(modLSM, Y, drawCB = TRUE, 
  colPl = rev(heat.colors(3))[smoking], main = "Smoking")
legend("topright", legend = c("No", "Moderate", "Serious"), 
  col = rev(heat.colors(3)), pch = 20, pt.cex = 2)

data(lazega, package = "sand")
Y <- as.matrix(igraph::get.adjacency(lazega))
A <- igraph::get.data.frame(lazega, what = "vertices")

Y1 <- RSiena::s501
Y2 <- RSiena::s502
Y3 <- RSiena::s503
Y123 <- list(Y1 = Y1, Y2 = Y2, Y3 = Y3)

modLSM1 <- lsm(Y1, D)
modLSM2 <- lsm(Y2, D)
modLSM3 <- lsm(Y3, D)

modLSJM123 <- lsjm(Y123, D)

Z123 <- modLSJM123$EZ

Z1 <- rotXtoY(modLSM1$lsmEZ, Z123)$Xrot
Z2 <- rotXtoY(modLSM2$lsmEZ, Z123)$Xrot
Z3 <- rotXtoY(modLSM3$lsmEZ, Z123)$Xrot

XYlimb <- range(Z123, Z1, Z2, Z3)
namesb <- paste('Network ', 1:3, sep ='')
colPl <- rainbow(50)

par(mfrow = c(1,3))
plotY(Y1, Ndata = 1, EZ = Z1, VZ = modLSM1$lsmVZ, main = namesb[1], xlim = XYlimb, ylim = XYlimb, colPl = colPl)
plotY(Y2, Ndata = 1, EZ = Z2, VZ = modLSM2$lsmVZ, main = namesb[2], xlim = XYlimb, ylim = XYlimb, colPl = colPl)
plotY(Y3, Ndata = 1, EZ = Z3, VZ = modLSM3$lsmVZ, main = namesb[3], xlim = XYlimb, ylim = XYlimb, colPl = colPl)

Zlsm <- list()
Zlsm[[1]] <- Z1
Zlsm[[2]] <- Z2
Zlsm[[3]] <- Z3

bpbLSM <- boxroc(Y123, 
	EZ = Zlsm,
	xiT = c(modLSM1$xiT, modLSM2$xiT, modLSM3$xiT), 
	Lroc = 150, 
	ROC = TRUE, 
	BOXPLOT = TRUE,
	labelsPlot = namesb
	)

plot(modLSJM123, Y123, drawCB = TRUE, plotZtilde = TRUE, colPl = colPl)

plot(modLSJM123, Y123, drawCB = TRUE, colPl = colPl, main = 'Multiple networks')

bpbLSJM <- boxroc(Y123, 
	EZ = modLSJM123$lsmEZ,	
	xiT = modLSJM123$xiT, 
	Lroc = 150, 
	ROC = TRUE, 
	BOXPLOT = TRUE,
	labelsPlot = namesb
	)

par(mfrow = c(1, 2))
plot(modLSJM123, Y123, drawCB = TRUE, colPl = rev(heat.colors(5))[RSiena::s50a[,3]], main = 'Multiple networks - Alcohol')
plot(modLSJM123, Y123, drawCB = TRUE, colPl = rev(heat.colors(3))[RSiena::s50s[,3]], main = 'Multiple networks - Smoking')

load("sampson.Rdata")
Y1 <- (sampson$SAMPES > 0) * 1 # Esteem
Y2 <- (sampson$SAMPIN > 0) * 1 # Positive Influence
Y3 <- (sampson$SAMPPR > 0) * 1 # Praise

data(noordin, package = "manet")

Y <- as.matrix(noordin)

heatmap(
  Y,
  Rowv = NA,
  Colv = NA,
  col = grey(c(0.95, 0.0)),
  scale = "none",
  margins = c(3, 3),
  xlab = "Event",
  ylab = "Terrorist"
  )

G <- 2:4 # is the number of groups
D <- 0:3 # is the dimension of the latent variable

mod.mlta <- mlta(Y, G = G, D = D, wfix = FALSE)
mod.mlta$BIC$`Table of BIC Results`

mod.mlta.wfix <- mlta(Y, G = G, D = 1:3, wfix = TRUE)
mod.mlta.wfix$BIC$`Table of BIC Results`

res <- mod.mlta.wfix[[1]]
plot(c(res$w), xlab = "Event", ylab = "w", pch = 19)
abline(h = 0)

par(mfrow = c(1, 2))

plot(c(res$b[1,]), xlab = "Event", ylab = "b", pch = 19, main = "Group 1")
abline(h = 0)

plot(c(res$b[2,]), xlab = "Event", ylab = "b", pch = 19, main = "Group 2")
abline(h = 0)

plot(res$z[,1], pch = 19, 
  xlab = "Sender node",
  ylab = "Probability to belong to group 1")
abline(h = 0.5, col = "red")

pig0 <- 1 / ( 1 + exp(-res$b))

matplot(t(pig0), type = "l", 
  ylim = c(0, 1), ylab = expression(paste(pi[rg](0))),
  xlab = "Receiver node (r)", xaxt = "n",
  main = "Probability that the median sender node in group g\n has a link with receiver node r")
axis(1, at = 1:ncol(Y))
legend("topright", paste("Group", 1:2, sep = " "), col = 1:2, lty = 1:2)

loglift <- log(lift(res, pdGH = 21))

par(mfrow = c(1, 2))
heatmap(
  loglift[,,1],
  Rowv = NA,
  Colv = NA,
  col = colorspace::diverge_hsv(20),
  breaks = seq(-10, 10, by = 1),
  revC = TRUE,
  scale = "none",
  xlab = "Event",
  ylab = "Event",
  main = "Log-Lift for Group 1"
  )

heatmap(
  loglift[,,2],
  Rowv = NA,
  Colv = NA,
  col = colorspace::diverge_hsv(20),
  breaks = seq(-10, 10, by = 1),
  revC = TRUE,
  scale = "none",
  xlab = "Event",
  ylab = "Event",
  main = "Log-Lift for Group 2"
  )

mod.lca <- lca(Y, G = 2:4)

mod.lta <- lta(Y, D = 1:3)

data(davis, package = "latentnet")
Y <- network::as.matrix.network(davis)

