#==============================================================================#
#                                                                              #
#                           Modelos N-mixtos                                   #
#                       Jose Jimenez. CSIC-IREC                                #
#                        28/02/2020 10:14:57                                   #
#                          MASTER TOLEDO 2023                                  #
#                                                                              #
#==============================================================================#

# Script extraido del libro:
# Kery, M., & Schaub, M. (2012). Bayesian population analysis using WinBUGS. 
# A hierarchical perspective. Academic Press/Elsevier  
# <http://doi.org/10.1016/B978-0-12-387020-9.00014-6>

## MODELO SIMPLE CON *UNMARKED*
##==============================
### Simulacion de datos
# Elegimos un tamaño de muestra y preparamos la matriz C
set.seed(24)                # Para repetir siempre la misma creacion de datos
M <- 150                    # Numero de sitios
J <- 4                      # Numero de relicas por sitios
C <- matrix(NA, nrow = M, ncol = J) # matriz de datos


# Valores de Parametros
lambda <- 14                # Abundancia esperada
p <- 0.3                    # Probabilidad de deteccion (por individuo)


# PROCESO ECOLOGICO
# Generamos datos de abundancia local
N <- rpois(n = M, lambda = lambda)


# PROCESO DE OBSERVACION
# Simulamos la observacion "J" veces de esos datos
for(j in 1:J){
   C[,j] <- rbinom(n = M, size = N, prob = p)
}


# Veamos los datos creados
# Echemos un vistazo a los datos
# La realidad ....
table(N)                    # Verdadera distribucion de la abundancia
sum(N)                      # Tamaño total de poblacion en los M sitios
sum(N>0)                    # Numero real en los sitios ocupados
mean(N)                     # Verdadera abundancia media (estima de lambda)

#... y lo que observamos
table(apply(C, 1, max))     # distribucion de la abundancia observada (max conteos)
sum(apply(C, 1, max))       # Tamaño de la poblacion observada en los M sitios
sum(apply(C, 1, max)>0)     # Numero de sitios observados como ocupados
mean(apply(C, 1, max))      # Media observada de la abundancia

head(cbind(N=N, count1=C[,1], count2=C[,2], count3=C[,3])) # Veamos los 6 primeros sitios


### Ejecucion del modelo
library(unmarked)                         # Cargamos la libreria
umf <- unmarkedFramePCount(y = C)         # preparamos datos
summary(umf)                              # ... y vemos como son

# Ajustamos los modelos
fm1 <- pcount(~1 ~1, data = umf, K=150, mixture="P")
fm2 <- pcount(~1 ~1, data = umf, K=150, mixture="NB")
fm3 <- pcount(~1 ~1, data = umf, K=150, mixture="ZIP")

f1 <- fitList('Poisson'            =fm1,
              'Negative Binomial'  =fm2,
              'ZIP'                =fm3)

modSel(f1)
##                   nPars     AIC delta AICwt cumltvWt
## Poisson               2 2505.01  0.00  0.57     0.57
## Negative Binomial     3 2506.86  1.84  0.23     0.79
## ZIP                   3 2507.01  2.00  0.21     1.00

# Obtenemos resultados
backTransform(fm1, "state")               # y lo obtenemos a escala natural
backTransform(fm1, "det")


#... y podemos ver el resultado para cada sitio utilizando métodos bayesianos
# empíricos
re<- ranef(fm1)
plot(re, xlim = c(0,20))[sort(sample(1:100, 30))]
# Para obtener los intervalos de confianza
(ppd <- posteriorSamples(re, nsims=10000))
res<-data.frame(Media=apply(ppd@samples,1,mean),
                SD=apply(ppd@samples,1,sd),
                lower=apply(ppd@samples,1,quantile,0.025),
                upper=apply(ppd@samples,1,quantile,0.975))

# Contrastamos la realidad y el resultado
head(cbind(N,res))


# Test de bondad del ajuste
fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2)
    chisq <- sum((observed - expected)^2 / expected)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }

(pb <- parboot(fm1, fitstats, nsim=1000, report=1))

# chisquare test statistic
(c.hat <- pb@t0[2] / mean(pb@t.star[,2])) # c-hat as ratio of observed/expected

par(mfrow=c(2,2))
hist(pb@t.star[,1], xlab='SSE', col="grey90", breaks=15, main="")
abline(v=pb@t0[1],lty=2, lwd=3)
hist(pb@t.star[,2], xlab='Chisq', col="grey90", breaks=15, main="")
abline(v=pb@t0[2],lty=2, lwd=3)
hist(pb@t.star[,3], xlab='freemanTukey', col="grey90", breaks=15, main="")
abline(v=pb@t0[3],lty=2, lwd=3)



## MODELO SIMPLE CON *BUGS*
##==============================
# La estructura del modelo N-mixto es muy sencilla, y similar al modelo de 
# ocupacion, pero en vez de una estructura Binomial-Binomial, es Poisson-Binomial.

y<-C

cat(file = "model.jags",
"
model {

  # Informacion a priori
  lambda ~ dgamma(0.005, 0.005) # Informacion previa para lambda
  p ~ dunif(0, 1)

  # Probabilidad
  # Modelo ecologico para la verdadera abundancia
    for (i in 1:R) {
    N[i] ~ dpois(lambda)
    # Modelo de observacion para los conteos replicados
    for (j in 1:T) {
      y[i,j] ~ dbin(p, N[i])
    } # j
  } # i
}

")


# Preparamos datos
win.data <- list(y = y, R = nrow(y), T = ncol(y))

# Valores de inicio
Nst <- apply(y, 1, max, na.rm=TRUE) + 1	# Importante!!
inits <- function() list(N = Nst)

# Parametros a monitorizar
params <- c("lambda", "p")

# Configuracion MCMC
ni <- 25000
nt <- 10
nb <- 5000
nc <- 3


### Ejecucion del modelo
# Llamamos a JAGS desde R
library(jagsUI)
out <- jags(win.data, inits, params, "model.jags", n.chains = nc, n.thin = nt, 
  n.iter = ni, n.burnin = nb, parallel=FALSE)

# Resumimos resultados
print(out, dig = 2)

# Vemos la convergencia de las cadenas de Markov
xyplot(out$samples)

