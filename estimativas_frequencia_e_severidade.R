##### Parte 1: Estudo da freqüência

### Lendo os dados de frequência
frequencias.entrada <- read.table("C:\\Users\\andre\\OneDrive\\Documentos\\Atuaria\\TeoriaDoRisco\\perdascarteiras\\1frequencias.txt",header=T)
frequencias <- frequencias.entrada[,1]

## Função de verossimilhança da Poisson
vero.pois <- function(dados,param){
  -sum(log(dpois(dados,param)))
}
estima.num <- nlminb(100,vero.pois,dados=frequencias,lower=0)
estima.num$par

library(fitdistrplus)

summary(frequencias)

fpoisMLE <- fitdist(frequencias, "pois", method="mle")

fpoisMLE

summary(fpoisMLE)

