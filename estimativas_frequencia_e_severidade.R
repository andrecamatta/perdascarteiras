##### Parte 1: Estudo da freqüência

### Lendo os dados de frequência
frequencias.entrada <- read.table("C:\\Users\\andre\\OneDrive\\Documentos\\Atuaria\\TeoriaDoRisco\\perdascarteiras\\1frequencias.txt",header=T)
frequencias <- frequencias.entrada[,1]

## Função de verossimilhança da Poisson
vero.pois <- function(dados,param){
  -sum(log(dpois(dados,param)))
}
estima.num <- nlminb(100,vero.pois,dados=frequencias,lower=0)
lambda<-estima.num$par

library(fitdistrplus)

summary(frequencias)

fpoisMLE <- fitdist(frequencias, "pois", method="mle")

fpoisMLE

summary(fpoisMLE)

par(mfrow=c(1,1))
seq=min(frequencias):max(frequencias)
points(dpois(seq,lambda),col="red",lwd=2)

frequencias.df<-data.frame(frequencias)
frequencias.df["ocorrencias"]<-rep("Ocorrencias",length(frequencias))
library(ggplot2)
ggplot(data = frequencias.df) +
  geom_boxplot(aes(x = ocorrencias, y = frequencias)) +
  coord_flip() +
  theme_bw()


