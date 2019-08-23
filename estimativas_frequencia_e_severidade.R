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
frequencias.df["sinistros"]<-rep("Sinistros",length(frequencias))
library(ggplot2)
frequencias.boxplot<-ggplot(data = frequencias.df) +
  geom_boxplot(aes(x = sinistros, y = frequencias)) +
  coord_flip() 

poisson.iqr<-qpois(0.75,lambda)-qpois(0.25,lambda)

frequencias.boxplot + 
  geom_hline(yintercept = lambda, color="red", size=1) + 
  geom_hline(yintercept = qpois(0.75,lambda), linetype="dashed", color="red", size=1) +
  geom_hline(yintercept = qpois(0.25,lambda), linetype="dashed", color="red", size=1) +
  geom_hline(yintercept = qpois(0.75,lambda), linetype="dashed", color="red", size=1) +
  geom_hline(yintercept = qpois(0.25,lambda), linetype="dashed", color="red", size=1) +
  geom_hline(yintercept = qpois(0.75,lambda)+1.5*poisson.iqr, linetype="dotdash", color="red", size=1) +
  geom_hline(yintercept = qpois(0.25,lambda)-1.5*poisson.iqr, linetype="dotdash", color="red", size=1)

vero.nbinom <- function(dados,param){
  -sum(log(dnbinom(dados,size=param[1],prob=param[2])))
}
estima.num.nbinom <- nlminb(c(1000,0.5),vero.nbinom,dados=frequencias,lower=c(0,0),upper=c(Inf,1))
estima.num.nbinom$par

fnbinomMLE <- fitdist(frequencias, "nbinom", method="mle")

fnbinomMLE

summary(fnbinomMLE)


prob.nbinomMLE <- fnbinomMLE$estimate["size"]/(fnbinomMLE$estimate["mu"]+fnbinomMLE$estimate["size"])
prob.nbinomMLE

nbinom.size<-fnbinomMLE$estimate["size"]
nbinom.prob<-prob.nbinomMLE
nbinom.iqr<-qnbinom(0.75,nbinom.size,nbinom.prob)-qnbinom(0.25,nbinom.size,nbinom.prob)

  
  




 



