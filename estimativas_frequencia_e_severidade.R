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

library(Countr)
breaks_ <- c(190, 200, 210, 220, 230, 240, 250, 260, 270)
count_table(count = frequencias, breaks=breaks_, formatChar=TRUE)
  
  
fnbinomMLE$aic
fpoisMLE$aic



gofstat(list(fnbinomMLE, fpoisMLE), fitnames = c("nbinom", "pois"))


sinistros.arquivo <- "C:\\Users\\andre\\OneDrive\\Documentos\\Atuaria\\TeoriaDoRisco\\perdascarteiras\\1sinistros.txt"
severidades.entrada <- read.table(sinistros.arquivo, header=T)
severidades <- severidades.entrada[,1]

slnormMLE <- fitdist(severidades, "lnorm", method="mle")
summary(slnormMLE)

sgammaMLE <- fitdist(severidades, "gamma", method="mle", lower=c(0,0))
summary(sgammaMLE)

library(actuar)

sparetoMLE <- fitdist(severidades, "pareto", method="mle")
summary(sparetoMLE)

sinvgaussMLE <- fitdist(severidades, "invgauss", method="mle", start = list(mean = mean(severidades), shape = 1))
summary(sinvgaussMLE)

sweibullMLE <- fitdist(severidades, "weibull", method="mle")
summary(sweibullMLE)

sburrMLE <- fitdist(severidades, "burr", method="mle", start = list(shape1=0.9,shape2=1,scale=1))
summary(sburrMLE)



par(mfrow=c(2,3))
hist(severidades,nclass=1000,prob=T, xlim=c(0,20000), main="Ajuste Log-Normal")
curve(dlnorm(x,meanlog=slnormMLE$estimate["meanlog"],sdlog=slnormMLE$estimate["sdlog"]),add=T,col="blue",lwd=2)
hist(severidades,nclass=1000,prob=T, xlim=c(0,20000),main="Ajuste Gamma")
curve(dgamma(x,shape=sgammaMLE$estimate["shape"],rate=sgammaMLE$estimate["rate"]),add=T,col="red",lwd=2)
hist(severidades,nclass=1000,prob=T, xlim=c(0,20000), main="Ajuste Pareto")
curve(dpareto(x,shape=sparetoMLE$estimate["shape"],scale=sparetoMLE$estimate["scale"]),add=T,col="green",lwd=2)
hist(severidades,nclass=1000,prob=T, xlim=c(0,20000), main="Ajuste Gaussiana Inversa" )
curve(dinvgauss(x,mean=sinvgaussMLE$estimate["mean"],shape=sinvgaussMLE$estimate["shape"]),add=T,col="purple",lwd=2)
hist(severidades,nclass=1000,prob=T, xlim=c(0,20000), main="Ajuste Weibull" )
curve(dweibull(x,shape=sweibullMLE$estimate["shape"],scale=sweibullMLE$estimate["scale"]),add=T,col="orange",lwd=2)
hist(severidades,nclass=1000,prob=T, xlim=c(0,20000), main="Ajuste Burr" )
curve(dburr(x,shape1=sburrMLE$estimate["shape1"],shape2=sburrMLE$estimate["shape2"],scale=sburrMLE$estimate["scale"]),add=T,col="pink",lwd=2)


par(mfrow = c(1, 1))
plot.legend <- c("invgauss", "lognormal", "burr")
cdfcomp(list(sinvgaussMLE, slnormMLE, sburrMLE), legendtext = plot.legend, xlogscale=TRUE)


quantile(severidades,0.99)
quantile(slnormMLE,0.99)
quantile(sgammaMLE,0.99)
quantile(sparetoMLE,0.99)
quantile(sinvgaussMLE,0.99)
quantile(sweibullMLE,0.99)
quantile(sburrMLE,0.99)

g<-gofstat(list(slnormMLE, sgammaMLE, sparetoMLE, sinvgaussMLE, sweibullMLE, sburrMLE), fitnames = c("lnorm", "gamma", "pareto", "invgauss", "weibull", "burr"))

bootdist.sinvgaussMLE<-bootdist(sinvgaussMLE, niter=1001)
summary(bootdist.sinvgaussMLE)
plot(bootdist.sinvgaussMLE)


sinvgaussMLE.AD2L <- fitdist(severidades, "invgauss", method="mge", gof="ADL", start = list(mean = mean(severidades), shape = 1))
summary(sinvgaussMLE.AD2L)

par(mfrow = c(1, 1))
plot.legend <- c("invgauss", "lognormal", "burr", "invgauss.AD2L")
cdfcomp(list(sinvgaussMLE, slnormMLE, sburrMLE,sinvgaussMLE.AD2L), legendtext = plot.legend, xlogscale=TRUE)


denscomp(list(slnormMLE, sgammaMLE, sparetoMLE, sinvgaussMLE, sweibullMLE, sburrMLE), legendtext = c("lnorm", "gamma", "pareto", "invgauss", "weibull", "burr"))

g$kstest


ks.test(severidades,"plnorm",meanlog=slnormMLE$estimate["meanlog"],sdlog=slnormMLE$estimate["sdlog"])

Var(X)=(E[X]^2)

