# Roger's new R script to calculate reproductive isolation


#########################################################################################################

# order:  ARC__ARC, ARC__SAX,SAX__ARC,SAX__SAX

# Wales:
rm(list=ls())
total <- c(91,168,126,246)
mate <- c(26,6,2,16)


fem <- c("A","A","S","S")
mal <- c("A","S","A","S")

matlm <- glm(cbind(mate,total-mate)~fem+mal, family=binomial(link=logit))
summary(matlm)
anova(matlm)
fitted(matlm)*total

# Jmating like calculations

Jmate <- function(tot,mat){
  t <- sum(mat)
  TT <- sum(tot)
  Am <- tot[1]+tot[3]
  Af <- tot[1]+tot[2]
  Bm <- tot[2]+tot[4]
  Bf <- tot[3]+tot[4]
  S <- Am*Af+Am*Bf+Bm*Af+Bm*Bf
  AB <- c(Af*Af,Af*Bm,Bf*Am,Bf*Bm)

  PTI <- (mat*S)/(t*AB)
  PSI <- (mat*t)/((mat[c(1,1,1,3)]+mat[c(2,2,3,4)])*(mat[c(1,2,3,2)]+mat[c(3,4,4,4)]))
  PSS <- ((mat[c(1,1,1,3)]+mat[c(2,2,3,4)])*(mat[c(1,2,3,2)]+mat[c(3,4,4,4)])*S)/(t^2*AB)

  I_psi <- (PSI[1]+PSI[4]-PSI[2]-PSI[3])/sum(PSI)
  
  return(l = list(PSI,PSS,I_psi))
}

Jmate(total,mate)

#downsample to minimum of total for all combinations
# and bootstrap with n_boot repeats

AA <- c(rep(0,total[1]-mate[1]),rep(1,mate[1]))
AS <- c(rep(0,total[2]-mate[2]),rep(1,mate[2]))
SA <- c(rep(0,total[3]-mate[3]),rep(1,mate[3]))
SS <- c(rep(0,total[4]-mate[4]),rep(1,mate[4]))


n_boot <- 100

Ipsi_ds <- rep(NA,n_boot)
PSSboot <- matrix(NA,n_boot,4)
PSIboot <- matrix(NA,n_boot,4)

for (i in 1:n_boot){
  AAr <- AA[trunc(runif(min(total),max=length(AA))+1)]
  ASr <- AS[trunc(runif(min(total),max=length(AS))+1)]
  SAr <- SA[trunc(runif(min(total),max=length(SA))+1)]
  SSr <- SS[trunc(runif(min(total),max=length(SS))+1)]
  totds <- rep(min(total),4)
  matds <- c(sum(AAr),sum(ASr),sum(SAr),sum(SSr))
  
  JM <- Jmate(totds,matds)
  Ipsi_ds[i] <- JM[[3]]
  PSSboot[i,] <- JM[[2]]
  PSIboot[i,] <- JM[[1]]
}

hist(Ipsi_ds)

print(c("Mean I_psi",mean(Ipsi_ds)),quote=F)
print(c("CI I_psi",quantile(Ipsi_ds,probs=c(0.025,0.975))),quote=F)

print("Mean PSI [AA,AS,SA,SS; female:male]",quote=F)
colMeans(PSIboot)
print("CI_2.5%",quote=F)
c(quantile(PSIboot[,1],prob=0.025,names=F),quantile(PSIboot[,2],prob=0.025,names=F),quantile(PSIboot[,3],prob=0.025,names=F),quantile(PSIboot[,4],prob=0.025,names=F))
print("CI_97.5%",quote=F)
c(quantile(PSIboot[,1],prob=0.975,names=F),quantile(PSIboot[,2],prob=0.975,names=F),quantile(PSIboot[,3],prob=0.975,names=F),quantile(PSIboot[,4],prob=0.975,names=F))

print("Mean PSS [AA,AS,SA,SS; female:male]",quote=F)
colMeans(PSSboot)
print("CI_2.5%",quote=F)
c(quantile(PSSboot[,1],prob=0.025,names=F),quantile(PSSboot[,2],prob=0.025,names=F),quantile(PSSboot[,3],prob=0.025,names=F),quantile(PSSboot[,4],prob=0.025,names=F))
print("CI_97.5%",quote=F)
c(quantile(PSSboot[,1],prob=0.975,names=F),quantile(PSSboot[,2],prob=0.975,names=F),quantile(PSSboot[,3],prob=0.975,names=F),quantile(PSSboot[,4],prob=0.975,names=F))


###############################################################################################################

# England:
rm(list=ls())

total <- c(128,138,159,323)
mate <- c(48,7,12,45)

fem <- c("A","A","S","S")
mal <- c("A","S","A","S")

matlm <- glm(cbind(mate,total-mate)~fem+mal, family=binomial(link=logit))
summary(matlm)
anova(matlm)
fitted(matlm)*total

# Jmating like calculations

Jmate <- function(tot,mat){
  t <- sum(mat)
  TT <- sum(tot)
  Am <- tot[1]+tot[3]
  Af <- tot[1]+tot[2]
  Bm <- tot[2]+tot[4]
  Bf <- tot[3]+tot[4]
  S <- Am*Af+Am*Bf+Bm*Af+Bm*Bf
  AB <- c(Af*Af,Af*Bm,Bf*Am,Bf*Bm)
  
  PTI <- (mat*S)/(t*AB)
  PSI <- (mat*t)/((mat[c(1,1,1,3)]+mat[c(2,2,3,4)])*(mat[c(1,2,3,2)]+mat[c(3,4,4,4)]))
  PSS <- ((mat[c(1,1,1,3)]+mat[c(2,2,3,4)])*(mat[c(1,2,3,2)]+mat[c(3,4,4,4)])*S)/(t^2*AB)
  
  I_psi <- (PSI[1]+PSI[4]-PSI[2]-PSI[3])/sum(PSI)
  
  return(l = list(PSI,PSS,I_psi))
}

Jmate(total,mate)

#downsample to minimum of total for all combinations
# and bootstrap with n_boot repeats

AA <- c(rep(0,total[1]-mate[1]),rep(1,mate[1]))
AS <- c(rep(0,total[2]-mate[2]),rep(1,mate[2]))
SA <- c(rep(0,total[3]-mate[3]),rep(1,mate[3]))
SS <- c(rep(0,total[4]-mate[4]),rep(1,mate[4]))


n_boot <- 100

Ipsi_ds <- rep(NA,n_boot)
PSSboot <- matrix(NA,n_boot,4)
PSIboot <- matrix(NA,n_boot,4)

for (i in 1:n_boot){
  AAr <- AA[trunc(runif(min(total),max=length(AA))+1)]
  ASr <- AS[trunc(runif(min(total),max=length(AS))+1)]
  SAr <- SA[trunc(runif(min(total),max=length(SA))+1)]
  SSr <- SS[trunc(runif(min(total),max=length(SS))+1)]
  totds <- rep(min(total),4)
  matds <- c(sum(AAr),sum(ASr),sum(SAr),sum(SSr))
  
  JM <- Jmate(totds,matds)
  Ipsi_ds[i] <- JM[[3]]
  PSSboot[i,] <- JM[[2]]
  PSIboot[i,] <- JM[[1]]
}

hist(Ipsi_ds)

print(c("Mean I_psi",mean(Ipsi_ds)),quote=F)
print(c("CI I_psi",quantile(Ipsi_ds,probs=c(0.025,0.975))),quote=F)

print("Mean PSI [AA,AS,SA,SS; female:male]",quote=F)
colMeans(PSIboot)
print("CI_2.5%",quote=F)
c(quantile(PSIboot[,1],prob=0.025,names=F),quantile(PSIboot[,2],prob=0.025,names=F),quantile(PSIboot[,3],prob=0.025,names=F),quantile(PSIboot[,4],prob=0.025,names=F))
print("CI_97.5%",quote=F)
c(quantile(PSIboot[,1],prob=0.975,names=F),quantile(PSIboot[,2],prob=0.975,names=F),quantile(PSIboot[,3],prob=0.975,names=F),quantile(PSIboot[,4],prob=0.975,names=F))

print("Mean PSS [AA,AS,SA,SS; female:male]",quote=F)
colMeans(PSSboot)
print("CI_2.5%",quote=F)
c(quantile(PSSboot[,1],prob=0.025,names=F),quantile(PSSboot[,2],prob=0.025,names=F),quantile(PSSboot[,3],prob=0.025,names=F),quantile(PSSboot[,4],prob=0.025,names=F))
print("CI_97.5%",quote=F)
c(quantile(PSSboot[,1],prob=0.975,names=F),quantile(PSSboot[,2],prob=0.975,names=F),quantile(PSSboot[,3],prob=0.975,names=F),quantile(PSSboot[,4],prob=0.975,names=F))


#########################################################################################################

# England and Wales combained

rm(list=ls())

total <- c(219, 306, 285, 569)
mate <- c(74, 13, 14, 61)

fem <- c("A","A","S","S")
mal <- c("A","S","A","S")

matlm <- glm(cbind(mate,total-mate)~fem+mal, family=binomial(link=logit))
summary(matlm)
anova(matlm)
fitted(matlm)*total

# Jmating like calculations

Jmate <- function(tot,mat){
  t <- sum(mat)
  TT <- sum(tot)
  Am <- tot[1]+tot[3]
  Af <- tot[1]+tot[2]
  Bm <- tot[2]+tot[4]
  Bf <- tot[3]+tot[4]
  S <- Am*Af+Am*Bf+Bm*Af+Bm*Bf
  AB <- c(Af*Af,Af*Bm,Bf*Am,Bf*Bm)
  
  PTI <- (mat*S)/(t*AB)
  PSI <- (mat*t)/((mat[c(1,1,1,3)]+mat[c(2,2,3,4)])*(mat[c(1,2,3,2)]+mat[c(3,4,4,4)]))
  PSS <- ((mat[c(1,1,1,3)]+mat[c(2,2,3,4)])*(mat[c(1,2,3,2)]+mat[c(3,4,4,4)])*S)/(t^2*AB)
  
  I_psi <- (PSI[1]+PSI[4]-PSI[2]-PSI[3])/sum(PSI)
  
  return(l = list(PSI,PSS,I_psi))
}

Jmate(total,mate)

#downsample to minimum of total for all combinations
# and bootstrap with n_boot repeats

AA <- c(rep(0,total[1]-mate[1]),rep(1,mate[1]))
AS <- c(rep(0,total[2]-mate[2]),rep(1,mate[2]))
SA <- c(rep(0,total[3]-mate[3]),rep(1,mate[3]))
SS <- c(rep(0,total[4]-mate[4]),rep(1,mate[4]))


n_boot <- 100

Ipsi_ds <- rep(NA,n_boot)
PSSboot <- matrix(NA,n_boot,4)
PSIboot <- matrix(NA,n_boot,4)

for (i in 1:n_boot){
  AAr <- AA[trunc(runif(min(total),max=length(AA))+1)]
  ASr <- AS[trunc(runif(min(total),max=length(AS))+1)]
  SAr <- SA[trunc(runif(min(total),max=length(SA))+1)]
  SSr <- SS[trunc(runif(min(total),max=length(SS))+1)]
  totds <- rep(min(total),4)
  matds <- c(sum(AAr),sum(ASr),sum(SAr),sum(SSr))
  
  JM <- Jmate(totds,matds)
  Ipsi_ds[i] <- JM[[3]]
  PSSboot[i,] <- JM[[2]]
  PSIboot[i,] <- JM[[1]]
}

hist(Ipsi_ds)

print(c("Mean I_psi",mean(Ipsi_ds)),quote=F)
print(c("CI I_psi",quantile(Ipsi_ds,probs=c(0.025,0.975))),quote=F)

print("Mean PSI [AA,AS,SA,SS; female:male]",quote=F)
colMeans(PSIboot)
print("CI_2.5%",quote=F)
c(quantile(PSIboot[,1],prob=0.025,names=F),quantile(PSIboot[,2],prob=0.025,names=F),quantile(PSIboot[,3],prob=0.025,names=F),quantile(PSIboot[,4],prob=0.025,names=F))
print("CI_97.5%",quote=F)
c(quantile(PSIboot[,1],prob=0.975,names=F),quantile(PSIboot[,2],prob=0.975,names=F),quantile(PSIboot[,3],prob=0.975,names=F),quantile(PSIboot[,4],prob=0.975,names=F))

print("Mean PSS [AA,AS,SA,SS; female:male]",quote=F)
colMeans(PSSboot)
print("CI_2.5%",quote=F)
c(quantile(PSSboot[,1],prob=0.025,names=F),quantile(PSSboot[,2],prob=0.025,names=F),quantile(PSSboot[,3],prob=0.025,names=F),quantile(PSSboot[,4],prob=0.025,names=F))
print("CI_97.5%",quote=F)
c(quantile(PSSboot[,1],prob=0.975,names=F),quantile(PSSboot[,2],prob=0.975,names=F),quantile(PSSboot[,3],prob=0.975,names=F),quantile(PSSboot[,4],prob=0.975,names=F))