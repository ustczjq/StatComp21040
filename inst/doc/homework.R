## -----------------------------------------------------------------------------
set.seed(7725)
U <- runif(1000,0,1)          # U ~ U[0,1]
V <- runif(1000,0,1)          # v ~ U[0,1]

## -----------------------------------------------------------------------------
m<-sqrt(-2*log(U))*cos(2*pi*V)
n<-sqrt(-2*log(U))*sin(2*pi*V)
par(mfrow=c(1,2))
hist(m,main = "Histogram of Z1")
hist(n,main = "Histogram of Z2")


## -----------------------------------------------------------------------------
r<- rnorm(1000,0,1)
par(mfrow=c(1,2))
qqplot(r,m,ann="FALSE")
mtext("Normal Q-Q plots for U",adj=0.5,line=1)
qqplot(r,n,ann="FALSE")
mtext("Normal Q-Q plots for V",adj=0.5,line=1)


## -----------------------------------------------------------------------------
library(nortest)
cvm.test(m)
cvm.test(n)

## -----------------------------------------------------------------------------
knitr::kable(head(trees))     

## -----------------------------------------------------------------------------
lmtrees <- lm(Volume~Girth+Height,data=trees)
summary(lmtrees)$coef
 

## -----------------------------------------------------------------------------
plot(lmtrees)  
 

## -----------------------------------------------------------------------------
set.seed(2021)

U <- runif(1000,0,1)        
C1<-1                                 # sigma = 1
x<- sqrt(-2*C1^2*log(1-U)) 
hist(x,probability=TRUE,main=expression(f(x)==x%.%exp(-frac(x^2,2))))
z1<-seq(0,4,0.01)
lines(z1,z1*exp(-z1^2/2),col="red")

C2<-5                                 # sigma = 5
y<- sqrt(-2*C2^2*log(1-U)) 
hist(y,probability=TRUE,main=expression(f(x)==frac(x,25)%.%exp(-frac(x^2,50))))
z2<-seq(0,20,0.05)
lines(z2,z2/25*exp(-z2^2/50),col="red")


C3<-15                                # sigma = 15
w<- sqrt(-2*C3^2*log(1-U)) 
hist(w,probability=TRUE,main=expression(f(x)==frac(x,225)%.%exp(-frac(x^2,500))))
z3<-seq(0,55,0.1)
lines(z3,z3/225*exp(-z3^2/500),col="red")

## -----------------------------------------------------------------------------
set.seed(2021)
n<-1000
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)

#Use sample method to simulate p1=0.75
r1<- sample(c(1,0),n,replace = TRUE,prob = c(0.75,0.25))
Z1<-r1*X1+(1-r1)*X2
hist(Z1,prob=TRUE,xlim = c(-5,9),ylim = c(0,0.28),main="p1=0.75,p2=0.25")
lines(density(Z1),lwd = 2, col = "blue")


## -----------------------------------------------------------------------------
r2<- sample(c(1,0),n,replace = TRUE,prob = c(0.25,0.75))
Z2<-r2*X1+(1-r2)*X2
hist(Z2,prob=TRUE,xlim =  c(-5,9),ylim = c(0,0.3),main="p1=0.25,p2=0.75")
lines(density(Z2),lwd = 2, col = "blue")

## -----------------------------------------------------------------------------

r3<- sample(c(1,0),n,replace = TRUE,prob = c(0.5,0.5))
Z3<-r3*X1+(1-r3)*X2
hist(Z3,prob=TRUE,xlim = c(-5,9),ylim = c(0,0.25),main="p1=0.5,p2=0.5")
lines(density(Z3),lwd = 2, col = "blue")


## -----------------------------------------------------------------------------

# Compound Poisson–Gamma process

CPG<-function(n,t,lamb,alpha,beta){
N<-rpois(n,lamb*t)
Z<-c()
for (i in 1:n){
    Z[i]<-sum(rgamma(N[i],alpha,beta));
}


print(paste("The mean of X is ",mean(Z),",which is very close to the theoretical value ",lamb*t*alpha/beta))
print(paste("The variance of X is ",var(Z),",which is very close to the theoretical value ",lamb*t*((alpha/beta)^2+(alpha/beta^2))))
}


#Different parameter value combinations:
CPG(1000,10,1,2,2)          
CPG(2000,10,3,4,4)
CPG(1500,10,0.5,1,1)

## -----------------------------------------------------------------------------
n<-10000
x<-seq(from=0.1,to=0.9,by=0.1) 
set.seed(2021)
u <- numeric(n)
E<- numeric(9)

for (i in 1:9){
u <- runif(n,0,x[i])
E[i]<- round(mean(x[i]*30*u^2*(1-u)^2),5)   #
}

F<-pbeta(x,3,3)
print(paste("The Monte Carlo estimate of Beta(3,3) cdf F(x) for x=",x," is ",E,",while the value returned by the pbeta function in R is ",F))


## -----------------------------------------------------------------------------
set.seed(2021)
MCR<- function(s=1,N=1000,antithetic=TRUE) {   # s means sigma
  u <- runif(N/2,0,1)              # half random numbers from U[0,1]
  if (!antithetic) 
    v <- runif(N/2,0,1)     # antithetic=False,i.e.the original situation
  else 
    v <- 1 - u
  
  u <- c(u, v)
  X<-mean(s*sqrt(-2*log(1-u)))
  X
}


m <- 1000
MCR1<- MCR2<- numeric(m)


for (i in 1:m) {
  MCR1[i]<-MCR(s=2,N=1000,anti=FALSE)    #Original method
  MCR2[i]<-MCR(s=2,N=1000,anti=TRUE)     #Antithetic method
}

VarianceRatio<-round(var(MCR2)/var(MCR1),4)
print(paste("The ratio of variance between two methods is",VarianceRatio,"and the percent reduction in variance is about",round((var(MCR1)-var(MCR2))/var(MCR1),4)*100,"%"))
  


## -----------------------------------------------------------------------------
 g<-function(x) 
   {x^2/sqrt(2*pi)*exp(-x^2/2)*(x>=1)}  
 f1<-function(x) 
   {exp(-x+1)}     
 f2<-function(x) 
 {4/pi/(1+x^2)}          

 
 x<-seq(1,5,0.01)
gf1<-c(expression(g(x)==frac(x^2,sqrt(2*pi))*exp(-frac(x^2,2))),expression(f1(x)==exp(-(x-1))),expression(f2(x)==frac(4,pi*(1+x^2))))

 plot(x,g(x), type = "l", ylab = "",ylim = c(0,1), lwd = 0.3,col="red",main='The plots of g and f')
 lines(x, f1(x), lty = 2, lwd = 0.3,col="blue")
 lines(x, f2(x), lty = 3, lwd = 0.3,col="green")
 legend("topright", legend = gf1, lty = 1:3, col=c("red","blue","green"),lwd = 0.2,cex=1.2)

## -----------------------------------------------------------------------------
gf2<-c(expression(g(x)/f1(x)),expression(g(x)/f2(x)))

 plot(x,g(x)/f1(x), type = "l", ylab = "",ylim = c(0,1), lwd = 0.3,col="blue",main='The plots of g/f')
 lines(x, g(x)/f2(x), lty = 2, lwd = 0.3,col="green")
 legend("topright", legend = gf2, lty = 1:2, col=c("blue","green"),lwd = 0.2,cex=1.2)

## -----------------------------------------------------------------------------
set.seed(100)
N<-10000
u<-runif(N,0,1)
x<-1-log(1-u)
y<-tan(pi/4*(u+1))
integral1<-mean(g(x)/f1(x))
integral2<-mean(g(y)/f2(y))
Real<-integrate(g,1,Inf)
print(paste("The estimate of the integral using f1 is ",round(integral1,5),"and the real value is ",round(Real$value,5)))
print(paste("The estimate of the integral using f2 is ",round(integral2,5),"and the real value is ",round(Real$value,5)))
print(paste("The standard variance of estimation using f1 is ",round(sd(g(x)/f1(x)),5),"and the one using f2 is ",round(sd(g(x)/f2(x)),5)))

## -----------------------------------------------------------------------------
 g<-function(x) 
   {x^2/sqrt(2*pi)*exp(-x^2/2)}
 f3<-function(x) 
   {exp(-x)/(1-exp(-1))}
 f4<-function(x) 
 {x/(1-exp(1)^(-1/2))*exp(-x^2/2)}
 
  x<-seq(0,1,0.001)
  
 gf3<-c(expression(g(x)==frac(x^2,sqrt(2*pi))*exp(-frac(x^2,2))),expression(f1(x)==frac(exp(-x),(1-exp(-1)))),expression(f2(x)==frac(x*exp(-frac(x^2,2)),(1-e^(-frac(1,2))))))

 plot(x,g(x), type = "l", ylab = "",ylim = c(0,1.5), lwd = 0.3,col="red",main='The plots of g and f')
 lines(x, f3(x), lty = 2, lwd = 0.3,col="blue")
 lines(x, f4(x), lty = 3, lwd = 0.3,col="green")
 legend("topleft", legend = gf3, lty = 1:3, col=c("red","blue","green"),lwd = 0.2,cex=0.9)

## -----------------------------------------------------------------------------
gf4<-c(expression(g(x)/f1(x)),expression(g(x)/f2(x)))

 plot(x,g(x)/f3(x), type = "l", ylab = "",ylim = c(0,0.2), lwd = 0.3,col="blue",main='The plots of g/f')
 lines(x, g(x)/f4(x), lty = 2, lwd = 0.3,col="green")
 legend("topleft", legend = gf4, lty = 1:2, col=c("blue","green"),lwd = 0.2,cex=1.2)

## -----------------------------------------------------------------------------
set.seed(2021)
N<-10000
u<-runif(N,0,1)
z<--log(1-u*(1-exp(-1)))
integral3<-1/2-mean(g(z)/f3(z))

x<-sqrt(-2*log(1-u))         #Generate Rayleigh random numbers 
y<-x[x<=1]                   #Select the random numbers locating on [0,1]
integral4<-1/2-(1-exp(-1/2))/sqrt(2*pi)*mean(y)

print(paste("The estimate of the integral using f1 is ",round(integral3,5),"which is close to the real value",round(Real$value,5)))
print(paste("The estimate of the integral using f2 is ",round(integral4,5),"which is close to the real value",round(Real$value,5)))

## -----------------------------------------------------------------------------
#t-interval
set.seed(2021)
n<-20
alpha <- 0.05
N<-1e3
Lower<-Upper<-numeric(N)

for(i in 1:N)
{
    x<-rchisq(n, 2) 
    Lower[i]=mean(x)-qt(alpha/2,df=n-1,lower.tail = FALSE)*sd(x)/sqrt(n) 
    #lower.tail = FALSE calculates the upper quantile  
    Upper[i]=mean(x)+qt(alpha/2,df=n-1,lower.tail = FALSE)*sd(x)/sqrt(n)
}

print(paste("The estimation of the coverage probability of the t-interval is",mean(Lower<2 & Upper>2),"while the nominal level is ",1-alpha))  

## -----------------------------------------------------------------------------
#variance interval
set.seed(2021)
n<-20
alpha <- 0.05
Upper2<-numeric(N)
for(i in 1:N)
{
    x <- rchisq(n, 2) 
    Upper2[i] <- (n-1) * var(x) / qchisq(alpha, df=n-1)
}
print(paste("The estimation of the coverage probability of variance interval is",mean(Upper2>4),"while the nominal level is ",1-alpha))  

## -----------------------------------------------------------------------------
set.seed(2021)
N<-10000
num<-c(20,50,250,500)
alpha<-0.05  
mu<-1
M<-numeric(N)

results <- matrix(0,nrow = length(num),ncol = 4)
colnames(results) <- c("sample size","chi square(1)","U[0,2]","Exp(1)")

for (n in num){
  quant<-qt(1-alpha/2,n-1)
for (i in 1:N) {
 X<-rchisq(n,1)
 T<-abs(sqrt(n)*(mean(X)-mu)/sd(X))
 M[i]<-as.integer(T>quant)
 
}
results[which(num == n),2]<-mean(M)


for (i in 1:N) {
 X<-runif(n,0,2)
T<-abs(sqrt(n)*(mean(X)-mu)/sd(X))
 M[i]<-as.integer(T>quant)
}
results[which(num == n),3]<-mean(M)

for (i in 1:N) {
 X<-rexp(n,1)
 T<-abs(sqrt(n)*(mean(X)-mu)/sd(X))
 M[i]<-as.integer(T>quant)
}
results[which(num == n),4]<-mean(M)
}

for (p in 1:length(num)){
  results[p,1]<-num[p]
  }


library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
set.seed(2021)
N<-10000
n<-100
alpha<-0.05
p1<-0.651
p2<-0.676
quant<-qt(1-alpha/2,n-1)
M<-numeric(N)
for (i in 1:N) {
  x<-rbinom(n,1,p1)
  y<-rbinom(n,1,p2)
  T<-(mean(x-y))/sqrt((var(x-y))/n)
  M[i]<-as.integer(T>quant)
}
print(paste("the p-value is",2*mean(M),"while the significance level is",alpha))

## -----------------------------------------------------------------------------
library(MASS)
set.seed(2021)
d<-3                       #dimension 
num<-c(10,50,500)          #sample size
m<-500                    #simulation times
alpha<-0.05                #significance level
R <- diag(rep(1,d))        #required in function mvnorm() as Sigma 
skk<- numeric(m) 
criticalvalue<-qchisq(1-alpha,d*(d+1)*(d+2)/6)

#function to get b1d.
skew<-function(X) {
Xbar<-colMeans(X)
Sigma<-cov(X)
b<-numeric(n)
  for (i in 1:n) {
    b[i]<-mean(((X[i,]-Xbar)%*%ginv(Sigma)%*%t((X-Xbar)))^3)
  }

return(mean(b)*n/6)

}
#The first loop is to produce results for different sample size,and the second is used for m simulations.
for (n in num){
for (k in 1:m) {
  X<-mvrnorm(n,mu=rep(0, d),Sigma = R)
  
  #Or generate X,Y in this way:
  #X<-matrix(rnorm(n*d),nrow=n,ncol=d)
 
  
  skk[k] <- as.integer(skew(X)>=criticalvalue)
}
p<-mean(skk)
print(paste("The critical probability for n=",n,"is",p))
}


## -----------------------------------------------------------------------------
library(MASS)
set.seed(2021)
alpha<-0.05    #significance level
n<-30          #sample size  
m<-300         #simulation times
d<-3           #dimension 
criticalvalue2<-qchisq(1-alpha,d*(d+1)*(d+2)/6)
eps<-c(seq(0,0.15,0.01), seq(0.15,1,0.05))
N<-length(eps)
power<-numeric(N)    
skk<-numeric(m)

#function to get b1d.
skew<-function(X) {
Xbar<-colMeans(X)
Sigma<-cov(X)
b<-numeric(n)
  for (i in 1:n) {
    b[i]<-mean(((X[i,]-Xbar)%*%ginv(Sigma)%*%t((X-Xbar)))^3)
  }
return(mean(b)*n/6)
}


#for each epsilon,produce the contaminated probability.
for (j in 1:N) { 
e <- eps[j]

#The loop is produced for m simulations.
for (i in 1:m) { 

 sigma<-sample(c(1, 10), replace = TRUE,size = n*d, prob = c(1-e, e))
 X<-matrix(rnorm(n*d,0,sigma),nrow=n,ncol=d)
 skk[i]<-as.integer(skew(X)>=criticalvalue2)
}
 power[j] <- mean(skk)
}


#plot power vs epsilon
plot(eps, power, type = "b",
xlab = bquote(epsilon), ylim = c(0,1))
abline(h=alpha,lty =3,col="red")
se<-sqrt(power*(1-power)/m)     
lines(eps,power+se,lty=3)
lines(eps,power-se,lty=3)

## -----------------------------------------------------------------------------
set.seed(2021)
library(bootstrap)     #obtain scor data
library(boot)          #use function boot()
lambdaesti<-eigen(cov(scor))$values
thetaesti<-lambdaesti[1]/sum(lambdaesti)
N<-400                #number of bootstrap samples
SCOR<-cbind(scor$mec,scor$vec,scor$alg,scor$ana,scor$sta)


#Bootstrap
bs<-function(data,index){
x<-data[index,]
lambda<-eigen(cov(x))$values
theta<-lambda[1]/sum(lambda)
return(theta)
}

bsr<-boot(data=SCOR,statistic=bs,R=N)

thetabs<-bsr$t
biasbs<-mean(thetabs)-thetaesti
sebs<-sd(thetabs)


results<-matrix(0,nrow=1,ncol=2)
rownames(results)<-c("Bootstrap")
colnames(results)<-c("bias","standard error")
results[1,1]<-biasbs
results[1,2]<-sebs

library(knitr)
knitr::kable(results,align=rep('c', 5))


## -----------------------------------------------------------------------------
# Jackknife
n<-nrow(scor)           #sample size
thetajk<-rep(0,n)
for (i in 1:n) {
x<-scor[-i,]
lambda<-eigen(cov(x))$values
thetajk[i]<-lambda[1]/sum(lambda)
}

biasjk<-(n-1)*(mean(thetajk)-thetaesti)
sejk<-(n-1)*sqrt(var(thetajk)/n)


results2<-matrix(0,nrow=1,ncol=2)
rownames(results2)<-c("Jackknife")
colnames(results2)<-c("bias","standard error")
results2[1,1]<-biasjk
results2[1,2]<-sejk


library(knitr)
knitr::kable(results2,align=rep('c',5))


## -----------------------------------------------------------------------------
set.seed(2021)
bs<-function(data,index){
x<-data[index,]
lambda<-eigen(cov(x))$values
theta<-lambda[1]/sum(lambda)
return(theta)
}

bsr<-boot(data=SCOR,statistic=bs,R=1000)

bsCI<-boot.ci(bsr,conf=0.95,type=c("perc", "bca"))
print(bsCI)

## -----------------------------------------------------------------------------
library(boot)
set.seed(2021)
sk<-function(data,index){
x<-data[index]
xbar<-mean(x)
e3<-mean((x-xbar)^3)
e2<-mean((x-xbar)^2)
return(e3/(e2^1.5))
}
n<-50

s1<-0
s2<-4/sqrt(10)     #THe skewness of χ2(5) distributions
B<-500

normCI1<-basicCI1<-percentCI1<-matrix(0,B,2)
normCI2<-basicCI2<-percentCI2<-matrix(0,B,2)
for (i in 1:B){
  X<-rnorm(n)
  Y<-rchisq(n,5)
  boot.obj1<-boot(data=X,statistic=sk,R=500)
  boot.obj2<-boot(data=Y,statistic=sk,R=500)
  CI1<-boot.ci(boot.obj1,conf=0.95,type =c("norm","basic","perc"))
  CI2<-boot.ci(boot.obj2,conf=0.95,type=c("norm","basic","perc"))
normCI1[i,]<-CI1$norm[2:3]
basicCI1[i,]<-CI1$basic[4:5]
percentCI1[i,]<-CI1$percent[4:5]

normCI2[i,]<-CI2$norm[2:3]
basicCI2[i,]<-CI2$basic[4:5]
percentCI2[i,]<-CI2$percent[4:5]
}


## -----------------------------------------------------------------------------
results3<-matrix(0,nrow=2,ncol=3)
rownames(results3)<-c("normal","chi-suqared(5)")
colnames(results3)<-c("standard normal ","basic","percentile")
results3[1,1]<-mean(normCI1[,1]<=s1 & normCI1[,2]>=s1)
results3[1,2]<-mean(basicCI1[,1]<=s1 & basicCI1[,2]>=s1)
results3[1,3]<-mean(percentCI1[,1]<=s1 & percentCI1[,2]>=s1)
results3[2,1]<-mean(normCI2[,1]<=s2 & normCI2[,2]>=s2)
results3[2,2]<-mean(basicCI2[,1]<=s2 & basicCI2[,2]>=s2)
results3[2,3]<-mean(percentCI2[,1]<=s2 & percentCI2[,2]>=s2)

## -----------------------------------------------------------------------------
library(knitr)
knitr::kable(results3,align=rep('c',5))

## -----------------------------------------------------------------------------

results4<-matrix(0,nrow=2,ncol=3)
rownames(results4)<-c("normal","chi-suqared(5)")
colnames(results4)<-c("standard normal","basic","percentile")
results4[1,1]<-mean(normCI1[,1]>s1)
results4[1,2]<-mean(basicCI1[,1]>s1)
results4[1,3]<-mean(percentCI1[,1]>s1)
results4[2,1]<-mean(normCI2[,1]>s2)
results4[2,2]<-mean(basicCI2[,1]>s2)
results4[2,3]<-mean(percentCI2[,1]>s2)

## -----------------------------------------------------------------------------
knitr::kable(results4,align=rep('c',5))

## -----------------------------------------------------------------------------

results5<-matrix(0,nrow=2,ncol=3)
rownames(results5)<-c("normal","chi-suqared(5)")
colnames(results5)<-c("standard normal bootstrap CI ","basic bootstrap CI","percentile bootstrap CI")
results5[1,1]<-mean(normCI1[,2]<s1)
results5[1,2]<-mean(basicCI1[,2]<s1)
results5[1,3]<-mean(percentCI1[,2]<s1)
results5[2,1]<-mean(normCI2[,2]<s2)
results5[2,2]<-mean(basicCI2[,2]<s2)
results5[2,3]<-mean(percentCI2[,2]<s2)

## -----------------------------------------------------------------------------
knitr::kable(results5,align=rep('c',5))

## -----------------------------------------------------------------------------
set.seed(2021)
n<-200
X<-runif(n,0,2)
Y<-rexp(n,1)


independence<-function(z,index) {
#dims contains dimensions of x and y
x<-z[,1]       
y<-z[index,2] 
return(cor(x,y,method="spearman"))
}


library(boot)
z<-matrix(cbind(X,Y),ncol=2)
boot.obj<-boot(data=z,statistic=independence,R=499,sim = "permutation")
tb<-c(boot.obj$t0, boot.obj$t)

print(paste("The proportion of insignificant permutation is",mean(abs(tb)>=abs(boot.obj$t0))))
print(paste("The p-value obtained by function cor.test() is",round(cor.test(X,Y,alternative="two.sided",method="spearman")$p.value,3)))


## -----------------------------------------------------------------------------
set.seed(2021)
n<-200
X<-runif(n,0,2)
Z<-2-X


z<-matrix(cbind(X,Z),ncol=2)
boot.obj<-boot(data=z,statistic=independence,R=499,sim = "permutation")
tb<-c(boot.obj$t0, boot.obj$t)

print(paste("The proportion of insignificant permutation is",mean(abs(tb)>=abs(boot.obj$t0))))
print(paste("The p-value obtained by function cor.test() is",round(cor.test(X,Z,alternative="two.sided",method="spearman")$p.value,3)))


## -----------------------------------------------------------------------------
m<-500               #Permutation times
d<-2                 #The dimension of data
n1<-n2<-30           #The sample size of X and Y
R<-499               #The number of bootstrap replicates.
k<-2                 #The Nearest Neighbor parameter
n<-n1+n2             #Total sample size
N=c(n1,n2)
p.values<-p.values2<-p.values3<-p.values4<-matrix(NA,m,3)   #contains p-value for 3 methods.
alpha<-0.05          #Significance level

#NN method
Tn<-function(z,ix,sizes,k){
  n1<-sizes[1]; n2<-sizes[2]; n<-n1+n2
  if(is.vector(z))z<-data.frame(z,0);
  z<-z[ix,];
  NN<-nn2(data=z,k=k+1)
  block1<-NN$nn.idx[1:n1,-1]
  block2<-NN$nn.idx[(n1+1):n,-1]
  i1<-sum(block1<n1+.5)
  i2<-sum(block2>n1+.5)
  return((i1+i2)/(k*n))
}

eqdist.nn<-function(z,sizes,k){
  boot.obj<-boot(data=z,statistic=Tn,R=R,sim="permutation",sizes=sizes,k=k)
  ts<-c(boot.obj$t0,boot.obj$t)
  p.value<-mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

set.seed(2021)
library(RANN)
library(energy)
library(Ball)

#Unequal variances and equal expectations
for(i in 1:m){
  X<-matrix(rnorm(n1*d,mean=0,sd=1),ncol=d)
  Y<-matrix(rnorm(n2*d,mean=0,sd=2),ncol=d)
  z<-rbind(X,Y)
  #NN method
  p.values[i,1]<-eqdist.nn(z,N,k)$p.value
  #energy methods
  p.values[i,2]<-eqdist.etest(z,sizes=N,R=R)$p.value
  #Ball method
  p.values[i,3]<-bd.test(x=X,y=Y,num.permutations=R,seed=i*2021)$p.value
}

power1<-colMeans(p.values<alpha)
results<-matrix(0,nrow=1,ncol=3)
colnames(results)<-c("NN method","energy method","Ball method")
rownames(results)<-"power"
results[1,1]<-power1[1]
results[1,2]<-power1[2]
results[1,3]<-power1[3]


library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
for(i in 1:m){
  X2<-matrix(rnorm(n1*d,mean=0,sd=1),ncol=d)
  Y2<-matrix(rnorm(n2*d,mean=0.5,sd=1.5),ncol=d)
  z2<-rbind(X2,Y2)
  #NN method
  p.values2[i,1]<-eqdist.nn(z2,N,k)$p.value
  #energy methods
  p.values2[i,2]<-eqdist.etest(z2,sizes=N,R=R)$p.value
  #Ball method
  p.values2[i,3]<-bd.test(x=X2,y=Y2,num.permutations=R,seed=i*2021)$p.value
}


power2<-colMeans(p.values2<alpha)
results<-matrix(0,nrow=1,ncol=3)
colnames(results)<-c("NN method","energy method","Ball method")
rownames(results)<-"power"
results[1,1]<-power2[1]
results[1,2]<-power2[2]
results[1,3]<-power2[3]


library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
for(i in 1:m){
  X3<-matrix(rt(n1*d,df=1),ncol=d)
  y1=rnorm(n2*d,mean=0,sd=1);
  y2=rnorm(n2*d,mean=1,sd=2)
  r<-sample(c(1,0),n,replace =TRUE,prob=c(0.5,0.5))
  Y3<-matrix(r*y1+(1-r)*y2,ncol=d)
  z3<-rbind(X3,Y3)
  #NN method
  p.values3[i,1]<-eqdist.nn(z3,N,k)$p.value
  #energy methods
  p.values3[i,2]<-eqdist.etest(z3,sizes=N,R=R)$p.value
  #Ball method
  p.values3[i,3]<-bd.test(x=X3,y=Y3,num.permutations=R,seed=i*2021)$p.value
}


power3<-colMeans(p.values3<alpha)
results<-matrix(0,nrow=1,ncol=3)
colnames(results)<-c("NN method","energy method","Ball method")
rownames(results)<-"power"
results[1,1]<-power3[1]
results[1,2]<-power3[2]
results[1,3]<-power3[3]


library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
N2=c(n1,2*n2)
for(i in 1:m){
  X4<-matrix(rnorm(n1*d,mean=0,sd=1),ncol=d)
  Y4<-matrix(rnorm(2*n2*d,mean=0.5,sd=1),ncol=d)
  z4<-rbind(X4,Y4)

  #NN method
  p.values4[i,1]<-eqdist.nn(z4,N2,k)$p.value
  #energy methods
  p.values4[i,2]<-eqdist.etest(z4,sizes=N2,R=R)$p.value
  #Ball method
  p.values4[i,3]<-bd.test(x=X4,y=Y4,num.permutations=R,seed=i*2021)$p.value
}


power4<-colMeans(p.values4<alpha)
results<-matrix(0,nrow=1,ncol=3)
colnames(results)<-c("NN method","energy method","Ball method")
rownames(results)<-"power"
results[1,1]<-power4[1]
results[1,2]<-power4[2]
results[1,3]<-power4[3]


library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
set.seed(2021)

#The standard Cauchy distribution's density function

f<-function(x) {
  return(1/(pi*(1+x^2)))
}


m <- 1000                     #the number of samples
x <- numeric(m)                #the vector that stores samples
x[1]<-runif(1,-1,1)            #the first sample
k<-0              #k records the number of rejected candidate points.
burntime<-100    #discard the first 1000 samples
u<-runif(m)       #useful for acceptance probability. 


#generate the chain
for (i in 2:m) {
  xt<-x[i-1]
  y<-xt+runif(1,min=-1,max=1)
  num<-f(y)*dnorm(xt,mean=y,sd=1)          #the numerator
  den<-f(xt)*dnorm(y,mean=xt,sd=1)         #the denominator
  if (u[i]<=num/den) x[i]<-y 
  else {
    x[i]<-xt
    k<-k+1               #y is rejected 
  } 
}
print(paste("The rate of candidate points rejected is ",k/m))



## -----------------------------------------------------------------------------
X<-x[burntime+1:m]             #discard the first 1000 samples.
index<-burntime+1:m
plot(index,X,type="l",main="trace plot",ylab="random number")
abline(h=-5,lty =2,col="red")     
#to highlight the range where most samples from t(1) fall in.
abline(h=5,lty =2,col="red")

## -----------------------------------------------------------------------------
hist(X,breaks="scott",main="histogram  of the generated samples",xlab="",freq=FALSE,xlim=c(-10,10),ylim=c(0,0.35))
lines(k,dcauchy(k),col="red",lwd=1.5)

## -----------------------------------------------------------------------------
Q<-quantile(X,probs=seq(0.1,0.9,0.1),na.rm=TRUE)
Qqc<-qcauchy(seq(0.1,0.9,0.1),loc=0,scale=1)
k<-seq(-10,10,0.01)


results<-matrix(0,nrow=9,ncol=2)
rownames(results)<-c("10%","20%","30%","40%","50%","60%","70%","80%","90%")
colnames(results)<-c("Generated observations","standard Cauchy distribution")
results[,1]<-Q
results[,2]<-Qqc

library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
    Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
        }

    Chain<-function(sigma,N,StartPoint) { 
#generates a Metropolis chain for t(1) with N(X[t],sigma) proposal distribution and deterimined starting value. 
        x <- rep(0, N)
        x[1] <-StartPoint
        u <- runif(N)

        for (i in 2:N) {
        xt<-x[i-1]
        y<-xt+runif(1,min=-1,max=1)
        num<-dt(y,1)*dnorm(xt,mean=y,sd=sigma)         
        den<-dt(xt,1)*dnorm(y,mean=xt,sd=sigma)
        s<-num/den
          if (u[i]<=s)    x[i]<-y 
          else            x[i]<-xt
        }
        return(x)
        }

    sigma<-1     #parameter of proposal distribution  N(X[t],sigma)
    k<-4         #number of chains to generate
    n<- 10000      #length of chains
    burntime <- 1000       #burn-in length

    #choose overdispersed initial values
    x0 <- c(-1,1,-4,4)

    #generate the chains
    set.seed(2021)
    X <- matrix(0, nrow=k, ncol=n)
    for (i in 1:k)
        X[i, ] <- Chain(sigma, n, x0[i])

    #compute diagnostic statistics : median
    ind<-1
    psi<-matrix(NA,nrow=k,ncol=n)
    #for a matrix 1 indicates rows, 2 indicates columns
    for (i in 1:k){
       for (ind in 1:n){
         Z<-as.matrix(X[i,1:ind])
     psi[i,ind]<-apply(Z,2,median)
     ind<-ind+1
       }
}
    #plot psi for the four chains
     for (i in 1:k)
      if(i==1){
        plot((burntime+1):n,psi[i, (burntime+1):n],ylim=c(-1,1), type="l",
            xlab='index', ylab=bquote(phi),main="trace plot")
      }else{
        lines((burntime+1):n,psi[i,(burntime+1):n],col=i)
    }

    

## -----------------------------------------------------------------------------
    
    #plot the sequence of R-hat statistics
    rhat<-rep(0,n)
    for (j in (burntime+1):n)
        rhat[j]<-Gelman.Rubin(psi[,1:j])
    plot(rhat[(burntime+1):n], type="l", xlab="",ylab="R",main = "the Gelman-Rubin statistic ")
    abline(h=1.2,lty=2,col="red")

## -----------------------------------------------------------------------------
set.seed(2021)
N<-1000                    #length of chain
burntime<-100             #burn-in length
X<-matrix(0,N,2)            #bivariate sample
a<-1
b<-1
#Beta(1,1) is equivalent to U[0,1]
n<-25 
X[1,]<-c(0,0)              #initialize starting point

for (i in 2:N) {
x2 <- X[i-1, 2]
X[i, 1] <- rbinom(1,size=n,prob=x2)
x1<- X[i,1]
X[i, 2] <- rbeta(1,x1+a,n-x1+b)
}

## -----------------------------------------------------------------------------
x<-X[(burntime+1):N,]
plot(x[,1],main="trace plot for X and Y",type='l',lwd=0.7,xlab='index',ylab='Random numbers',col="red")
lines(x[,2],lwd=0.7,col="blue")
legend('topright',c("X","Y"),lwd=1,col=c("red","blue"),cex=1.2)


## -----------------------------------------------------------------------------
plot(x, main="Bivariate chain generated by the Gibbs sampler",cex=0.5,xlab="X",ylab="Y",ylim=range(x[,2]))

## -----------------------------------------------------------------------------
print(paste("The mean of X and Y is ",round(colMeans(x)[1],4),"and ",round(colMeans(x)[2],4),",respectively."))

## -----------------------------------------------------------------------------
set.seed(2021)
a<-1
b<-1
n<-25
k<-3                     #number of chains to generate
N<- 1000                #length of chains
burntime <- 100         #burn-in length
X1<-X2<-X3<-matrix(0,N,2)         #bivariate sample

 #initialize starting point
X1[1,]<-c(0,0)             
X2[1,]<-c(0.4,0.2)
X3[1,]<-c(0.3,0.6)

 Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
        }

for (i in 2:N) {
x2 <- X1[i-1, 2]
X1[i, 1] <- rbinom(1,size=n,prob=x2)
x1<- X1[i,1]
X1[i, 2] <- rbeta(1,x1+a,n-x1+b)
}

for (i in 2:N) {
x2 <- X2[i-1, 2]
X2[i, 1] <- rbinom(1,size=n,prob=x2)
x1<- X2[i,1]
X2[i, 2] <- rbeta(1,x1+a,n-x1+b)
}

for (i in 2:N) {
x2 <- X3[i-1, 2]
X3[i, 1] <- rbinom(1,size=n,prob=x2)
x1<- X3[i,1]
X3[i, 2] <- rbeta(1,x1+a,n-x1+b)
}

Z1<-matrix(0,nrow=k,ncol=N)
Z2<-matrix(0,nrow=k,ncol=N)

## -----------------------------------------------------------------------------
Z1[1,]<-X1[,1]
Z1[2,]<-X2[,1]
Z1[3,]<-X3[,1]

psi<- t(apply(Z1, 1, cumsum))
    for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))
 
rhat <- rep(0,N)
    for (j in (burntime+1):N)
        rhat[j] <- Gelman.Rubin(psi[,1:j])
    plot(rhat[(burntime+1):N], type="l", xlab="", ylab="R",main="the Gelman-Rubin statistic for X")
    abline(h=1.2, lty=2,col="red")

## -----------------------------------------------------------------------------
Z2[1,]<-X1[,2]
Z2[2,]<-X2[,2]
Z2[3,]<-X3[,2]

psi2<- t(apply(Z2, 1, cumsum))
    for (i in 1:nrow(psi2))
        psi2[i,] <- psi2[i,] / (1:ncol(psi2))
 
rhat2 <- rep(0,N)
    for (j in (burntime+1):N)
        rhat2[j] <- Gelman.Rubin(psi2[,1:j])
    plot(rhat2[(burntime+1):N], type="l", xlab="", ylab="R",main="the Gelman-Rubin statistic for Y")
    abline(h=1.2, lty=2,col="red")


## -----------------------------------------------------------------------------

#compute the Euclidean norm of a vector
norm_vec<-function(x,k)    sqrt(sum(x*x)^(2*k+2))
#compute the kth term of the series
kth<-function(vector,d,k){
  stopifnot(d>=1) 
  term<-(-1)^k*exp(lgamma((d+1)/2)-lgamma(k+d/2+1)+lgamma(k+3/2)-lgamma(k+1))*norm_vec(a,k)/(2^k)/(2*k+1)/(2*k+2)

return (term)
}


## -----------------------------------------------------------------------------

#compute the sum of the series
series<-function(vector,d,eps){
  n<-0                           #the nth term
  total<-0                       #to store the sum
  repeat{
    total<-total+kth(vector,d,n)
    n<-n+1
    s<-abs(kth(vector,d,n))      #when   
  
  if (s<eps){
    break 
  }
  }
  return(total)
}


## -----------------------------------------------------------------------------
eps<-.Machine$double.eps^0.25
a<-c(1,2)
d<-2

print(paste("The sum of the series is",round(series(a,d,eps),5),"when a=(1,2)T"))


## -----------------------------------------------------------------------------
K<-c(4:10)
N<-length(K)
root<-numeric(N)
eps<-.Machine$double.eps^0.25

for (n in 1:N) {
  k<-K[n]

  f<-function(a){
  z1<-sqrt(a^2*(k-1)/(k-a^2))
  z2<-sqrt(a^2*k/(k+1-a^2))
  f1<-function(x){
exp(lgamma(k/2)-lgamma((k-1)/2))*(1+(z1*x)^2 /(k-1))^(-k/2) *z1/sqrt((k-1)*pi)-exp(lgamma((k+1)/2)-lgamma(k/2))*(1+(z2*x)^2 /k)^(-(k+1)/2) *z2/sqrt(k*pi) }
  return(integrate(f1,eps,1)$value)
  }

 times<-0
 b0<-0
 b1<-sqrt(k)-eps
 r<-seq(b0, b1, length=3)
 y<-c(f(r[1]),f(r[2]),f(r[3]))
 if (y[1]*y[3] > 0)
stop("f does not have opposite sign at endpoints")
 while(times<1000 && abs(y[2])>eps) {
    times<-times + 1
if (y[1]*y[2]<0) {
  r[3]<-r[2]
  y[3]<-y[2]
        }
 else {
   r[1]<-r[2]
  y[1]<-y[2]
}
r[2]<-(r[1]+r[3])/2
y[2]<-f(r[2])

}

 root[n]<-r[2]
}

## -----------------------------------------------------------------------------
# results<-matrix(0,nrow=7,ncol=1)
# rownames(results)<-4:10
# colnames(results)<-c("Solution")
# results[,1]<-root
# 
# library(knitr)
# knitr::kable(results,align=rep('c', 5))


## -----------------------------------------------------------------------------
set.seed(2021)
K1<-c(4:25,100,500,1000)
N1<-length(K1)
eps<-.Machine$double.eps^0.25
root2<-numeric(N1)

#root-finding function,with regard to a
findinter<-function(a) {
  c1<-sqrt(a^2*(k-1)/(k-a^2))
  c2<-sqrt(a^2*k/(k+1-a^2))
  difference<-pt(c1,df=k-1)-pt(c2,df=k)
  #lower.tail=TRUE by default,computing P(X<=x)
  return (difference)
 }

for (n in 1:N1) {
  k<-K1[n]
  b0<-eps
  b1<-sqrt(k)-eps

 root2[n]<-uniroot(findinter,interval=c(b0, b1))$root
}

results2<-matrix(0,nrow=25,ncol=1)
rownames(results2)<-K1
colnames(results2)<-c("Intersection point")
results2[,1]<-root2

library(knitr)
knitr::kable(results2,align=rep('c', 5))

## -----------------------------------------------------------------------------
Y<-c(0.54,0.48,0.33,0.43,0.91,0.21,0.85,1,1,1)

#log-likelihood for observed data 
loglike<- function(lambda) {
        -7*log(lambda)-sum(Y)/lambda
}

#MLE for observed data 
lambdahat<-optimize(loglike,c(0,3),maximum = TRUE)$maximum

prior<-posterior<-lambdahat     #initialize lambda^(k)
eps<-1                

#log-likelihood for complete data 
conditional<-function(lambda) {
        -10*log(lambda)-(sum(Y)+3*prior)/lambda
}

#E-M iteration:
  while (eps>.Machine$double.eps^0.25) {
    prior<-posterior
    posterior<-optimize(conditional,c(0,2),maximum = TRUE)$maximum
    eps<-abs(posterior-prior)
  }


#Present the result
results3<-matrix(0,nrow=1,ncol=2)
rownames(results3)<-c("lambda")
colnames(results3)<-c("MLE for observed data","E-M result")
results3[,1]<-lambdahat
results3[,2]<-posterior

library(knitr)
knitr::kable(results3,align=rep('c', 5))


## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
#For conciseness,use unlist() to convert the output from a list to a vector.
#the trim option in mean() is the fraction (0 to 0.5) of observations to be trimmed from each end of x before the mean is computed. 
unlist(lapply(trims, function(trim) mean(x, trim = trim)))
unlist(lapply(trims, mean, x = x))


## -----------------------------------------------------------------------------
# Exercise 3
data("mtcars")
set.seed(2021)

#using the formulas stored in this list
formulas <- list( 
  mtcars$mpg ~ mtcars$disp,
  mtcars$mpg ~ I(1 / mtcars$disp),
  mtcars$mpg ~ mtcars$disp + mtcars$wt,
  mtcars$mpg ~ I(1 / mtcars$disp) + mtcars$wt
)
n<-length(formulas)

# loop
result<-list(n)
for (i in 1:n){
   result[[i]]<-summary(lm(formulas[[i]]))$r.squared
}
unlist(result)

# lapply
unlist(lapply(formulas,function(x){
  summary(lm(x))$r.squared
}))


## -----------------------------------------------------------------------------
# Exercise 4

bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

# loop
result2<-list(10)
for(i in 1:10){
  result2[[i]]<-summary(lm(bootstraps[[i]]$mpg~bootstraps[[i]]$disp))$r.squared
}
unlist(result2)

# lapply
unlist(lapply(bootstraps,function(x){
      summary(lm(x$mpg~x$disp))$r.squared}))

## -----------------------------------------------------------------------------

set.seed(2021)
x<-runif(1000,0,1)
y<-rexp(1000,1)
df<-data.frame(x,y)
vapply(df,sd,numeric(1))

## -----------------------------------------------------------------------------
df2<-data.frame(x=1:5,y=c("a","b","c","d","e"))
vapply(df2[vapply(df2,is.numeric,logical(1))],sd,numeric(1))

## -----------------------------------------------------------------------------
# Error: processing vignette 'homework.Rmd' failed with diagnostics:
#  4 simultaneous processes spawned,so don't present the result
# library(parallel)
# options(warn=-1)           #neglect all the warnings
# cores<-detectCores()
# cluster<-makePSOCKcluster(cores)
# 
# mcsapply<-function(x,f){
#        out<-parSapply(cluster,x,f)
#        simplify2array(out)
#    }
# 
# #Take an easy problem for example.
# system.time(mcsapply(x=1:1e6,f=mean))
 

## -----------------------------------------------------------------------------
# library(parallel)
# options(warn=-1)           #neglect all the warnings
# cores<-detectCores()
# cluster<-makePSOCKcluster(cores)
# 
# #FUN.VALUE=1.0 means the return value should be double(numeric) type.
# mcvapply<-function(x,f,FUN.VALUE) {
#        out_list<-parSapply(cluster,x,f)
#        out<-matrix(rep(FUN.VALUE,length(x)),nrow=length(x))
#        for (i in seq_along(x)){
#            out[i,]<-out_list[[i]]
#        }
#        out
# }
# 
# #Take an easy problem for example.
# system.time(mcvapply(x=1:1e6,f=mean,FUN.VALUE = 1.0))


## -----------------------------------------------------------------------------
options(warn=-1)
library(Rcpp)
library(microbenchmark)

#Gibbs Sampler method using R
GibbsR<-function (a,b,N) {
X<-matrix(0,N,2)            
X[1,]<-c(0,0)                 #initialize       
for (i in 2:N) {
x2<-X[i-1,2]
X[i,1]<-rbinom(1,size=n,prob=x2)
x1<-X[i,1]
X[i,2]<-rbeta(1,x1+a,n-x1+b)
}
return (X)
}

## -----------------------------------------------------------------------------
#Gibbs Sampler method using C

cppFunction(
  'NumericMatrix GibbsC (int a,int b,int N){
  NumericMatrix X(N,2); 
  X(0,0)=0;
  X(0,1)=0;
  double n=25;
  double x1=0;
  double x2=0;
  for (int i=1;i<N;i++){ 
    x2=X(i-1,1);
    X(i,0)=rbinom(1,n,x2)[0];
    x1=X(i,0);
    X(i,1)=rbeta(1,x1+a,n-x1+b)[0];
  }
  return X;
}')

## -----------------------------------------------------------------------------
a<-1       
b<-1
n<-25
N<-10000                 #length of chain

GR<-GibbsR(a,b,N)  
GC<-GibbsC(a,b,N)
#Binomial condition function
qqplot(GR[,1],GC[,1],xlab = "R function",ylab="C++ function",main= paste("a=",a,",b=",b,",n=",n,",Condition function:Binomial"))
#Beta condition function
qqplot(GR[,2],GC[,2],xlab = "R function",ylab="C++ function",main= paste("a=",a,",b=",b,",n=",n,",Condition function:Beta"))


## -----------------------------------------------------------------------------
ts<-microbenchmark(GibbsUsingR=GibbsR(a,b,N), GibbsUsingC=GibbsC(a,b,N))

library(knitr)
knitr::kable(summary(ts)[,c(1,3,4,5,6)],align=rep('c', 5))


