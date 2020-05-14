 #Define the inverse logit function
 expit<-function(x){
 exp(x)/(1+exp(x))
 }
#####################Simulated Data With Binary Exposure###
 n<-10000
 set.seed(24092019)
 id<-seq(1,n,by=1)
 U<-rnorm(n,0,1)
 L.1<-rnorm(n,1+U,1)
 a.1<-expit(1+0.1*L.1)
 A.1<-rbinom(n,1,a.1)

 L.2<-rnorm(n,1+(A.1/3)+L.1+U,1)
 a.2<-expit(1+0.1*L.2+0.1*A.1)
 A.2<-rbinom(n,1,a.2)


 L.3<-rnorm(n,1+(A.2/2)+L.2+U,1)
 a.3<-expit(1+0.1*L.3+0.1*A.2)
 A.3<-rbinom(n,1,a.3)

 Y<-rnorm(n,1+(A.1/3)+(A.2/2)+A.3+L.1+L.2+L.3+U,1)

 #Set data as long format
 dw<-as.data.frame(cbind(id,Y,A.1,A.2,A.3,L.1,L.2,L.3,U))
 dl<-reshape(dw,direction="long",varying=c("A.1","A.2","A.3",
 "L.1","L.2","L.3"))
 #Order data by time and ID
 dl<-dl[order(dl$id,dl$time),]
 data<-dl

 #SNMM type 1
 gest(data=data,Yn="Y",An="A",Ybin=FALSE,Abin=TRUE,
 Lny=c("L","U"),Lnp=("L"),z=c(1),type=1,
 timevarying=FALSE,Cn=NA,LnC=NA)

 
 #SNMM type 2
 gest(data=data,Yn="Y",An="A",Ybin=FALSE,Abin=TRUE,
 Lny=c("L","U"),Lnp=("L"),z=c(1,"L"),type=2,
 timevarying=FALSE,Cn=NA,LnC=NA)

 
 #SNMM type 3
 gest(data=data,Yn="Y",An="A",Ybin=FALSE,Abin=TRUE,
 Lny=c("L","U"),Lnp=("L"),z=c(1),type=3,
 timevarying=TRUE,Cn=NA,LnC=NA)

 #SNMM type 3 bootstrap
 gest.boot(data=data,gestfunc=gest,Yn="Y",An="A",Ybin=FALSE,Abin=TRUE,
 Lny=c("L","U"),Lnp=c("L"),z=c(1),type=3,
 timevarying=TRUE,Cn=NA,LnC=NA,cutoff=NA,bn=5,alpha=0.05)
 
 #SNMM type 3 gestmult
 gestmult(data=data,Yn="Y",An="A",Ybin=FALSE,Abin=TRUE,
 Lny=c("L","U"),Lnp=("L"),z=c(1),type=3,
 timevarying=TRUE,Cn=NA,LnC=NA,cutoff=NA)


################Simulated Data With Categorical Exposure###############
 n<-10000
 set.seed(24092019)
 id<-seq(1,n,by=1)
 U<-rnorm(n,0,1)
 L.1<-rnorm(n,1+U,1)
 a.1<-expit(1+0.1*L.1)

 A.1<-as.vector(NULL)
 for (i in 1:n){
 A.1[i]<-sample(letters[1:3],1, replace=TRUE, prob=c(1-(3*a.1[i])/5,a.1[i]/5,2*(a.1[i]/5)))
 }
 A.1<-as.factor(A.1)
 A.1par<-as.vector(NULL)
 A.1par[A.1==letters[1]]<-0
 A.1par[A.1==letters[2]]<-1
 A.1par[A.1==letters[3]]<-2


 L.2<-rnorm(n,1+A.1par+L.1+U,1)
 a.2<-expit(1+0.1*L.2+A.1par)
 A.2<-as.vector(NULL)
 for (i in 1:n){
 A.2[i]<-sample(letters[1:3],1, replace=TRUE, prob=c(1-(3*a.2[i]/5),a.2[i]/5,2*(a.2[i]/5)))
 }
 A.2<-as.factor(A.2)
 A.2par<-as.vector(NULL)
 A.2par[A.2==letters[1]]<-0
 A.2par[A.2==letters[2]]<-1
 A.2par[A.2==letters[3]]<-2



 L.3<-rnorm(n,1+A.2par+L.2+U,1)
 a.3<-expit(1+0.1*L.3+A.2par)
 
 A.3<-as.vector(NULL)
 for (i in 1:n){
 A.3[i]<-sample(letters[1:3],1, replace=TRUE, prob=c(1-(3*a.3[i]/5),a.3[i]/5,2*(a.3[i]/5)))
 }
 A.3<-as.factor(A.3)
 A.3par<-as.vector(NULL)
 A.3par[A.3==letters[1]]<-0
 A.3par[A.3==letters[2]]<-1
 A.3par[A.3==letters[3]]<-2

 Y<-rnorm(n,1+A.1par+A.2par+A.3par+L.1+L.2+L.3+U,1)

 dw<-data.frame(id,Y,A.1,A.2,A.3,L.1,L.2,L.3,U)
 dl<-reshape(dw,direction="long",varying=c("A.1","A.2","A.3",
 "L.1","L.2","L.3"))
 #Order data by time and ID
 dl<-dl[order(dl$id,dl$time),]
 data<-dl
 dat$A<-as.factor(dat$A)

 #SNMM type 1
 gest.cat(data=data,Yn="Y",An="A",Ybin=FALSE,
 Lny=c("L","U"),Lnp=("L"),z=c(1),type=1,
 timevarying=FALSE,Cn=NA,LnC=NA)

 #SNMM type 3
 gest.cat(data=data,Yn="Y",An="A",Ybin=FALSE,
 Lny=c("L","U"),Lnp=("L"),z=c(1),type=3,
 timevarying=TRUE,Cn=NA,LnC=NA)

 #SNMM type 3 bootstrap
 gest.boot(data=data,gestfunc=gest.cat,Yn="Y",An="A",Ybin=FALSE,
 Lny=c("L","U"),Lnp=c("L"),z=c(1),type=3,
 timevarying=TRUE,Cn=NA,LnC=NA,cutoff=NA,bn=100,alpha=0.05)
 
 #SNMM type 3 gestmult
 gestmult.cat(data=data,Yn="Y",An="A",Ybin=FALSE,
 Lny=c("L","U"),Lnp=("L"),z=c(1),type=3,
 timevarying=TRUE,Cn=NA,LnC=NA,cutoff=NA)

