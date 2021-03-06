
gest.cat<-function(data,Yn,An,Ybin,Lny,Lnp,z,type=NA,timevarying=FALSE,Cn=NA,LnC=NA,...){
#########################Shared content################################
#Get the data, and create lagged and leading variables
T<-max(data$time)

if(is.na(type)==TRUE){
type<-5}

d<-data

d<-slide(data=d,Var=An,GroupVar="id",slideBy=-1,NewVar=paste("L",1,An,sep=""),reminder=F)


#Get lag and lead names and place them into vectors

Alagn<-paste("L",1,An,sep="") 

d[is.na(d[,Alagn]),][,Alagn]<-1

d[,Alagn]<-as.factor(d[,Alagn])
d$int<-1

#Propensity Scores
lmp<-reformulate(termlabels=c(Alagn,"as.factor(time)",Lnp),response=An)

modp<-multinom(lmp,family="binomial",data=d)


props<-predict(modp,type="probs",newdata=d)
d$prs<-props[,-1]





#Censoring Scores

if(is.na(Cn)==TRUE){
d$w<-1} else{

lmc<-reformulate(termlabels=c(LnC,"as.factor(time)"),response=Cn)

modc<-glm(lmc,family="binomial",data=d)



cps<-1-predict(modc,type="response",newdata=d)
d$cps<-cps

#Indicator Function of whether C=0
nm<-paste(Cn,"0",sep="")
d[,nm]<-as.integer(!d[,Cn])

#Setting up the denominator product of censoring weight "cprod" and weights "w"
d$cprod<-d$cps

##Set cprod to 1 where it is missing to allow for censoring weights
d[is.na(d$cprod)==TRUE,"cprod"]<-1


d$w<-d[,nm]/d$cps
}


dcom<-d[complete.cases(d),]

d$H<-0
d[d$time==T,"H"]<-d[d$time==T,Yn]
dcom$H<-0
dcom[dcom$time==T,"H"]<-dcom[dcom$time==T,Yn]

######################Perform Algorithm###########
#Get adjusted Outcome Model

#Create types

if (type==1){
z<-c(1)
timevarying<-FALSE
}else if (type==2){
z<-c(1,Lny[1])
timevarying<-FALSE
}else if (type==3){
z<-c(1)
timevarying<-TRUE
}else if (type==4){
z<-c(1,Lny[1])
timevarying<-TRUE
}

#Get correct family for outcome models
if(Ybin==TRUE){
family<-Gamma(link="log")
} else {
family<-gaussian
}


z[z==1]<-"int"
int1<-paste(eval(An),eval(z),sep=":")
int1c<-paste(eval(An),"int",sep=":")
int1[int1==int1c]<-paste(eval(An))
int2<-paste("prs",eval(z),sep=":")
int2c<-paste("prs","int",sep=":")
int2[int2==int2c]<-paste("prs")



if (timevarying==FALSE){
lmy<-reformulate(termlabels=c(int1,int2,Alagn,Lny),response=Yn)
lmH<-reformulate(termlabels=c(int1,int2,Alagn,Lny),response="H")
dt<-dcom[dcom$time==T,]
out<-summary(geem(terms(lmy),family=family,id=dt$id,data=dt,weights=dt$w))


nam1<-paste(An,levels(d[,An])[-1],sep="")

nam2<-apply(expand.grid(nam1,z[-1]), 1, paste, collapse=":")

Acoef<-c(nam1,nam2)

psi0<-out$beta[match(Acoef,out$coefnames)]
names(psi0)<-Acoef
psimar<-out$beta[match(z[-1],out$coefnames)]

psicat<-as.list(NULL)

for (l in 2:nlevels(d[,An])){
psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
}
psicat[[1]]<-rep(0,length(psicat[[2]]))

#Counterfactual Step
i<-(T-1)
while(i>=1){
###Obtain psiZA for appropriate times

j<-T
d$psiZA<-0

while(j>=(i+1)){
if(length(z)==1){
for (l in 1:nlevels(d[,An])){
d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),"psiZA"]<-psicat[[l]]
}
}else{
for (l in 1:nlevels(d[,An])){
d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),"psiZA"]<-rowSums(
sweep(d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),z],2,psicat[[l]],"*"))
}
}
j<-j-1
}


	k<-(T-1)
	while(k>=i){

	if(is.na(Cn)==FALSE){
	d[d$time==k,"cprod"]<-d[d$time==k+1,"cprod"]*d[d$time==k,"cps"]
	
	d[d$time==k,"w"]<-d[d$time==T,nm]/d[d$time==k,"cprod"]
	
		}
	if (Ybin==FALSE){	
	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	d[d$time==(k+1),"psiZA"]
	}
	if (Ybin==TRUE){
	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	exp(-d[d$time==(k+1),"psiZA"])
	}


	k<-k-1
	}



dt<-d[d$time %in% seq(i,T,by=1),]
dtcom<-dt[complete.cases(dt),]
out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))

psi0<-out$beta[match(Acoef,out$coefnames)]
	names(psi0)<-Acoef
	psimar<-out$beta[match(z[-1],out$coefnames)]
	psicat<-as.list(NULL)

	for (l in 2:nlevels(d[,An])){
	psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
	}
	psicat[[1]]<-rep(0,length(psicat[[2]]))

i<-i-1

}
results<-unlist(psicat[-1])

return(list(psi=results))
}
#############Time varying value##############
else if (timevarying==TRUE){
lmy<-reformulate(termlabels=c(int1,int2,Alagn,Lny),response=Yn)
lmH<-reformulate(termlabels=c(int1,int2,Alagn,Lny),response="H")
lmH1<-reformulate(termlabels=c(int1,int2,Lny),response="H")
dt<-dcom[dcom$time==T,]
out<-summary(geem(terms(lmy),family=family,id=dt$id,data=dt,weights=dt$w))


nam1<-paste(An,levels(d[,An])[-1],sep="")

nam2<-apply(expand.grid(nam1,z[-1]), 1, paste, collapse=":")

Acoef<-c(nam1,nam2)

psi0<-out$beta[match(Acoef,out$coefnames)]
names(psi0)<-Acoef


psicat<-as.list(NULL)

for (l in 2:nlevels(d[,An])){
psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
}
psicat[[1]]<-rep(0,length(psicat[[2]]))
psicatlist<-as.list(NULL)
psicatresult<-as.list(NULL)
psicatlist[[T]]<-psicat
psicatresult[[T]]<-psicat[-1]

#Counterfactual Step
i<-(T-1)
while(i>=1){
###Obtain psiZA for appropriate times

j<-T
d$psiZA<-0

while(j>=(i+1)){
if(length(z)==1){
for (l in 1:nlevels(d[,An])){
d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),"psiZA"]<-psicatlist[[j]][[l]]
}
}else{
for (l in 1:nlevels(d[,An])){
d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),"psiZA"]<-rowSums(
sweep(d[d$time==j & d[,An]==levels(d[,An])[l] & !is.na(d[,An]),z],2,psicatlist[[j]][[l]],"*"))
}
}
j<-j-1
}


	k<-(T-1)
	while(k>=i){

	if(is.na(Cn)==FALSE){
	d[d$time==k,"cprod"]<-d[d$time==k+1,"cprod"]*d[d$time==k,"cps"]
	
	d[d$time==k,"w"]<-d[d$time==T,nm]/d[d$time==k,"cprod"]
	
		}
	if (Ybin==FALSE){	
	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	d[d$time==(k+1),"psiZA"]
	}
	if (Ybin==TRUE){
	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	exp(-d[d$time==(k+1),"psiZA"])
	}


	k<-k-1
	}



dt<-d[d$time %in% i,]
dtcom<-dt[complete.cases(dt),]

if (i==1){
out<-summary(geem(terms(lmH1),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
}else{
out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
}

psi0<-out$beta[match(Acoef,out$coefnames)]
	names(psi0)<-Acoef
	psimar<-out$beta[match(z[-1],out$coefnames)]
	psicat<-as.list(NULL)

	for (l in 2:nlevels(d[,An])){
	psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
	}
	psicat[[1]]<-rep(0,length(psicat[[2]]))
	
psicatresult[[i]]<-psicat[-1]
psicatlist[[i]]<-psicat

i<-i-1

}

nam<-as.vector(NULL)
for (p in 1:T){
nam[p]<-paste("t=",p,sep="")
}
names(psicatresult)<-nam
results<-unlist(psicatresult)
return(list(psi=results))


}
}
