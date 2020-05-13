gestmult.cat<-function(data,Yn,An,Ybin,Lny,Lnp,z,type=NA,timevarying=FALSE,Cn=NA,LnC=NA,cutoff=NA,...){
#Get the data, and create lagged and leading variables

T<-max(data$time)

if(is.na(cutoff)==TRUE){
cutoff<-T}

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
d$w<-1
} else{

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

##Set up 1 step counterfactuals
d$H<-d[,Yn]
dc<-d
dc$cntstep<-1

###Set up required augmented data subset dc, with an indicator of the step length for the 
###counterfactual



for (i in 2:cutoff){
d2<-d[d$time %in% seq(1,T-(i-1),by=1),]
d2$cntstep<-i
dc<-rbind(dc,d2)
}
dc<-dc[order(dc$id,dc$time),]

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

#Collect correct values for psi
out<-summary(geem(terms(lmy),family=family,id=dcom$id,data=dcom,weights=dcom$w))

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


i<-2
	while(i<=cutoff && i<=T){
	
#Need to obtain psi*Z*A for all relevant data points before calculating H
#call this psiZA
#i=Maximum Counterfactual step
#j=Current Counterfactual step
#l<-category level
j<-1
dc$psiZA<-0
while(j<=(i-1)) {

if (length(z)==1){
	for (l in 1:nlevels(dc[,An])){
	dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),"psiZA"]<-psicat[[l]]}
	}else{
	for (l in 1:nlevels(dc[,An])){
	dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),"psiZA"]<-rowSums(
	sweep(dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),z],2,psicat[[l]],"*"))
	}
	}


j<-j+1
}

#Now calculate the counterfactuals and censoring weights
	j<-2

	while(j<=i) {
	for(k in 1:(T-(j-1))){ 
	#Censoring Weights
	
	if(is.na(Cn)==FALSE){

	dc[dc$cntstep==j & dc$time==k,"cprod"]<-dc[dc$cntstep==(j-1) & dc$time==k,"cps"]*
	dc[dc$cntstep==(j-1) & dc$time==(k+1),"cprod"]
	
	dc[dc$cntstep==j & dc$time==k,"w"]<-d[d$time==(k+j-1),nm]/dc[dc$cntstep==j & dc$time==k,"cprod"]
	}

	#Counterfactuals
	if (Ybin==FALSE){
	dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]-
	dc[dc$cntstep==(j-1) & dc$time==(k+1),"psiZA"]
	}
	if (Ybin==TRUE){
	dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]*
	exp(-dc[dc$cntstep==(j-1) & dc$time==(k+1),"psiZA"])
	}
	
	}

	j<-j+1
	}
	#Update estimates
	dt<-dc[dc$cntstep %in% seq(1,i,by=1),]
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
	i<-i+1
	}
#Alternate parametrisation
psicatsum<-as.list(NULL)
for (m in 1:length(psicat)){
a<-psicat[[m]][1]
b<-sum(psicat[[m]])
psicatsum[[m]]<-c(a,b)
}

results<-unlist(psicat[-1])
return(list(psi=results))
	}
else if (timevarying==TRUE){
lmy<-reformulate(termlabels=c(int1,int2,Alagn,Lny),response=Yn)
lmH<-reformulate(termlabels=c(int1,int2,Alagn,Lny),response="H")
lmH1<-reformulate(termlabels=c(int1,int2,Lny),response="H")

out<-summary(geem(terms(lmy),family=family,id=dcom$id,data=dcom,weights=dcom$w))

nam1<-paste(An,levels(d[,An])[-1],sep="")

nam2<-apply(expand.grid(nam1,z[-1]), 1, paste, collapse=":")

Acoef<-c(nam1,nam2)

psi0<-out$beta[match(Acoef,out$coefnames)]
names(psi0)<-Acoef
psimar<-out$beta[match(z[-1],out$coefnames)]

#Overall time varying list
psicatlist<-as.list(NULL)
psicatresult<-as.list(NULL)
psicat<-as.list(NULL)

#First time varying values
for (l in 2:nlevels(d[,An])){
psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
}
psicat[[1]]<-rep(0,length(psicat[[2]]))

psicatlist[[1]]<-psicat
psicatresult[[1]]<-psicat[-1]

i<-2
	while(i<=cutoff && i<=T){
	
#Need to obtain psi*Z*A for all relevant data points before calculating H
#call this psiZA
#i=Maximum Counterfactual step
#j=Current Counterfactual step
#l<-category level
j<-1
dc$psiZA<-0
while(j<=(i-1)) {

if (length(z)==1){
	for (l in 1:nlevels(dc[,An])){
	dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),"psiZA"]<-psicatlist[[j]][[l]]}
	}else{
	for (l in 1:nlevels(dc[,An])){
	dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),"psiZA"]<-rowSums(
	sweep(dc[dc$cntstep==j & dc[,An]==levels(dc[,An])[l] & !is.na(dc[,An]),z],2,psicatlist[[j]][[l]],"*"))
	}
	}


j<-j+1
}

#Now calculate the counterfactuals and censoring weights
	j<-2

	while(j<=i) {
	for(k in 1:(T-(j-1))){ 
	#Censoring Weights
	
	if(is.na(Cn)==FALSE){

	dc[dc$cntstep==j & dc$time==k,"cprod"]<-dc[dc$cntstep==(j-1) & dc$time==k,"cps"]*
	dc[dc$cntstep==(j-1) & dc$time==(k+1),"cprod"]
	
	dc[dc$cntstep==j & dc$time==k,"w"]<-d[d$time==(k+j-1),nm]/dc[dc$cntstep==j & dc$time==k,"cprod"]
	}

	#Counterfactuals
	if (Ybin==FALSE){
	dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]-
	dc[dc$cntstep==(j-1) & dc$time==(k+1),"psiZA"]
	}
	if (Ybin==TRUE){
	dc[dc$cntstep==j & dc$time==k,"H"]<-dc[dc$cntstep==(j-1) & dc$time==(k+1),"H"]*
	exp(-dc[dc$cntstep==(j-1) & dc$time==(k+1),"psiZA"])
	}
	
	}

	j<-j+1
	}
	#Update estimates
	dt<-dc[dc$cntstep %in% i,]
	dtcom<-dt[complete.cases(dt),]

	if (i==T){
	out<-summary(geem(terms(lmH1),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	
	psi0<-out$beta[match(Acoef,out$coefnames)]
	psimar<-out$beta[match(z[-1],out$coefnames)]

	}else{
	out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	
	psi0<-out$beta[match(Acoef,out$coefnames)]
	psimar<-out$beta[match(z[-1],out$coefnames)]
	}
	names(psi0)<-Acoef
	psicat<-as.list(NULL)

	for (l in 2:nlevels(d[,An])){
	psicat[[l]]<-psi0[grep(levels(d[,An])[l],Acoef)]
	}
	psicat[[1]]<-rep(0,length(psicat[[2]]))
	psicatlist[[i]]<-psicat
	psicatresult[[i]]<-psicat[-1]
	i<-i+1
	}
#Alternate parametrisation
psicatsum<-as.list(NULL)
for (m in 1:length(psicat)){
a<-psicat[[m]][1]
b<-sum(psicat[[m]])
psicatsum[[m]]<-c(a,b)
}

nam<-as.vector(NULL)
for (p in 1:T){
nam[p]<-paste("s-",p,sep="")
}
names(psicatresult)<-nam

results<-unlist(psicatresult)


return(list(psi=results))

	}
}
