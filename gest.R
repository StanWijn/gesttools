
gest<-function(data,Yn,An,Ybin,Abin,Lny,Lnp,z,type=NA,timevarying=FALSE,Cn=NA,LnC=NA,...){
#########################Shared content################################
#Get the data, and create lagged and leading variables
T<-max(data$time)

if(is.na(type)==TRUE){
type<-5}

d<-data

d<-slide(data=d,Var=An,GroupVar="id",slideBy=-1,NewVar=paste("L",1,An,sep=""),reminder=F)


#Get lag and lead names and place them into vectors

Alagn<-paste("L",1,An,sep="") 

d[is.na(d[,Alagn]),][,Alagn]<-0

d$int<-1

#Propensity Scores
lmp<-reformulate(termlabels=c(Alagn,"as.factor(time)",Lnp),response=An)

if (Abin==TRUE){
modp<-glm(lmp,family="binomial",data=d)
} else {
modp<-glm(lmp,family="gaussian",data=d)
}

props<-predict(modp,type="response",newdata=d)
d$prs<-props



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
psi<-out$beta[match(int1,out$coefnames)]
names(psi)<-int1
#se.robust<-out$se.robust[match(int1,out$coefnames)]

#Counterfactual Step
i<-(T-1)
while(i>=1){

	k<-(T-1)
	while(k>=i){

	if(is.na(Cn)==FALSE){
	d[d$time==k,"cprod"]<-d[d$time==k+1,"cprod"]*d[d$time==k,"cps"]
	
	d[d$time==k,"w"]<-d[d$time==T,nm]/d[d$time==k,"cprod"]
	
		}
	if (Ybin==FALSE){
	if (length(z)==1){
	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	d[d$time==(k+1),An]*psi*d[d$time==(k+1),z]

	}else{

	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	rowSums(d[d$time==(k+1),An]*sweep(d[d$time==(k+1),z],2,psi,"*"))
	}
	}

	if (Ybin==TRUE){
	if (length(z)==1){
	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	exp(-d[d$time==(k+1),An]*psi*d[d$time==(k+1),z])

	}else{

	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	exp(-rowSums(d[d$time==(k+1),An]*sweep(d[d$time==(k+1),z],2,psi,"*")))
	}
	}


	k<-k-1
	}



dt<-d[d$time %in% seq(i,T,by=1),]
dtcom<-dt[complete.cases(dt),]
out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
psi<-out$beta[match(int1,out$coefnames)]
names(psi)<-int1
#se.robust<-out$se.robust[match(int1,out$coefnames)]
i<-i-1

}

return(list(psi=psi))
}
#############Time varying value##############
else if (timevarying==TRUE){
lmy<-reformulate(termlabels=c(int1,int2,Alagn,Lny),response=Yn)
lmH<-reformulate(termlabels=c(int1,int2,Alagn,Lny),response="H")
lmH1<-reformulate(termlabels=c(int1,int2,Lny),response="H")

dt<-dcom[dcom$time==T,]

out<-summary(geem(terms(lmy),family=family,id=dt$id,data=dt,weights=dt$w))

psi<-out$beta[match(int1,out$coefnames)]
names(psi)<-int1
se.robust<-out$se.robust[match(int1,out$coefnames)]

psilist<-as.list(NULL)
se.robust<-as.list(NULL)
psilist[[T]]<-psi
se.robust[[T]]<-out$se.robust[match(int1,out$coefnames)]

i<-(T-1)
while(i>=1){

	k<-(T-1)
	while(k>=i){

	if(is.na(Cn)==FALSE){
	d[d$time==k,"cprod"]<-d[d$time==k+1,"cprod"]*d[d$time==k,"cps"]
	
	d[d$time==k,"w"]<-d[d$time==T,nm]/d[d$time==k,"cprod"]
	
		}
	if (Ybin==FALSE){
	if (length(z)==1){
	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	d[d$time==(k+1),An]*psilist[[k+1]]*d[d$time==(k+1),z]

	}else{

	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]-
	rowSums(d[d$time==(k+1),An]*sweep(d[d$time==(k+1),z],2,psilist[[k+1]],"*"))
	}
	}

	if (Ybin==TRUE){
	if (length(z)==1){
	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	exp(-d[d$time==(k+1),An]*psilist[[k+1]]*d[d$time==(k+1),z])

	}else{

	d[d$time==k,"H"]<-d[d$time==(k+1),"H"]*
	exp(-rowSums(d[d$time==(k+1),An]*sweep(d[d$time==(k+1),z],2,psilist[[k+1]],"*")))
	}
	}


	k<-k-1
	}



	dt<-d[d$time %in% i,]
	dtcom<-dt[complete.cases(dt),]
	if (i==1){
	out<-summary(geem(terms(lmH1),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	psi<-out$beta[match(int1,out$coefnames)]
	names(psi)<-int1
	psilist[[i]]<-psi
	se.robust[[i]]<-out$se.robust[match(int1,out$coefnames)]

	}else{
	out<-summary(geem(terms(lmH),family=family,id=dtcom$id,data=dtcom,weights=dtcom$w))
	psi<-out$beta[match(int1,out$coefnames)]
	names(psi)<-int1
	psilist[[i]]<-psi
	se.robust[[i]]<-out$se.robust[match(int1,out$coefnames)]

	}
	i<-i-1
	}
nam<-as.vector(NULL)
for (p in 1:T){
nam[p]<-paste("t=",p,sep="")
}
names(psilist)<-nam
results<-unlist(psilist)
return(list(psi=results))
}

}


