gest.boot<-function(data,gestfunc,Yn,An,Ybin,Abin,Lny,Lnp,
z,type=NA,timevarying=FALSE,Cn=NA,LnC=NA,cutoff=NA,bn,alpha,...){

T<-max(data$time)
t0<-gestfunc(data=data,Yn=Yn,An=An,Ybin=Ybin,Abin=Abin,Lny=Lny,Lnp=Lnp,z=z,type=type,
timevarying=timevarying,Cn=Cn,LnC=LnC,cutoff=cutoff)$psi
nams<-names(t0)
#Create tibble data based on ID
Data<-data %>% nest_legacy(-id)

bs <- bootstraps(Data, times = bn)
#bs

#start<-Sys.time()
#Get bootstrapped samples
results1<-as.list(NULL)
results<-as.list(NULL)
for (j in 1:bn){
tryCatch({
b<-as.data.frame(as_tibble(bs$splits[[j]]) %>% unnest_legacy())

#########Treat each set of T as a new unique ID
b$id<-sort(rep(1:length(unique(data$id)),T))

b<-b[order(b$id,b$time),]

results1[[j]]<-gestfunc(data=b,Yn=Yn,An=An,Ybin=Ybin,Abin=Abin,Lny=Lny,Lnp=Lnp,z=z,type=type,
timevarying=timevarying,Cn=Cn,LnC=LnC,cutoff=cutoff)$psi
results[[j]]<-unlist(results1[[j]])
},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

mean<-colMeans(do.call(rbind,results))

sd<-colStdevs(do.call(rbind,results))

bias<-mean-unlist(t0)


low<-unlist(t0)-(qnorm(1-(alpha/2)))*sd

upp<-unlist(t0)+(qnorm(1-(alpha/2)))*sd

conf<-cbind(low,upp)

lowb<-unlist(t0)-(qnorm(1-((alpha/2)/length(unlist(t0)))))*sd

uppb<-unlist(t0)+(qnorm(1-((alpha/2)/length(unlist(t0)))))*sd

confb<-cbind(lowb,uppb)
return(list(original=t0,conf=conf,conf.Bonferroni=confb,
mean=mean,s.e=sd))

}

