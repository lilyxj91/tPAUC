
#' Title
#'
#' @param beta 
#' @param tt 
#' @param nf 
#'
#' @return baseline PAUC
#' @export
#'
#' @examples
auc.pred.baseline <- function(beta,tt,nf=7){ 
  n0 <- length(beta)-nf-1
  res <- rep(0,length(tt))
  for (i in 1:length(tt)){
    t <- tt[i]
    polytt <- c(1,t^(-2),t^(-1),t^(-.5),log(t),t^(.5),t,t^2,rep(0,n0))
    res[i] <- exp(c(polytt %*% beta)) / (1 + exp(c(polytt %*% beta)) )
  }
  return(res)
}


#' Title
#'
#' @param n 
#' @param vt 
#' @param vc 
#' @param vm 
#' @param vxs 
#' @param u 
#'
#' @return reshaped data for easy optimization
#' @export
#'
#' @examples
data_crossingdc <- function(n,vt,vc,vm,vxs,u)
{
  #Y0,C0,M0,VXS
  #vt: vector of time: Y0
  #vc: vector of censoring: C0
  #vm: vector of biomarker: M0
  #vxs: matrix of xs: VXS
  nx<-ncol(vxs)
  nj<-length(vc)
  index<-which(vc==1)
  dat<-cbind(j=1:nj,vt,vm,vxs)
  dati<-dat[index,]%>%rename(k=j)
  nx<-ncol(vxs)
  names(dati)[1:3]<-c("k","yk","mk")
  names(dati)[4:(nx+3)]<-paste0(names(vxs),"k")
  names(dat)[1:3]<-c("j","yj","mj")
  names(dat)[4:(nx+3)]<-paste0(names(vxs),"j")
  h = 1*(u*n)^(-1/3)
  tableij<-crossing(dati,dat,.name_repair = "universal")
  tableij<-tableij%>%filter(yk<yj)%>%                               #filter ti<tj
    mutate(Ikj=as.numeric(mk>mj))%>%
    filter(xdk==xdj)%>%
    group_by(k)%>%arrange(k,mj)%>%mutate(dif=abs(xck-xcj),
                                         kh=0.75*(1-(dif/h)^2)/h*as.numeric(abs(dif/h)<1),
                                         rho01 = Ikj*kh,
                                         rho02 =(1-Ikj)*kh,
                                         rhoweight = rho01 + rho02,
    )%>%select(-c("mk","mj","dif","kh","yj","xcj"))%>%
    mutate(cumid=cumsum(rhoweight),
           nj=sum(rhoweight),
           #spe=cumid/nj,
           fp=1-cumid/nj)%>%
    select(-c("k","cumid","nj"))
  
  return(data.frame(tableij))
}

