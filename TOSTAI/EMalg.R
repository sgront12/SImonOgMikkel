#EM-algoritmen
#Funktion, som estimere parametrene i en bivariat normalfordeling vha. EM-algoritmen
bivnEM=function(x, tol){
  #x er en nx2 matrix bestaaende af realiseringer af en bivariat normalfordeling. tol er tolerancen 
  
  #Antal realiseringer
  n=dim(x)[1]
  #Identifikation af observeret samt manglende data
  na.ind=!is.na(x)

  #Oversigt
  x1.obs=na.ind[,1]
  x2.obs=na.ind[,2]
  both.obs=x1.obs & x2.obs
  x1.obs.x2.mis=x1.obs & !x2.obs
  x1.mis.x2.obs=x2.obs & !x1.obs
  both.mis=!x1.obs & !x2.obs
  
  #Vektor til theta_t1
  theta_t1=numeric(5)

  #Begyndelsesbetingelserne er bestemt ud fra de realiseringer, hvor baade x1 og x2 er kendt
  theta_t1[1]=sum(x[both.obs,1])/length(x[both.obs,1])
  theta_t1[2]=sum(x[both.obs,2])/length(x[both.obs,2])
  theta_t1[3]=sum(x[both.obs,1]^2)/length(x[both.obs,1])-theta_t1[1]^2
  theta_t1[4]=sum(x[both.obs,2]^2)/length(x[both.obs,2])-theta_t1[2]^2
  theta_t1[5]=sum(x[both.obs,1]*x[both.obs,2])/length(x[both.obs,2])-theta_t1[1]*theta_t1[2]
  
  #Vektor til theta_t
  theta_t=numeric(5)

  #Vektor til den sufficiente stikproevefunktion
  s_t=numeric(5)

  #Vektor til log-likelihood
  ll=c()
  
  #Variabel til bestemmelse af antal iterationer
  k=0
  
  #While-loekke, som udfoerer EM-algoritmen
  while(all((abs(theta_t1-theta_t)) > tol)){
    theta_t=theta_t1
    #E-trin
    s_t[1]=sum(x[x1.obs,1])+sum(theta_t[1]+theta_t[5]/theta_t[4]*(x[x1.mis.x2.obs,2]-theta_t[2]))+sum(theta_t[1]*both.mis)
    s_t[2]=sum(x[x2.obs,2])+sum(theta_t[2]+theta_t[5]/theta_t[3]*(x[x1.obs.x2.mis,1]-theta_t[1]))+sum(theta_t[2]*both.mis)
    s_t[3]=sum(x[x1.obs,1]^2)+sum(theta_t[3]-theta_t[5]^2/theta_t[4]+(theta_t[1]+theta_t[5]/theta_t[4]*(x[x1.mis.x2.obs,2]-theta_t[2]))^2)+sum((theta_t[3]+theta_t[1]^2)*both.mis)
    s_t[4]=sum(x[x2.obs,2]^2)+sum(theta_t[4]-theta_t[5]^2/theta_t[3]+(theta_t[2]+theta_t[5]/theta_t[3]*(x[x1.obs.x2.mis,1]-theta_t[1]))^2)+sum((theta_t[4]+theta_t[2]^2)*both.mis);
    s_t[5]=sum(x[both.obs,1]*x[both.obs,2])+sum((theta_t[1]+theta_t[5]/theta_t[4]*(x[x1.mis.x2.obs,2]-theta_t[2]))*x[x1.mis.x2.obs,2])+sum(x[x1.obs.x2.mis,1]*(theta_t[2]+theta_t[5]/theta_t[3]*(x[x1.obs.x2.mis,1]-theta_t[1])))+sum((theta_t[5]+theta_t[1]*theta_t[2])*both.mis)
    
    #Log-likelihood
    ll[k+1] = -(length(which(both.obs))/2)*log(theta_t[3]*theta_t[4]-theta_t[5]^2)-1/(2*(theta_t[3]*theta_t[4]-theta_t[5]^2))*sum((x[both.obs,1]-theta_t[1])^2*theta_t[4]+(x[both.obs,2]-theta_t[2])^2*theta_t[3]-2*theta_t[5]*(x[both.obs,1]-theta_t[1])*(x[both.obs,2]-theta_t[2]))-
      (length(which(x1.obs.x2.mis))/2)*log(theta_t[3])-1/(2*theta_t[3])*sum((x[x1.obs.x2.mis,1]-theta_t[1])^2)-
      (length(which(x1.mis.x2.obs))/2)*log(theta_t[4])-1/(2*theta_t[4])*sum((x[x1.mis.x2.obs,2]-theta_t[2])^2)
    
    #M-trin
    theta_t1[1]=s_t[1]/n
    theta_t1[2]=s_t[2]/n
    theta_t1[3]=s_t[3]/n-theta_t1[1]^2
    theta_t1[4]=s_t[4]/n-theta_t1[2]^2
    theta_t1[5]=s_t[5]/n-theta_t1[1]*theta_t1[2]
  
    k=k+1
  }
  
  #Den sidste log-likelihood
  ll[k+1] = (-(length(which(both.obs))/2)*log(theta_t1[3]*theta_t1[4]-theta_t1[5]^2)-1/(2*(theta_t1[3]*theta_t1[4]-theta_t1[5]^2))*sum((x[both.obs,1]-theta_t1[1])^2*theta_t1[4]+(x[both.obs,2]-theta_t1[2])^2*theta_t1[3]-2*theta_t1[5]*(x[both.obs,1]-theta_t1[1])*(x[both.obs,2]-theta_t1[2]))
             -(length(which(x1.obs.x2.mis))/2)*log(theta_t1[3])-1/(2*theta_t1[3])*sum((x[x1.obs.x2.mis,1]-theta_t1[1])^2)
             -(length(which(x1.mis.x2.obs))/2)*log(theta_t1[4])-1/(2*theta_t1[4])*sum((x[x1.mis.x2.obs,2]-theta_t1[2])^2))

  #De estimerede parametre returneres
  return(list(mean=c(theta_t1[1],theta_t1[2]), cov=matrix(c(theta_t1[3], theta_t1[5], theta_t1[5], theta_t1[4]), nrow=2), LogLike=ll, iterationer=k))	
}

#Anvendelse af ovenstaaende funktion
load("~/Desktop/misdat.RData")
tol=1.e-10 #Tolerance
X=misdat #Datasaet

EM1=bivnEM(X,tol)
EM1$mean; EM1$cov; EM1$iterationer
plot(EM1$LogLike[-1],ylim=c(-140,-133))

#Datasaet uden realiseringer, hvor baade x1 og x2 mangler
X1=matrix(c(X[1:150,1],X[1:150,2]),ncol=2)

EM2=bivnEM(X1,tol)
EM2$mean; EM2$cov; EM2$iterationer
plot(EM2$LogLike[-1],ylim=c(-140,-133))
