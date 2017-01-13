mnEM=function(x, tol){
  #x er en nx2 matrix bestaaende af realiseringer af en multinomiel fordeling.
  #tol er tolerancen. k er antallet af iterationer.
  k=0;
  #Antal realiseringer
  n=dim(x)[1]
  #Opstil tabel
  counts = table(x, useNA="always")
  
  #E og M vektorer
  #Vektor til pi_ij => (pi_00, pi_01, pi_10, pi_11)
  pi_t1=c(1/4,1/4,1/4,1/4)
  pi_t=numeric(4)
  
  #Vektor til E-trin
  s_t=numeric(4)
  #Likelihood-funktionsværdi-vektoren initialiseres
  like = c(0)
  
  #EM-Algoritmen
  while(all((abs(pi_t1-pi_t)) > tol) || k>1000){
    pi_t = pi_t1
    ##E-trin
    #SSH for bl.a. n_+0^J er pi_00+pi_10.
    #F.eks. Dermed er N_00^J = N_+0^J*pi_00/(pi_00+pi_10).
    s_t[1] = counts[1,1]+counts[1,3]*pi_t[1]/(pi_t[1]+pi_t[2])+counts[3,1]*pi_t[1]/(pi_t[1]+pi_t[3])
    s_t[2] = counts[1,2]+counts[1,3]*pi_t[2]/(pi_t[2]+pi_t[1])+counts[3,2]*pi_t[2]/(pi_t[2]+pi_t[4])
    s_t[3] = counts[2,1]+counts[2,3]*pi_t[3]/(pi_t[3]+pi_t[4])+counts[3,1]*pi_t[3]/(pi_t[3]+pi_t[1])
    s_t[4] = counts[2,2]+counts[2,3]*pi_t[4]/(pi_t[4]+pi_t[3])+counts[3,2]*pi_t[4]/(pi_t[4]+pi_t[2])

    #M-trin
    pi_t1[1] = s_t[1]/sum(counts)
    pi_t1[2] = s_t[2]/sum(counts)
    pi_t1[3] = s_t[3]/sum(counts)
    pi_t1[4] = s_t[4]/sum(counts)
    k <- k+1;
    
    #Indsæt ny information i likelihood-funktionen.
    like[k] = s_t[1]*log(pi_t1[1])+s_t[2]*log(pi_t1[2])+s_t[3]*log(pi_t1[3])+s_t[4]*log(pi_t1[4])
    
  }
  return(list(tal = counts, Estep = s_t, itt = k, pi = pi_t1, like=like));
}

#Anvendelse af ovenstaaende funktion
library(readr)
IJdata <- read_csv("IJdata.csv", col_types = cols())
View(IJdata)
tol=1.e-10 #Tolerance
X=IJdata #Datasaet

EM1=mnEM(X,tol) #Kører EM-algortimen
sum(EM1$Estep) #Ser, at data summer til 200.
sum(EM1$pi) #Ser, at parametre summer til 1.
EM1$pi #Estimerede parametre 
EM1$itt #Antallet af iterationer
plot(EM1$like) #Plot af likelihood-værdierne

#Sammenlign komplette data-estimater med ukomplette data.
tabel = table(X, useNA="no")
pi_right = c(0.2,0.075,0.3,0.425)
pi_comp = c(tabel[1,1]/sum(tabel),tabel[1,2]/sum(tabel),tabel[2,1]/sum(tabel),tabel[2,2]/sum(tabel))

#Vi finder differencen mellem genereringssandsynlighederne og estimaterne og sammenligner.
pi_difEM = abs(pi_right-EM1$pi)
pi_difcomp = abs(pi_right-pi_comp)
pi_difEM
pi_difcomp
#Der er umiddelbart ingen klar fordel i at vælge EM-algoritmen frem for MLE i dette tilfædle.