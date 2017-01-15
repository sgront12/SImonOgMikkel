#Bayesianske netvaerk
library(gRim)
library(Rgraphviz)
X=read.csv("math.csv") #Datasaet

#Hver variabel i data inddeles i 3 niveauer
V=lapply(X,cut,breaks=3,labels=c("L","M","H"))
Y=matrix(0,nrow(X),ncol(X))
for(i in 1:ncol(X)){
  Y[,i]=V[[i]]
}
colnames(Y)=names(X)
Y[Y==1]="L"
Y[Y==2]="M"
Y[Y==3]="H"

#Frekvenstabel over data
ftab=as.table(ftable(xtabs(~me+ve+al+an+st,data=Y)))

#Test af betinget uafhaengighed
ciTest(ftab,set=c("an","st","al","me","ve"))
ciTest(ftab,set=c("an","ve","al","me"))
ciTest(ftab,set=c("ve","me","al"))
ciTest(ftab,set=c("st","me","al","ve"))
ciTest(ftab,set=c("an","me","al"))
ciTest(ftab,set=c("st","ve","al"))

#DAG
DAG=dag(~me|al + ve|al + an|al + st|al)
plot(DAG)

#Moraliseret graf
MG=ug(~me|al + ve|al + an|al + st|al)
plot(MG)

#Nummereret graf
NG=ug(~2|1 + 3|1 + 4|1 + 5|1)
plot(NG)

#Fordelingerne
al=armarg(ftab,~al); al
me_al=armarg(ftab,~me + al); me_al
ve_al=armarg(ftab,~ve + al); ve_al
an_al=armarg(ftab,~an + al); an_al
st_al=armarg(ftab,~st + al); st_al

#Normaliserede fordelinger
p.al=arnormalize(al, "first"); p.al
p.me_al=arnormalize(me_al, "first"); p.me_al
p.ve_al=arnormalize(ve_al, "first"); p.ve_al
p.an_al=arnormalize(an_al, "first"); p.an_al
p.st_al=arnormalize(st_al, "first"); p.st_al

#P(V) inddeles i 4 funktioner
q1.me.al=arprod(p.al,p.me_al); q1.me.al
q2.ve.al=p.ve_al; q2.ve.al
q3.an.al=p.an_al; q3.an.al
q4.st.al=p.st_al; q4.st.al

#Collect evidence
q4.al=armarg(q4.st.al,"al"); q4.al
q4.st.al=ardiv(q4.st.al,q4.al); q4.st.al
q3.an.al=armult(q3.an.al, q4.al); q3.an.al

q3.al=armarg(q3.an.al,"al"); q3.al
q3.an.al=ardiv(q3.an.al,q3.al); q3.an.al
q2.ve.al=armult(q2.ve.al, q3.al); q2.ve.al

q2.al=armarg(q2.ve.al,"al"); q2.al
q2.ve.al=ardiv(q2.ve.al,q2.al); q2.ve.al
q1.me.al=armult(q1.me.al, q2.al); q1.me.al

#Distribute evidence
q1.al=armarg(q1.me.al,"al"); q1.al
q2.ve.al=armult(q2.ve.al,q1.al); q2.ve.al

q2.al=armarg(q2.ve.al,"al"); q2.al
q3.an.al=armult(q3.an.al,q2.al); q3.an.al

q3.al=armarg(q3.an.al,"al"); q3.al
q4.st.al=armult(q4.st.al,q3.al); q4.st.al

#Klikemarginalerne
p.me.al=q1.me.al; p.me.al
p.ve.al=q2.ve.al; p.ve.al
p.an.al=q3.an.al; p.an.al
p.st.al=q4.st.al; p.st.al

