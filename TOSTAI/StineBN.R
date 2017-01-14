#Opgave beskrivelse
dg <- dag(~ smoke + lung|smoke + xray|lung:smoke + bronc|smoke + dysp|bronc:lung )
plot(dg)

data(chestSim1000, package="gRbase")
head(chestSim1000)

#Opgave 2.1:
#CPT'er lavet ud fra data
s <- xtabs(~smoke, chestSim1000); s
x_ls <- xtabs(~xray+lung+smoke, chestSim1000); x_ls
l_s <- xtabs(~lung+smoke, chestSim1000); l_s
d_lb <- xtabs(~dysp+lung+bronc, chestSim1000); d_lb
b_s <- xtabs(~bronc+smoke, chestSim1000); b_s
#normaliser
p.s <- arnormalize(s, "first"); p.s
p.x_ls <- arnormalize(x_ls, "first"); p.x_ls
p.l_s <- arnormalize(l_s, "first"); p.l_s
p.d_lb <- arnormalize(d_lb, "first"); p.d_lb
p.b_s <- arnormalize(b_s, "first"); p.b_s
#Netværk laves ud fra CPT'erne
cpt.list <- compileCPT(list(p.s, p.x_ls, p.l_s, p.d_lb, p.b_s))
cpt.list

bn <- grain(cpt.list)


#Opgave 2.2
#p.s_d
qgrain(bn, nodes = c("smoke", "dysp"), type="conditional")

#Opgave 2.3
#p.l_s
qgrain(bn, nodes = c("lung", "smoke"), type="conditional")
#p.l_sb
qgrain(bn, nodes = c("lung", "smoke", "bronc"), type="conditional")
#Ud fra dette kan jeg se at bronc ikke giver mere viden om lung, da fordelingen ikke ændrer sig. 
#Altså er lung og brunc betinget uafhængige (når der betinges med smoke).
#Dette kan bekræftes hvis vi regner CPR = den ene diag / den anden diag = 1.


#Opgave 2.4
#p.l_sd
qgrain(bn, nodes = c("lung", "smoke", "dysp"), type="conditional")
#p.l_sdb
qgrain(bn, nodes = c("lung", "smoke", "dysp", "bronc"), type="conditional")
#Ud fra dette kan jeg se at vi får mere information om lung hvis vi kender bronc. 
#Så bronc og lung er ikke betinget uafhængige givet både smoke og dysp.

