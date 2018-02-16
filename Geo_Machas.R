##########################################################
# ANALISIS GEOSTADISTICO DE TAGELUS DOMBEII
##########################################################


rm(list=ls(all=TRUE))

#install.packages("geoRglm")
#install.packages("lattice")
#install.packages("splancs")
#install.packages("PBSmapping")
#install.packages("polyclip")
#install.packages("geoR")
#install.packages("coda")

#Llama las librerias necesarias para ejecutar este anlisis
library(geoRglm) # A Package for Generalized linear spatial models
library(coda)    # MCMC diagnostics package
library(splancs) # Spatial point pattern analysis
library(PBSmapping) # Data mapping, polygon creation and LL to UTM conversion
library(polyclip) # Necessary for buffering the polygon
library(geoR) 


setwd( "C:/Users/mauricio.mardones/Documents/IFOP/Prog_Seg_PM_Bentonicos_2016/Geo_Machas" );dir()
dat <- read.csv("Evadir_Godoy.csv", header=TRUE, sep=";")
names(dat)
max(dat$Tdom)
#pasar datos de densidad a 1m^2
#dat$tagelus<-(dat$tagelus/0.25)*1 ;head(dat$tagelus)
#dat$X<-dat$X/1000
#dat$Y<-dat$Y/1000

tag0<-subset(dat,Tdom>0&Tdom<=5)
tag5<-subset(dat,Tdom>5&Tdom<=10)
tag15<-subset(dat,Tdom>10&Tdom<=20)
tag45<-subset(dat,Tdom>20)

x11()
plot(dat$X,dat$Y,ylab="Latitud (km)",xlab="Longitud (km)",main="Densidad de Godoy",cex.main=0.8,cex=1.3)
#points(dat$X,dat$Y,ylab="Latitud (km)",xlab="Longitud (km)",pch=19,cex=1.5)
points(tag0$X,tag0$Y,ylab="Latitud (km)",xlab="Longitud (km)",pch=19,cex=1,col=4)
points(tag5$X,tag5$Y,ylab="Latitud (km)",xlab="Longitud (km)",pch=19,cex=1,col=3)
points(tag15$X,tag15$Y,ylab="Latitud (km)",xlab="Longitud (km)",pch=19,cex=1,col="orange")
points(tag45$X,tag45$Y,ylab="Latitud (km)",xlab="Longitud (km)",pch=19,cex=1,col=2)
text(611600,5394000,"N° individuos/m2",cex=0.8)
legend(611600,5394000,c("1 - 5","5 - 10","10 - 20",">20"),pch=19,col=c("blue","green","orange","red"),cex=0.8,bty="n",pt.cex=1.5)
#=====================================================================================
#POLIGONO DEL AREA DE EVALUACIÓN

#=====================================================================================

# 2.	Poligono para la estimacion

#=====================================================================================

datos<-data.frame(PID=rep(1,length(dat[,1])),POS=seq(1,length(dat[,1]),1),cuadrante=dat$Tdom,E=dat$X,N=dat$Y)
datos.pos<-datos[datos$cuadrante>0,,]


ntot   <-length(datos[,1]);ntot #numero de observaciones totales
npos   <-length(datos.pos[,1]);npos #numero de observaciones positivas
meanpos<-mean(datos.pos$cuadrante);meanpos #densdad media solo datos positivos
#====================================================================================
# Poligono de datos totales
#====================================================================================
source("C:/Users/mauricio.mardones/Documents/IFOP/Prog_Seg_PM_Bentonicos_2016/Geo_Machas/generapolyJC.R") #llama a función que está en la carpeta Geostat
polig   <- genera.polyJC(datos$E,datos$N,850,20)
polig.1 <- data.frame(X=polig[,1],Y=polig[,2])
convexH <- polig.1
# Obtains buffered convex hull para generar bordes
polyH            <- list(x=convexH$X, y=convexH$Y)
aux              <- polyoffset(polyH, .1, jointype='round')
poly.H           <- data.frame(cbind(aux[[1]]$x, aux[[1]]$y));colnames(poly.H) <- c('X', 'Y')
write.table(poly.H, file='blanca.poly.csv', col.names=TRUE, sep=',')

aux.1 <- 1:length(poly.H[,1])
aux.2 <- rep(1,length(poly.H[,1]))
poly.H <- cbind(aux.2, aux.1, poly.H);colnames(poly.H) <- c('PID','POS','X','Y')
poly.H <- as.PolySet(poly.H); class(poly.H)
poly.H <- fixBound(poly.H, tol=0.00001)
poly.H <- closePolys(poly.H)

#====================================================================================
# Poligono de datos positivos
#====================================================================================
poligp   <- genera.polyJC(datos.pos$E,datos.pos$N,150,12)
poligp.1 <- data.frame(X=poligp[,1],Y=poligp[,2])
convexHp <- poligp.1

# Obtains buffered convex hull
polyHposp            <- list(x=convexHp$X, y=convexHp$Y)
aux6                 <- polyoffset(polyHposp, .1, jointype='round')
poly.Hposp           <- data.frame(cbind(aux6[[1]]$x, aux6[[1]]$y));colnames(poly.Hposp) <- c('X', 'Y')
write.table(poly.Hposp, file='blanca.polypos.csv', col.names=TRUE, sep=',')

aux.1p <- 1:length(poly.Hposp[,1])
aux.2p <- rep(1,length(poly.Hposp[,1]))
poly.Hposp <- cbind(aux.2p, aux.1p, poly.Hposp);colnames(poly.Hposp) <- c('PID','POS','X','Y')
poly.Hposp <- as.PolySet(poly.Hposp); class(poly.Hposp)
poly.Hposp <- fixBound(poly.Hposp, tol=0.00001)
poly.Hposp <- closePolys(poly.Hposp)

#*************************************************************************************
# Estimacion del area de cada poliono
#*************************************************************************************
poligono      <- read.table(file= 'blanca.poly.csv', header=TRUE, sep=',')
poligono.pos  <- read.table(file= 'blanca.polypos.csv', header=TRUE, sep=',')
poltot        <- data.frame(PID=rep(1,length(poligono[,1])),POS=seq(1,length(poligono[,1]),1),X=poligono$X,Y=poligono$Y)
polefect      <- data.frame(PID=rep(1,length(poligono.pos[,1])),POS=seq(1,length(poligono.pos[,1]),1),X=poligono.pos$X,Y=poligono.pos$Y)

# Area total de muestreo y Area efectiva en Km^2
Areatotal<-calcArea(poltot);Areatotal[,2] 
Areaefect<-calcArea(polefect);Areaefect[,2] 


#*************************************************************************************
# Gráfico de poligonos
#*************************************************************************************

x11()
plot(convexH$X ,convexH$Y, type='l', xlab='Easting (km)', ylab='Northing (km)')
polygon(poly.H$X,poly.H$Y , col='grey')
polygon(poly.Hposp$X,poly.Hposp$Y,col=3)
points(datos$E,datos$N)
points(datos.pos$E,datos.pos$N,pch=19,cex=1.2)

#=====================================================================================

# 3.- Estimacion de la probabilidad de observar el stock

#=====================================================================================
#Establece los rangos de Easting y Northing para confeccionar una grilla
#La grilla es necesaria para agrupar los datos. 
#Esto se hace porque el analisis puntual (proceso bernoulli) es computacionalmente intenso
#Se requiere el tamaÃ±o de la celda en Km^2

head(datos)  # datos de densidad
head(poly.H) # poligono

celda <- 100

min_east  <- floor(min(datos[4]))-celda;min_east                         
max_east  <- ceiling(max(datos[4]))+celda;max_east 
min_north <- floor(min(datos[5]))-celda;min_north
max_north <- ceiling(max(datos[5]))+celda;max_north

# Genera los rangos de Easting y Northing de la grilla
east  <- seq(from = min_east, to =max_east, by= celda) 
north <- seq(from = min_north, to =max_north, by= celda)

# Crea una matriz que contendrÃ¡ la grilla
grilla <- matrix(NA,((length(east)-1)*(length(north)-1)), 7)
head(grilla)
#-----------------------------------------------------------------------------------------------------------
# Llena la matriz con los datos de Easting y Northing generando 
#simultaneamente la grilla
h <- 1

for (i in 1: (length(east)-1)) {
	for (j in 1: (length(north)-1)) {
		grilla[h,]<-c(east[i], east[i+1],north[j], north[j+1],(east[i]+east[i+1])/2,(north[j]+ north[j+1])/2,h)  
		h <- h + 1
	}
}
head(grilla)
#-----------------------------------------------------------------------------------------------------------
# Visualiza la grilla
grilla <- as.data.frame(grilla)
names(grilla) <- c("Easting_ini", "Easting_fin", "Northing_ini", "Northing_fin","Easting_avg","Northing_avg","ID")
#-----------------------------------------------------------------------------------------------------------
#Recorre la grilla contruyendo una tabla con el numero de observaciones 
#totales (unidades) y el número de observaciones positivas(data)que ocurren por cuadricula 
frec <- mat.or.vec(length(grilla[,1]),2)

for (i in 1: length(grilla[,1])){
	for (j in 1: length(datos[,1])){         
		if(datos[j,4] > grilla[i,1] & datos[j,4] <= grilla[i,2]){ 
			if (datos[j,5] > grilla[i,3] & datos[j,5] <= grilla[i,4]){
				frec[i,1] <- frec[i,1] + 1
				if (datos[j,3] > 0) frec[i,2] <-  frec[i,2] + 1	}           
		} }}
#-----------------------------------------------------------------------------------------------------------
dim(grilla)
head(grilla)
head(frec)

bindata.1<-cbind(grilla[,5],grilla[,6], frec[,1],frec[,2])
length(datos[,1]) - sum(datos[,3] == 0) # numero de observaciones positivas
sum(bindata.1[,4]) # numero de observaciones positivas
head(bindata.1)
tail(bindata.1)

#Elimina celdas no observadas
eval      <-bindata.1[,3]>0
bindata.2 <- bindata.1[eval,]
length(bindata.2[,3])

#Crea datos clase geodata y caracteriza los datos
bindata.3 <-as.geodata(bindata.2,coords.col= 1:2, units.m.col= 3,data.col=4)
summary(bindata.3)

#resultados
n.trials  <- sum(bindata.3$units.m);n.trials
n.succ    <- sum(bindata.3$data);n.succ
naive_est <- n.succ/n.trials;naive_est

#-------------------------------------------------------------------------------------------------------------
#Ploteo de datos binarios:
x11()
par(mfrow=c(2,2))
par(mar=c(3.8, 5, 1.1, 2.))
par(mgp=c(2.2,1,0))

plot(datos$E,datos$N, pch=20, xlab='Easting (km)', ylab='Northing (km)')
#plot(grilla$Easting_avg,grilla$Northing_avg, col="red", xlab="Easting (km)", ylab="Northing (km)", type='n')
plot(poly.H$X, poly.H$Y, xlab="Easting (km)", ylab="Northing (km)", type='n')
abline(v=(seq(min(east),max(east),.5)), col="lightgray", lty="dotted")
abline(h=(seq(min(north), max(north),.5)), col="lightgray", lty="dotted")
points(grilla$Easting_avg, grilla$Northing_avg, pch=1, col='green3')
points(datos$E,datos$N, pch=20)
lines(convexH$X, convexH$Y)
lines(poly.H$X, poly.H$Y, col='red')
#Plotea cuadriculas con observaciones indicando el numero de ensayos
plot(bindata.3$coords[,1],bindata.3$coords[,2],type="n",xlab="Easting (km)",ylab="Northing (km)",main="Ensayos")
text(bindata.3$coords[,1],bindata.3$coords[,2],format(bindata.3$units.m),cex=.7)
abline(v=(seq(min(east),max(east),celda)), col="lightgray", lty="dotted")
abline(h=(seq(min(north), max(north),celda)), col="lightgray", lty="dotted")    
rect(bindata.3$coords[,1]-50,bindata.3$coords[,2]-50,bindata.3$coords[,1]+50,bindata.3$coords[,2]+50)
points(datos$E,datos$N, pch=20)
#Plotea cuadriculas con observaciones indicando el numero de exitos
plot(bindata.3$coords[,1],bindata.3$coords[,2],type="n",xlab="Easting (km)",ylab="",main="Exitos")
text(bindata.3$coords[,1],bindata.3$coords[,2],format(bindata.3$data),cex=.7)
abline(v=(seq(min(east),max(east),celda)), col="lightgray", lty="dotted")
abline(h=(seq(min(north), max(north),celda)), col="lightgray", lty="dotted")    
rect(bindata.3$coords[,1]-50,bindata.3$coords[,2]-50,bindata.3$coords[,1]+50,bindata.3$coords[,2]+50)
#-------------------------------------------------------------------------------------------------------------

# La funcion de verosimilitud aqui no tiene una forma matematica conocida, queda expresada como una doble
# integral lo que se debe solucionar mediante MCMC.
# El beta es el beta promedio del puntaje logit se debe dejar igual a uno.
# El rango debe ser parecido a (pero algo mayor) al estimado al calcular la densidad media con las observaciones positivas
# De acuerdo con Christensen y Ribeiro, 2006 (An introducition to the package geoRglm), el parametro nugget al 
# igual que los parametros de suvizado (kappa) y de anisotropia si los ubiese pueden ser ctes.
# De acuerdo con estos autores se debe sintonizar el algoritmo escalando la varianza propuesta 
# de modo que la tasa de acceptacion sea aproximadamente 60% (0.6) valor optimo para el algoritmo empleado.
# Esto se hace mediante ensayo y error
# Un thin=10 corresponde a 10000 iteraciones

bindata.3.spmod <- list(cov.pars=c(0.434,5.192454), beta=1.0, cov.model="matern", nugget=0.2170433, kappa=5, family="binomial", link="logit")
bindata.3.mcmc  <- mcmc.control(S.scale=0.68, thin=20, burn.in=1500, n.iter=20000)
bindata.3.tune  <- glsm.mcmc(bindata.3, model= bindata.3.spmod, mcmc.input= bindata.3.mcmc)

# Tambien debemos estudiar la convergencia de la cadena y que tan bien la cadena se esta mezclando
# Para ello empleamos las funciones del paquete CODA
# Se recomienda realizar graficos trace y de autocorrelacion
X11()
bindata.3.tune.coda <- create.mcmc.coda(bindata.3.tune, mcmc.input= bindata.3.mcmc)
autocorr.plot(bindata.3.tune.coda, ask=TRUE)

# MLE
bindata.3.pre.lf <- prepare.likfit.glsm(bindata.3.tune)
bindata.3.lf     <- likfit.glsm(bindata.3.pre.lf,ini.phi=0.1, hessian=TRUE, nugget.rel=0.2170433, fix.nugget.rel=TRUE, cov.model="matern", kappa=5)
summary(bindata.3.lf)
bindata.3.lf.p   <- exp(bindata.3.lf$beta)/(1+exp(bindata.3.lf$beta)) #??? revisar
exp(bindata.3.lf$parameters.summary[[2]][2])/(1+exp(bindata.3.lf$parameters.summary[[2]][2]))#??? revisar


###########################################################################################
# kriggin MEJORAR ESTA PARTE DEL CODIGO
###########################################################################################

#Prop.Tag<-data.frame(E=bindata.3$coords[,1],N=bindata.3$coords[,2],data=bindata.3$data)
#write.table(Prop.Tag, file='geo_proptag.csv', sep=',', row.names=FALSE) 
#P.tagelus<-read.geodata('geo_proptag.csv', sep=",", header= TRUE)


#Tp.pred.grid <- expand.grid(seq(634,639,l=500),seq(5579,5589,l=500))

# Ejecuta el kriging
#Tp.krig<-krige.conv(P.tagelus,locations=Tp.pred.grid, krige=krige.control(obj.m=bindata.3.lf), borders=poligono)

#x11()
# Plotea el kriging
#image(tagelus.krig, loc=tagelus.pred.grid,col=terrain.colors(100),add=T, xlab="Easting (km)",ylab="Northing (km)")
#contour(tagelus.krig, tagelus.pred.grid, levels = seq(0, 145, by = 15),	add = TRUE)

############################################################################################
############################################################################################
############################################################################################

# ANÁLISIS DE LAS OBSERVACIONES CON DENSIDAD POSITIVA

############################################################################################
############################################################################################
############################################################################################

# Escribe y lee los datos hacia y desde un archivo de texto y los prepara para un analisis geoestadistico con GeoR (creando un objeto geodata)
aux <- datos[,c(4,5,3)]
a1  <- aux$cuadrant==0; sum(a1)
aux <- aux[!a1,] #solo datos positivos
length(aux[,1])

write.table(aux, file='geo_blanca.csv', sep=',', row.names=FALSE) #guarda solo datos positivos para analisis
G.tagelus<-read.geodata('geo_blanca.csv', sep=",", header= TRUE)

# Revisa la posibilidad de que existan datos con coordenadas duplicadas
dup.coords(G.tagelus)

# Analisis Descriptivo
summary(G.tagelus)
x11()
plot(G.tagelus, bor=poligono) 

# Los datos no se distribuyen normales y se requiere una #transformacion de la familia Box-Cox lo que es controlado por el parametro "lamda"
plot(G.tagelus, bor=poligono, lambda= 0.1)

# Variograma de nube (el argumento max.dist se obtuvo del valor de distancia maxima que entrega el summary)
tagelus.cloud<-variog(G.tagelus,option="cloud",max.dist=2500,lambda=0.1)
x11()
 plot(tagelus.cloud)

# Variograma de bins

#seq(0,5): esto es el rango de distancia en el cual se calcularan los bins;
# l: es el numero de bins

# Aqui puedo jugar con el numero de bins y el rango de X
# ej: seq(0,5) o seq(0,6) o seq(0,75), para visualizar mejor
# los parametros del variograma

tagelus_vbins<-variog(G.tagelus,uvec=seq(0,2500,l=100),lambda=0.1)
plot(tagelus_vbins,yl="Semi-variance",xl="Distance (km)")


#Ajuste del variograma teorico al experimental mediante minimos cuadrados
# ini.cov.pars=c(30,1) es equivalente a c(Sill, Rango), nugget es nugget

tagelus_variofit <-variofit(tagelus_vbins, ini.cov.pars=c(2,3), 
	cov.model = "gauss",
	fix.nugget = FALSE, nugget = 1,
	simul.number = NULL, max.dist = tagelus_vbins$max.dist,
	limits = pars.limits())

summary(tagelus_variofit)

x11()
plot(tagelus_vbins,yl="Semi-variance",xl="Distance (km)")
lines.variomodel(cov.model=tagelus_variofit$cov.model,
	cov.pars=tagelus_variofit$cov.pars,nugget=tagelus_variofit$nugget,
	max.dist=tagelus_variofit$max.dist,lwd=1,col="red")


# La funcion variog4, genera de manera automatica variogramas  en cuatro direcciones arbitrarias,
# lo que permite verificar el supuesto de isotrop  los graficos deberan mostrar una tendencia 
# similar de variograma.
tagelus_vaniso<-variog4(G.tagelus,uvec=seq(0,2500,l=100),lambda=0.1)
plot(tagelus_vaniso,col=c("black","red","blue","green"),yl="Semi-variance", xl="Distance (km)")
#--------------------------------------------------------------------------------------------
# MLE
#--------------------------------------------------------------------------------------------
#Emplearemos los paraetros estimados en la seccion anterior para ejecutar una estimacion 
#maximo verosimil.
#kappa determina el suavizado. cuando es 0,5 es como exponencial; 
#cuando tiene un numero grande (5) esto es como una s, cuando es 1 es Gauss; 
#el Hessiano sera necesario para calcular la incertidumbre de los parametros. 
#Nota: par(10,1): par(sigma,rango), el parametro tau es equivalente al nugget

tagelus.lf <- likfit(G.tagelus,cov.model="matern", ini.cov.pars=c(0.434,5.192454),kappa=0.7, fix.kappa=FALSE,nugget=0.2170433,lambda=0.6,fix.lambda=FALSE,hessian=TRUE)
summary(tagelus.lf)
names(tagelus.lf)

#si no converge me paso a un modelo de correlacion espacial mas simple como el gaussiano
#tagelus.lf <-likfit(G.tagelus,cov.model="gaussian", ini.cov.pars=tagelus_variofit$cov.pars, lambda=0.1, fix.lambda=FALSE,hessian=TRUE,fix.psiA=T, fix.psiR=T)


#tagelus.lf <- likfit(G.tagelus,cov.model="gaussian",ini.cov.pars=c(0.8,3),
#	nugget=0.5,lambda=0.1,fix.lambda=FALSE,hessian=TRUE,fix.psiA=T, fix.psiR=T)
#summary(tagelus.lf)
#names(tagelus.lf)

#para obtener la media en la escala de la transformacion Box-Cox
tagelus.lf$beta

#para obtener la media en la escala original
tagelus.lf.mean<-mean(BCtransform(rnorm(5000,mean=tagelus.lf$beta, sd=sqrt(tagelus.lf$beta.var)), lambda=tagelus.lf$lambda,inverse=TRUE)$data)
tagelus.lf.mean

#Para obtener la varianza de la estimacion de la media del proceso espacial
#Ojo: en escala de la tranformacion Box-Cox
tagelus.lf$beta.var

#para obtener la varianza de estimacion de la media del proceso espacial 
#se ejecuta el comando siguiente. Este emplea la curavatura del perfil de verosimilitud
tagelus.lf.sd      <- sd(BCtransform(rnorm(5000,mean=tagelus.lf$beta, sd=sqrt(tagelus.lf$beta.var)), 
								      lambda=tagelus.lf$lambda,inverse=TRUE)$data);tagelus.lf.sd
#para obtener el coeficiente de variacion hacer lo siguiente:
tagelus.lf.beta.cv <- sqrt(tagelus.lf$beta.var)/tagelus.lf$beta; tagelus.lf.beta.cv
#para obtener las varianzas de los parametros necesito 
#la inversa de la matriz hessiana (las segundas derivadas parciales de los parametros de la funcion de verosimilitud) 
#uso el siguiente comando ...
solve(tagelus.lf$info.minimisation.function$hessian)
#para obtener las desviaciones estandares de los parametros. 
#En el modelo gaussiano la diagonal me entrega 
#en la posicion (1,1) la varianza del rango,
#en la posicion (2,2)la razon entre la varianza  sigma cuadrado(sill)sobre tau cuadrado y
#en la posicion (3,3) la varianza de lambda
sdphi     <- sqrt(solve(tagelus.lf$info.minimisation.function$hessian)[1,1]);sdphi
sdsig_tau <- sqrt(solve(tagelus.lf$info.minimisation.function$hessian)[2,2]);sdsig_tau

#==================================================================================================
# ATENCION
#la instruccion siguiente se usa solo cuando el modelo de correlacion espacial es el de matern
sdkappa<-sqrt(solve(tagelus.lf$info.minimisation.function$hessian)[3,3])
#sdkappa
#en ese caso la varianza de lambda esta en la posicion [4,4]
sdlambda<-sqrt(solve(tagelus.lf$info.minimisation.function$hessian)[4,4])
sdlambda
#==================================================================================================
# Construye una grilla para poder ejecutar el kriging
tagelus.pred.grid <- expand.grid(seq(602000,607500,l=50),seq(5396000,5400000,l=50))
#pred.grid.in <- locations.inside(tagelus.pred.grid,poligono)


#data.grill<-data.frame(dens_T=tagelus.krig$predict,X=pred.grid.in[,1],Y=pred.grid.in[,2])
tagelus.krig<-krige.conv(G.tagelus,locations=tagelus.pred.grid, krige=krige.control(obj.m=tagelus.lf), borders=poligono)
tagelus.krig2<-krige.conv(G.tagelus,locations=tagelus.pred.grid, krige=krige.control(obj.m=tagelus.lf), borders=poligono.pos)


# Densidad
x11()
# Plotea el kriging
image(tagelus.krig2, loc=tagelus.pred.grid ,col=gray(seq(1,0.1,l=5)), xlab="Easting (km)",ylab="Northing (km)")
contour(tagelus.krig2, tagelus.pred.grid, levels = seq(0, 145, by = 15),	add = TRUE)
#points(datos$E,datos$N,pch="+")
#polygon(poligono.pos)


image(tagelus.krig, loc=tagelus.pred.grid ,col=gray(seq(1,0.1,l=5)), xlab="Easting (km)",ylab="Northing (km)")
contour(tagelus.krig, tagelus.pred.grid, levels = seq(0, 145, by = 15),	add = TRUE)
#points(datos$E,datos$N,pch="+")
#polygon(poli
polygon(poligono)
# Densidad media transformada donde estaba el stock Se obtiene una aproximacion desde el kriging cuando el comando de trasformacion
# a la escala original de la densidad desde la transformacion Box-Cox no funciona.
tagelus.krig.mean <- mean(tagelus.krig$predict);tagelus.krig.mean
tagelus.krig.mean2 <- mean(tagelus.krig2$predict);tagelus.krig.mean2
# Error estandar de la densidad media transformada donde estaba el stock
tagelus.sd.mean <- tagelus.lf.beta.cv * tagelus.krig.mean;tagelus.sd.mean
tagelus.sd.mean2 <- tagelus.lf.beta.cv * tagelus.krig.mean2;tagelus.sd.mean2
# @@@ Extrapola densidad media al area total @@@



# Empleando la densidad media del krigging

biomasa1 <- tagelus.krig.mean * Areatotal[,2] * bindata.3.lf.p
biomasa1.sd <- tagelus.sd.mean * Areatotal[,2] * bindata.3.lf.p

biomasa12 <- tagelus.krig.mean2 * Areaefect[,2] * bindata.3.lf.p
biomasa12.sd <- tagelus.sd.mean2 * Areaefect[,2] * bindata.3.lf.p


# IC:
cbind(biomasa1 - 1.96 * biomasa1.sd, biomasa1, biomasa1 + 1.96 * biomasa1.sd)
cbind(biomasa12 - 1.96 * biomasa12.sd, biomasa12, biomasa12 + 1.96 * biomasa12.sd)


