# Este algoritmo calcula el poligono que encierra los puntos en el espacio
# Autor: Fernando Espindola R.
# Fecha: 10/10/2006
# Modificado (Optimizado): Juan Carlos Saavedra (24/08/2009)
# -------------------------------------------------------------------------
genera.polyJC <- function(long,lati,nn,dt){
  min_lon <- min(long)
  max_lon <- max(long)
  min_lat <- min(lati)
  max_lat <- max(lati)
  dl_lon <- diff(range(long))/nn
  dl_lat <- diff(range(lati))/nn
  latt <- min_lat
  a1 <- which(lati==latt)
  lonn <- long[a1]
  pos <- cbind(lonn,latt)
  for(i in 1:nn){
    line <- latt + dl_lat*i
    a2 <- which(lati>=line-dl_lat*dt & lati<=line+dl_lat*dt)
    if(length(a2)>0){
      temp <- cbind(long[a2],lati[a2])
      posi <- which(temp[,1]==min(long[a2]))
      a3 <- a2[posi]
      new <- cbind(long[a3],lati[a3])
      pos <- rbind(pos,new)
    }
  }
  for(j in 1:nn){
    lina <- max_lat-dl_lat*(j-1)
    a4 <- which(lati<=lina+dl_lat*dt & lati>=lina-dl_lat*dt)
    if(length(a4)>0){
      junk <- cbind(long[a4],lati[a4])
      posy <- which(junk[,1]==max(long[a4]))
      a5 <- a4[posy]
      new <- cbind(long[a5],lati[a5])
      pos <- rbind(pos,new)
    }
  }
  pos <- rbind(pos,pos[1,])
  return(pos)
}
