#install.packages("meta")
library(meta)
genotypes <- c("arr_arr", "arr_ahq", "arr_arq", "arr_vrq", 
               "ahq_ahq", "ahq_arq", "ahq_vrq", "vrq_arq",
               "vrq_vrq","t137_arq")
filenames <- c("arr_arr.txt", "arr_ahq.txt", "arr_arq.txt", "arr_vrq.txt", 
               "ahq_ahq.txt", "ahq_arq.txt", "ahq_vrq.txt", "vrq_arq.txt",
               "vrq_vrq.txt","t137_arq.txt")

feffect <- rep(0,10)
reffect <- rep(0,10)

upperf <- rep(0,10)
lowerf <- rep(0,10)
sef <- rep(0,10)

for (i in 1:10)  {
  inn<-read.table(file=filenames[i])
  til = inn$V3+inn$V2
  vid = inn$V4+inn$V5

  metam<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", method = "MH")
  feffect[i] <- metam$TE.common  
  reffect[i] <- metam$TE.random
  upperf[i] <- metam$upper.common
  lowerf[i] <- metam$lower.common
  sef[i] <- metam$seTE.common
}  

ut <- data.frame(genotypes, reffect, feffect, sef, lowerf, upperf)

write.table(ut,file="logoddsratio.txt",col.names = F, quote = F)




##Ekki í lykkju
  
arr_arr<-read.table(file='arr_arr.txt')

inn<-arr_arr
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_arr_arr<-metabin(event.e = inn$V2, n.e = til,
             event.c = inn$V4, n.c = vid,
             sm = "OR", method = "Inverse")
effects[1] <- m_arr_arr$TE.common

arr_ahq<-read.table(file='arr_ahq.txt')

inn<-arr_ahq
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_arr_ahq<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", method = "MH")


effects[2] <- m_arr_ahq$TE.common

#arr/arq
arr_arq<-read.table(file='arr_arq.txt')

inn<-arr_arq
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_arr_arq<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", method = "MH")

effects[3] <- m_arr_arq$TE.common

#arr/vrq
arr_vrq<-read.table(file='arr_vrq.txt')

inn<-arr_vrq
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_arr_vrq<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", method = "MH")

effects[4] <- m_arr_vrq$TE.common

#ahq/ahq
ahq_ahq<-read.table(file='ahq_ahq.txt')

inn<-ahq_ahq
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_ahq_ahq<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", method = "MH")

effects[5] <- m_ahq_ahq$TE.common

#ahq/arq
ahq_arq<-read.table(file='ahq_arq.txt')

inn<-ahq_arq
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_ahq_arq<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", method = "MH")

effects[6] <- m_ahq_arq$TE.common

#AHQ/VRQ
ahq_vrq<-read.table(file='ahq_vrq.txt')

inn<-ahq_vrq
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_ahq_vrq<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", nethod = "MH")

effects[7] <- m_ahq_vrq$TE.common

#VRQ/ARQ
vrq_arq<-read.table(file='vrq_arq.txt')

inn<-vrq_arq
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_vrq_arq<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", nethod = "MH")

effects[8] <- m_vrq_arq$TE.common

#VRQ/VRQ
vrq_vrq<-read.table(file='vrq_vrq.txt')

inn<-vrq_vrq
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_vrq_vrq<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", nethod = "MH")

effects[9] <- m_vrq_vrq$TE.common

#T137/ARQ
t137_arq<-read.table(file='t137_arq.txt')

inn <- t137_arq
til = inn$V3+inn$V2
vid = inn$V4+inn$V5

m_t137_arq<-metabin(event.e = inn$V2, n.e = til,
                   event.c = inn$V4, n.c = vid,
                   sm = "OR", nethod = "MH")

effects[10] <- m_t137_arq$TE.common

ut <- data.frame(genotypes, effects)

write.table(ut,file="logoddsratio.txt",col.names = F, quote = F)
