#Calculating correlation and regression for allele content prediction

#7 scenarios hr40 hr20 rand hr40x1 hr40x2 hr40y1 hr40y2
# 3 assumed heritabilities
# 5 alleles
# 2 groups


#Empty vectors
corr <- rep(0,42)
scen <- rep("tomt",42)
hop <- rep("tomt",42)
slope <- rep(0,42)
intercept <- rep(0,42)
trait <- rep("tomt",42)
her <- rep("tomt",42)
dir <- rep("")

#Effects from sampling to use in estRSS
teffect<- c(-1.18, 1.55, -3.07, -2.98, 0.14)   #1
teffect<- c(-1.29, 1.52, -3.02, -2.51, 0.16)   #2
teffect<- c(-1.18, 1.50, -2.44, -2.62, 0.14)   #3
teffect<- c(-1.27, 1.47, -1.97, -2.91, 0.17)   #4
teffect<- c(-1.26, 1.53, -2.96, -3.38, 0.16)   #5
teffect<- c(-1.47, 1.56, -1.93, -3.22, 0.20)   #6
teffect<- c(-1.52, 1.52, -1.88, -2.95, 0.21)   #7
teffect<- c(-1.16, 1.41, -2.98, -3.09, 0.15)   #8
teffect<- c(-1.47, 1.67, -3.34, -2.91, 0.19)   #9
teffect<- c(-1.01, 1.45, -2.67, -3.16, 0.12)   #10



teffect <- c(-1.27, 1.51, -2.51, -3.03, 0.17)   #used in trueRSS

effect <- c(-1.27, 1.51, -2.51, -3.03, 0.17)
effect <- matrix(effect)
teffect <- matrix(teffect)

filenames <- c("hr20/forcor_2021.txt", "hr20/forcor_2016.txt",
               "hr40/forcor_2021.txt", "hr40/forcor_2016.txt",
               "rand/forcor_2021.txt", "rand/forcor_2016.txt",
               "hr40x1/forcor_2021.txt", "hr40x1/forcor_2016.txt",
               "hr40x2/forcor_2021.txt", "hr40x2/forcor_2016.txt",
               "hr40y1/forcor_2021.txt", "hr40y1/forcor_2016.txt",
               "hr40y2/forcor_2021.txt", "hr40y2/forcor_2016.txt")

scenes <- c("hr20", "hr20", "hr40", "hr40",
             "rand", "rand", 
             "hr40x1", "hr40x1", "hr40x2", "hr40x2",
             "hr40y1", "hr40y1", "hr40y2", "hr40y2")

hopar <- c("2021", "2016_20", "2021", "2016_20", 
           "2021", "2016_20", "2021", "2016_20",
           "2021", "2016_20", "2021", "2016_20",
           "2021", "2016_20")



         
#Setting the working directory.
setwd("C:/Users/jonhjalti/OneDrive - Menntaský/Documents/riða/kynbo/runs10/rep1")
#Effects from sampling to use in estRSS_2 where 
teffect<- c(-1.01, 1.00, -2.76, -3.42, 0.15)   #1
teffect<- c(-2.49, 3.01, -2.90, -3.01, 0.30)   #2
teffect<- c(-0.71, 0.76, -2.81, -1.98, 0.11)   #3
teffect<- c(-2.12, 0.96, -4.21, -2.35, 0.36)   #4
teffect<- c(-0.81, 1.43, -2.90, -4.02, 0.09)   #5
teffect<- c(-2.62, 1.89, -1.72, -4.43, 0.40)   #6
teffect<- c(-1.30, 2.06, -0.80, -4.00, 0.13)   #7
teffect<- c(-0.78, 2.19, -3.00, -3.04, 0.04)   #8
teffect<- c(-0.82, -0.31, -6.67, -6.06, 0.20)   #9
teffect<- c(-1.32, 1.89, -5.24, -2.54, 0.15)   #10

a <- -3

#Main loop for each replication
for (i in 1:14){ 
  inn <- read.table(file=filenames[i])
  a <- a+3
  s <- a+1
  e <- a+3
  scen[s:e] <- scenes[i]
  hop[s:e] <- hopar[i]
  trait[s:e] <- c("logodds","logodds","logodds")
  her[s:e] <- c("90","95","99")
  
  rest <- 2-inn$V14-inn$V15-inn$V16-inn$V17
  mat <- cbind(inn$V14,inn$V15,inn$V16,inn$V17,rest)
  tbv <- mat%*%effect
  tbv <- tbv[,1] 
  
  rest <- 2-inn$V2-inn$V5-inn$V8-inn$V11
  mat <- cbind(inn$V2,inn$V5,inn$V8,inn$V11,rest)
  ebv90 <- mat%*%teffect
  ebv90 <- ebv90[,1] 
  
  rest <- 2-inn$V3-inn$V6-inn$V9-inn$V12
  mat <- cbind(inn$V3,inn$V6,inn$V9,inn$V12,rest)
  ebv95 <- mat%*%teffect
  ebv95 <- ebv95[,1] 
  
  rest <- 2-inn$V4-inn$V7-inn$V10-inn$V13
  mat <- cbind(inn$V4,inn$V7,inn$V10,inn$V13,rest)
  ebv99 <- mat%*%teffect
  ebv99 <- ebv99[,1] 
  
  corr[1+a] <- cor(ebv90,tbv)
  corr[2+a] <- cor(ebv95,tbv)
  corr[3+a] <- cor(ebv99,tbv)


  reg <- lm(tbv~ebv90)
  slope[a+1] <- reg$coefficients[2]
  intercept[a+1] <- reg$coefficients[1]
  reg <- lm(tbv~ebv95)
  slope[a+2] <- reg$coefficients[2]
  intercept[a+2] <- reg$coefficients[1]
  reg <- lm(tbv~ebv99)
  slope[a+3] <- reg$coefficients[2]
  intercept[a+3] <- reg$coefficients[1]

}


ut <- data.frame(scen,her,hop,trait,corr,slope,intercept)

#Writing out results for each replicate

#write.table(ut,file="result_logodds.txt",col.names = F, quote = F)
#File for all reps:
#write.table(ut,file="../result_logodds.txt",col.names = F, quote = F, append = TRUE)

#For writing out the estRSS scenario
#write.table(ut,file="result_wreff.txt",col.names = F, quote = F)
#File for all reps:
#write.table(ut,file="../result_wreff.txt",col.names = F, quote = F, append = TRUE)

#For writing out the estRSS2 scenario
write.table(ut,file="result_wreff2.txt",col.names = F, quote = F)
#File for all reps:
write.table(ut,file="../result_wreff2.txt",col.names = F, quote = F, append = TRUE)



#Change working directory
setwd("C:/Users/jonhjalti/OneDrive - Menntaský/Documents/riða/kynbo/runs10")


library(dbplyr)
library(ggplot2)
library(ggridges)
library(tidyverse)
library(cowplot)

#read the results from all replicates
ut <- read.table(file = "result_logodds.txt", colClasses = c("NULL", "factor", "factor", "factor", "factor", "numeric", "numeric", "numeric"))
ut <- read.table(file = "result_wreff.txt", colClasses = c("NULL", "factor", "factor", "factor", "factor", "numeric", "numeric", "numeric"))
ut <- read.table(file = "result_wreff2.txt", colClasses = c("NULL", "factor", "factor", "factor", "factor", "numeric", "numeric", "numeric"))


colnames(ut) = c("Scenario", "Heritability", "group", "trait", "correlation", "bias", "level bias")
#Making graps

utlo_16 <- filter(ut, group == "2016_20")
utlo_21 <- filter(ut, group == "2021")

utlo_16a <- utlo_16 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(correlation),
            se_corr =sd(correlation)/sqrt(10))


p1 = ggplot(utlo_16a,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.0,0.65))+
  ylab("Accuracy") +
  xlab(element_blank())+
  scale_fill_discrete(labels=c('0.90', '0.95', '0.99'))+
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  theme(legend.position = "none")


utlo_21a <- utlo_21 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(correlation),
            se_corr =sd(correlation)/sqrt(10))

p2 = ggplot(utlo_21a,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.0,0.65))+
  ylab(element_blank()) +
  xlab(element_blank())+
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  theme(legend.position = "none")

utlo_16b <- utlo_16 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(bias),
            se_corr =sd(bias)/sqrt(10))


pb1 = ggplot(utlo_16b,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.7,1.3))+
  ylab("Dispersion") +
  xlab(element_blank())+
  scale_fill_discrete(labels=c('0.90', '0.95', '0.99'))+
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  theme(legend.position = "none")


utlo_21b <- utlo_21 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(bias),
            se_corr =sd(bias)/sqrt(10))

pb2 = ggplot(utlo_21b,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.7,1.3))+
  ylab(element_blank()) +
  xlab(element_blank())+
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  theme(legend.position = "none")

legend <- get_legend(p1+
                       guides(fill = guide_legend(nrow = 1)) +
                       theme(legend.position = c(0.4,1), legend.direction = "horizontal"))

cplot <- plot_grid(p1, p2, pb1, pb2, nrow = 2, ncol = 2, 
                   labels = c("2016-2020", "2021",NULL, NULL), 
                   label_size = 13, label_x = 0.5)
plot_grid(cplot,legend, ncol=1, rel_heights = c(1, 0.1))


#Write the results to files
write.table(utlo_16a,file="result_acc_16r.txt",col.names = F, quote = F)
write.table(utlo_21a,file="result_acc_21r.txt",col.names = F, quote = F)

write.table(utlo_16b,file="result_bias_16r.txt",col.names = F, quote = F)
write.table(utlo_21b,file="result_bias_21r.txt",col.names = F, quote = F)

write.table(utlo_16a,file="result_acc_16v.txt",col.names = F, quote = F)
write.table(utlo_21a,file="result_acc_21v.txt",col.names = F, quote = F)

write.table(utlo_16b,file="result_bias_16v.txt",col.names = F, quote = F)
write.table(utlo_21b,file="result_bias_21v.txt",col.names = F, quote = F)
