#Calculating correlation and regression for allele content prediction

#7 scenarios hr40 hr20 rand hr40x1 hr40x2 hr40y1 hr40y2
# 3 assumed heritabilities
# 5 alleles
# 2 groups

#Empty vectors
corr <- rep(0,210)
scen <- rep("tomt",210)
hop <- rep("tomt",210)
slope <- rep(0,210)
intercept <- rep(0,210)
trait <- rep("tomt",210)
her <- rep("tomt",210)

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

#Groups
hopar <- c("2021", "2016_20", "2021", "2016_20", 
           "2021", "2016_20", "2021", "2016_20",
           "2021", "2016_20", "2021", "2016_20",
           "2021", "2016_20")

a <- -15

#Set the working directory for each replication
setwd("C:/Users/jonhjalti/OneDrive - Menntaský/Documents/riða/kynbo/runs10/rep10")


#Main loop for reading the files and calculating correlation and regression within replication
for (i in 1:14){ 
  inn <- read.table(file=filenames[i])
  a <- a+15
  s <- a+1
  e <- a+15
  scen[s:e] <- scenes[i]
  hop[s:e] <- hopar[i]
  trait[s:e] <- c("S1","S1","S1","S2","S2","S2","S3","S3","S3",
                  "S4","S4","S4","S5","S5","S5")
  her[s:e] <- c("90","95","99")
  corr[1+a] <- cor(inn$V2,inn$V14)
  corr[2+a] <- cor(inn$V3,inn$V14)
  corr[3+a] <- cor(inn$V4,inn$V14)
  corr[4+a] <- cor(inn$V5,inn$V15)
  corr[5+a] <- cor(inn$V6,inn$V15)
  corr[6+a] <- cor(inn$V7,inn$V15)
  corr[7+a] <- cor(inn$V8,inn$V16)
  corr[8+a] <- cor(inn$V9,inn$V16)
  corr[9+a] <- cor(inn$V10,inn$V16)
  corr[10+a] <- cor(inn$V11,inn$V17)
  corr[11+a] <- cor(inn$V12,inn$V17)
  corr[12+a] <- cor(inn$V13,inn$V17)

  tbv <- 2-inn$V14-inn$V15-inn$V16-inn$V17
  ebv90 <- 2-inn$V2-inn$V5-inn$V8-inn$V11
  corr[13+a] <- cor(tbv,ebv90)
  ebv95 <- 2-inn$V3-inn$V6-inn$V9-inn$V12
  corr[14+a] <- cor(tbv,ebv95)
  ebv99 <- 2-inn$V4-inn$V7-inn$V10-inn$V13
  corr[15+a] <- cor(tbv,ebv99)

  reg <- lm(inn$V14~inn$V2)
  slope[a+1] <- reg$coefficients[2]
  intercept[a+1] <- reg$coefficients[1]
  reg <- lm(inn$V14~inn$V3)
  slope[a+2] <- reg$coefficients[2]
  intercept[a+2] <- reg$coefficients[1]
  reg <- lm(inn$V14~inn$V4)
  slope[a+3] <- reg$coefficients[2]
  intercept[a+3] <- reg$coefficients[1]

  reg <- lm(inn$V15~inn$V5)
  slope[a+4] <- reg$coefficients[2]
  intercept[a+4] <- reg$coefficients[1]
  reg <- lm(inn$V15~inn$V6)
  slope[a+5] <- reg$coefficients[2]
  intercept[a+5] <- reg$coefficients[1]
  reg <- lm(inn$V15~inn$V7)
  slope[a+6] <- reg$coefficients[2]
  intercept[a+6] <- reg$coefficients[1]

  reg <- lm(inn$V16~inn$V8)
  slope[a+7] <- reg$coefficients[2]
  intercept[a+7] <- reg$coefficients[1]
  reg <- lm(inn$V16~inn$V9)
  slope[a+8] <- reg$coefficients[2]
  intercept[a+8] <- reg$coefficients[1]
  reg <- lm(inn$V16~inn$V10)
  slope[a+9] <- reg$coefficients[2]
  intercept[a+9] <- reg$coefficients[1]

  reg <- lm(inn$V17~inn$V11)
  slope[a+10] <- reg$coefficients[2]
  intercept[a+10] <- reg$coefficients[1]
  reg <- lm(inn$V17~inn$V12)
  slope[a+11] <- reg$coefficients[2]
  intercept[a+11] <- reg$coefficients[1]
  reg <- lm(inn$V17~inn$V13)
  slope[a+12] <- reg$coefficients[2]
  intercept[a+12] <- reg$coefficients[1]

  reg <- lm(tbv~ebv90)
  slope[a+13] <- reg$coefficients[2]
  intercept[a+13] <- reg$coefficients[1]
  reg <- lm(tbv~ebv95)
  slope[a+14] <- reg$coefficients[2]
  intercept[a+14] <- reg$coefficients[1]
  reg <- lm(tbv~ebv99)
  slope[a+15] <- reg$coefficients[2]
  intercept[a+15] <- reg$coefficients[1]

}

#Write the results to file
ut <- data.frame(scen,her,hop,trait,corr,slope,intercept)

write.table(ut,file="result_table.txt",col.names = F, quote = F)
#File for all results across reps
write.table(ut,file="../result_table.txt",col.names = F, quote = F, append = TRUE)



##The rest is for making graphs

library(dbplyr)
library(ggplot2)
library(ggridges)
library(tidyverse)
library(cowplot)

ut <- read.table(file = "result_table.txt", colClasses = c("NULL", "factor", "factor", "factor", "factor", "numeric", "numeric", "numeric"))

colnames(ut) = c("Scenario", "Heritability", "group", "trait", "correlation", "bias", "level bias")

#First for the A1 allele

uts1 <- filter(ut, trait %in% c("S1"))

uts1_16 <- filter(uts1, group == "2016_20")

uts1_16a <- uts1_16 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(correlation),
            se_corr =sd(correlation)/sqrt(10))
  

p1 = ggplot(uts1_16a,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.0,0.65))+
  ylab("Accuracy") +
  xlab(element_blank())+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"),
                    labels=c('0.90', '0.95', '0.99'))+  
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  theme(legend.position = "none")

uts1_21 <- filter(uts1, group == "2021")

uts1_21a <- uts1_21 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(correlation),
            se_corr =sd(correlation)/sqrt(10))

p2 = ggplot(uts1_21a,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.0,0.65))+
  ylab(element_blank()) +
  xlab(element_blank())+
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"))+
  theme(legend.position = "none")


uts1_16b <- uts1_16 %>%
  group_by(Scenario, Heritability) %>%
  summarise(m_bias = mean(bias),
            se_corr =sd(bias)/sqrt(10))


pb1 = ggplot(uts1_16b,
            aes(x = Scenario, y = m_bias, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=m_bias-se_corr, ymax=m_bias+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.7,1.3))+
  ylab("Dispersion") +
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"))+
  theme(legend.position = "none")


uts1_21b <- uts1_21 %>%
  group_by(Scenario, Heritability) %>%
  summarise(m_bias = mean(bias),
            se_corr =sd(bias)/sqrt(10))

pb2 = ggplot(uts1_21b,
            aes(x = Scenario, y = m_bias, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=m_bias-se_corr, ymax=m_bias+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.7,1.3))+
  ylab(element_blank()) +
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"))+
  theme(legend.position = "none")


legend <- get_legend(p1+
                       guides(fill = guide_legend(nrow = 1)) +
                       theme(legend.position = c(0.4,1), legend.direction = "horizontal"))

cplot <- plot_grid(p1, p2, pb1, pb2, nrow = 2, ncol = 2, 
                   labels = c("2016-2020", "2021",NULL, NULL), 
                   label_size = 13, label_x = 0.5)
plot_grid(cplot,legend, ncol=1, rel_heights = c(1, 0.1))

#Writing results to file
write.table(uts1_16a,file="result_s1acc_16r.txt",col.names = F, quote = F)
write.table(uts1_21a,file="result_s1acc_21r.txt",col.names = F, quote = F)

write.table(uts1_16b,file="result_s1bias_16r.txt",col.names = F, quote = F)
write.table(uts1_21b,file="result_s1bias_21r.txt",col.names = F, quote = F)

#Now for the rare allele A4 
                                                          
uts4 <- filter(ut, trait %in% c("S4"))
uts4_16 <- filter(uts4, group == "2016_20")

uts4_16a <- uts4_16 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(correlation),
            se_corr =sd(correlation)/sqrt(10))


p1 = ggplot(uts4_16a,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  #                colour = Sviðsmynd,
  #                linetype = Arfgerð)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.0,0.65))+
  ylab("Accuracy") +
  xlab(element_blank())+
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"),
                    labels=c('0.90', '0.95', '0.99'))+  
  theme(legend.position = "none")

#p1
uts4_21 <- filter(uts4, group == "2021")

uts4_21a <- uts4_21 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(correlation),
            se_corr =sd(correlation)/sqrt(10))

p2 = ggplot(uts4_21a,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.0,0.65))+
  ylab(element_blank()) +
  xlab(element_blank())+
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"))+
  theme(legend.position = "none")


uts4_16b <- uts4_16 %>%
  group_by(Scenario, Heritability) %>%
  summarise(m_bias = mean(bias),
            se_corr =sd(bias)/sqrt(10))


pb1 = ggplot(uts4_16b,
             aes(x = Scenario, y = m_bias, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=m_bias-se_corr, ymax=m_bias+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.7,1.3))+
  ylab("Dispersion") +
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"))+
  theme(legend.position = "none")


uts4_21b <- uts4_21 %>%
  group_by(Scenario, Heritability) %>%
  summarise(m_bias = mean(bias),
            se_corr =sd(bias)/sqrt(10))

pb2 = ggplot(uts4_21b,
             aes(x = Scenario, y = m_bias, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=m_bias-se_corr, ymax=m_bias+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.7,1.3))+
  ylab(element_blank()) +
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"))+
  theme(legend.position = "none")

legend <- get_legend(p1+
                       guides(fill = guide_legend(nrow = 1)) +
                       theme(legend.position = c(0.4,1), legend.direction = "horizontal"))

cplot <- plot_grid(p1, p2, pb1, pb2, nrow = 2, ncol = 2, 
          labels = c("2016-2020", "2021",NULL, NULL), 
          label_size = 13, label_x = 0.5)
plot_grid(cplot,legend, ncol=1, rel_heights = c(1, 0.1))

write.table(uts4_16a,file="result_s4acc_16r.txt",col.names = F, quote = F)
write.table(uts4_21a,file="result_s4acc_21r.txt",col.names = F, quote = F)

write.table(uts4_16b,file="result_s4bias_16r.txt",col.names = F, quote = F)
write.table(uts4_21b,file="result_s4bias_21r.txt",col.names = F, quote = F)


### Supplimentary graphs


uts4 <- filter(ut, trait %in% c("S5"))


uts4_16 <- filter(uts1, group == "2016_20")

uts4_16a <- uts4_16 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(correlation),
            se_corr =sd(correlation)/sqrt(10))


p1 = ggplot(uts4_16a,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.0,0.65))+
  ylab("Accuracy") +
  xlab(element_blank())+
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"),
                    labels=c('0.90', '0.95', '0.99'))+  
  theme(legend.position = "none")


uts4_21 <- filter(uts4, group == "2021")

uts4_21a <- uts4_21 %>%
  group_by(Scenario, Heritability) %>%
  summarise(accuracy = mean(correlation),
            se_corr =sd(correlation)/sqrt(10))

p2 = ggplot(uts4_21a,
            aes(x = Scenario, y = accuracy, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=accuracy-se_corr, ymax=accuracy+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.0,0.65))+
  ylab(element_blank()) +
  xlab(element_blank())+
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"))+
  theme(legend.position = "none")

uts4_16b <- uts4_16 %>%
  group_by(Scenario, Heritability) %>%
  summarise(m_bias = mean(bias),
            se_corr =sd(bias)/sqrt(10))


pb1 = ggplot(uts4_16b,
             aes(x = Scenario, y = m_bias, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=m_bias-se_corr, ymax=m_bias+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.7,1.3))+
  ylab("Dispersion") +
  #  ggtitle("Þróun arfgerða allt landið") +
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"))+
  theme(legend.position = "none")


uts4_21b <- uts4_21 %>%
  group_by(Scenario, Heritability) %>%
  summarise(m_bias = mean(bias),
            se_corr =sd(bias)/sqrt(10))

pb2 = ggplot(uts4_21b,
             aes(x = Scenario, y = m_bias, fill = Heritability))+ 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=m_bias-se_corr, ymax=m_bias+se_corr), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim = c(0.7,1.3))+
  ylab(element_blank()) +
  theme_cowplot(13)+
  theme(panel.grid.major.y = element_line(colour="gray"))+
  scale_fill_manual(values=c("darkblue",
                             "forestgreen",
                             "darkorange"))+
  theme(legend.position = "none")

legend <- get_legend(p1+
                       guides(fill = guide_legend(nrow = 1)) +
                       theme(legend.position = c(0.4,1), legend.direction = "horizontal"))

cplot <- plot_grid(p1, p2, pb1, pb2, nrow = 2, ncol = 2, 
                   labels = c("2016-2020", "2021",NULL, NULL), 
                   label_size = 13, label_x = 0.5)
plot_grid(cplot,legend, ncol=1, rel_heights = c(1, 0.1))

