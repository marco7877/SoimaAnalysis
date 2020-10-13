library(tidyverse)
library(dplyr)
library(raster)
library(ggbo)
library(rstatix)
library(ggpubr)
library(ggsignif)
################################################################################
#           FUNCTIONS                                                          #
################################################################################
mysummary = function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr,dplyr)
  
  length2 = function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  

  datac = ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  
  datac = rename(datac, c("mean" = measurevar))
  
  datac$se = datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  

  ciMult = qt(conf.interval/2 + .5, datac$N-1)
  datac$ci = datac$se * ciMult
  
  return(datac)
}
################################################################################
# VARIABLES                                                                    #
################################################################################
workingDirectoy="/media/marco/MarcoHDD/github/SOIMA/MMR_testingSet(indovidual)/"
setwd(workingDirectoy)
correlations=Sys.glob("*.csv")#Activation Correlations from SOIMA
input=c("Baba","Gaga","Dada","McGurk","Baga","Bada","Gada","Daba","Daga")
input1=c("Baba","Gaga","Dada","McGurk")
numberfiles=1:length(correlations)
colors= c("#003F5C", "#BC5090","#FFA600")
Comparisons=list(c("Ba","Ga"),c("Da","Ba"),c("Ga","Da"))
################################################################################
# MAIN                                                                         #
################################################################################
for (i in seq_along(numberfiles)){
  print(i)
  if (i==1){
    correlationResult=read.csv(correlations[i],header=F)
    correlationResult=data.frame(correlationResult)
    correlationResult=cbind(correlationResult,input1)
  }
  else{
    matrix1=correlationResult
    matrix2=read.csv(correlations[i],header=F)
    matrix2=data.frame((matrix2))
    matrix2=cbind(matrix2,input1)
    correlationResult=rbind(matrix1,matrix2)
  }
}
#correlationResult$input=as.factor(correlationResult$input)
Ba=correlationResult%>%filter(input1=="Baba")
Ba=Ba["V4"]%>%
  mutate(class="Ba")
Ga=filter(correlationResult,input1=="Gaga")
Ga=Ga["V4"]%>%
  mutate(class="Ga")
Da=filter(correlationResult,input1=="Dada")
Da=Da["V4"]%>%
  mutate(class="Da")
result=rbind(Ba,Ga)
result=rbind(result,Da)
summary(aov(V4~class,data=result))
as.numeric(result$V4)
as.factor(result$class)
t_paired_test_corrected_Mcgurk=pairwise.t.test(result$V4, 
             result$class, 
             p.adjust = "bonferroni")
print(t_paired_test_corrected_Mcgurk)
stat_data=mysummary(result,measurevar = "V4",groupvars = "class")
print("Plotting group boxplot")
ggboxplot(result, x = "class", y = "V4", 
                       color = "class",
                       ylab = "Correlation with McGurk Stimuli", xlab = "Congruent syllable")
print("Plotting group correlations as bar charts with error bar")

mybarplot=ggplot(stat_data, aes(x=class, y=V4  )) + 
  geom_bar(position=position_dodge(), stat="identity",
           fill=colors,
           colour="black", 
           size=.3,
           width = .5,) +      
  geom_errorbar(aes(ymin=V4-se, ymax=V4+se),
                size=.3,    
                width=.2,
                position=position_dodge(.9)) +
  theme(text = element_text(size = 20))  +
  xlab("Congruent Syllable") +
  ylab("Correlation to McGurk MMR activation pattern") +
  ggtitle("Comparition between McGurk stimuli and congruent stimuli") +
  theme_bw()
mybarplot+ theme(axis.text.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.title.x = element_text(size=16))
######



plot=ggplot(result, aes(x=class, y=V4)) + 
  geom_boxplot(fill=colors) +
  xlab("Congruent Syllable") +
  ylab("Correlation to McGurk MMR activation pattern") +
  ggtitle("Comparition between McGurk stimuli and congruent stimuli") +
  geom_signif(comparisons = Comparisons, 
              map_signif_level=TRUE)+
  theme_bw()

plot+ theme(axis.text.y = element_text(size = 16)) +
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.title.x = element_text(size=16))
