library(tidyverse)
library(dplyr)
library(raster)
library(ggbo)
library(rstatix)
library(ggpubr)
library(ggsignif)
library(lme4)
library(MASS)
library(RColorBrewer)
library(afex)
library(psych)
################################################################################
#           FUNCTIONS                                                          #
################################################################################
mysummary = function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                     conf.interval=.95, .drop=TRUE) {
  library(plyr,dplyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
################################################################################
# VARIABLES                                                                    #
################################################################################
workingDirectoy="/media/marco/MarcoHDD/github/SOIMA/Test_model__MMR_MI/"
workingDirectory2="/media/marco/MarcoHDD/github/SOIMA/BMU_labeled/"
setwd(workingDirectoy)
correlations=Sys.glob("*.csv")#Activation Correlations from SOIMA
input=c("Baba","Gaga","Dada","McGurk","Baga","Bada","Gada","Daba","Daga")
columns=c("Baba","Gaga","Dada","McGurk","Baga","Bada","Gada","Daba","Daga","Input","Model")
columns_MI=c("MI","input","condition","model")
McGurk_columns=c("condition","model","proportion")
congruentStimuli=c("Baba","Gaga","Dada")
numberfiles=1:length(correlations)
colors= c("#003F5C", "#BC5090","#FFA600")
Comparisons=list(c("Ba","Ga"),c("Da","Ba"),c("Ga","Da"))
################################################################################
# MAIN                                                                         #
################################################################################
for (i in seq_along(numberfiles)){
  #print(i)
  if (i==1){
    correlationResult=read.csv(correlations[i],header=F)
    correlationResult=data.frame(correlationResult)
    #correlationResult["input"]=input
    #correlationResult["modelo"]=i
  }
  else{
    matrix1=correlationResult
    matrix2=read.csv(correlations[i],header=F)
    matrix2=data.frame((matrix2))
    #matrix2["input"]=input
    #matrix2["modelo"]=i
    correlationResult=rbind(matrix1,matrix2)
  }
}
names(correlationResult)=columns_MI
correlationResult$input=as.factor(correlationResult$input)
correlationResult$condition=as.factor(correlationResult$condition)
correlationResult$model=as.factor(correlationResult$model)



print("Charging BMU answers from McGurk-like stimuli")
setwd(workingDirectory2)
BMUs=Sys.glob("*.csv")#
numberfiles=1:length(BMUs)
for (i in seq_along(numberfiles)){
  #print(i)
  if (i==1){
    BMUResult=read.csv(BMUs[i],header=F)
    BMUResult=data.frame(BMUResult)
    BMUResult$model=i
    #correlationResult["input"]=input
    #correlationResult["modelo"]=i
  }
  else{
    matrix1=BMUResult
    matrix2=read.csv(BMUs[i],header=F)
    matrix2=data.frame((matrix2))
    matrix2$model=i
    BMUResult=rbind(matrix1,matrix2)
  }
}
print("Calculating Square Chi for BMU responses for McGurk-like stimuli")
McGurk_x2=chisq.test(BMUResult$V1,BMUResult$model)
print("Exploring to which congruent stimuli McGurk-like stimuli were asociated with")
table_=table(BMUResult$V1,BMUResult$model)
table_=table_/181
McGurk_data=as.data.frame(table_)
names(McGurk_data)=McGurk_columns
exploratory=kruskal.test(proportion~condition,data=McGurk_data)
wilcox_test(proportion~condition,data=McGurk_data, paired=T,detailed=T,p.adjust = "bonferroni")
plotData=mysummary(McGurk_data,measurevar="proportion",groupvars="condition")
print("Ploting results from BMU from McGurk-like stimuli")
bmu_plot=ggplot(plotData, aes(x=condition, y=proportion  )) + 
  geom_bar(position=position_dodge(), stat="identity",
           fill=colo[1:3],
           colour="black", # Use black outlines,
           size=.3,
           width = .5,) +      # Thinner lines
  geom_errorbar(aes(ymin=proportion-se, ymax=proportion+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  theme(text = element_text(size = 20))  +
  xlab("Congruent Percept") +
  ylab("Ocurrence proportion") +
  ggtitle("Ocurrence Proportion of McGurk stimuli being labeled as congruent percepts") +
  #scale_y_continuous(limits=(0,1)) +
  scale_x_discrete(labels = c("Ba", "Ga", "Da"))+
  theme_bw()
bmu_plot+ theme(axis.text.y = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16))+ theme(axis.title.y = element_text(size=16))+
  theme(axis.title.x = element_text(size=16))
################################################################################
### LMM with corr coef as dependant, input and condition as fixed effects      #
##  and model as random effect                                                 #
################################################################################
print("Contrast coding average mean")
contrast_i=contr.sdif(6)
contrast_c=contr.sdif(3)
#s=mixed(corr.coeficient~Input*Condition+(1|Model),data=correl_coef_lf)
#summary(s)
myt= lmer(MI~input*condition+(1|model),data=correlationResult,contrasts=list(input=contrast_i, condition=contrast_c))
summary(myt)
myd=aov(MI ~ input * condition + Error(1/(input*condition)), 
        data = correlationResult,observed="input",contrasts=list(input=contrast_i, condition=contrast_c))
summary(myd)

my_aov=aov(MI~input*condition+Error(input/condition),data=correlationResult,contrasts=list(input=contrast_i, condition=contrast_c))
summary(my_aov)




#correlationResult$input=as.factor(correlationResult$input)
s.dif()
Ba=correlationResult%>%filter(input=="Baba")
Ba=Ba["V4"]%>%
  mutate(class="Ba")
Ga=filter(correlationResult,input=="Gaga")
Ga=Ga["V4"]%>%
  mutate(class="Ga")
Da=filter(correlationResult,input=="Dada")
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
colo= brewer.pal(6,"Dark2")

an1=mysummary(correlationResult,measurevar="MI",groupvars="input")


ggplot(data = correlationResult, mapping = aes(x = input, y = MI,color=condition)) +
  geom_boxplot() +
  geom_rug(aes(col=condition),alpha=0.8)+
  #geom_jitter(position=position_jitter(0.1), shape = 16, alpha = 1, size = 0.5) +
  labs(x = "Tipo de estímulo Incongruente", y = "Información mutua (log)", size = 30) + 
  #stat_summary( geom = "point", position=position_dodge(width=0.9), shape=18, size=7) +
  scale_y_continuous(breaks = seq(0.075, 0.110, 0.01), limits=c(0.075, 0.110)) +
   theme(axis.title.x = element_text(color="black", size=12),
                          axis.title.y = element_text(color="black", size=20),
                          axis.ticks.length = unit(0.3,"cm"),
                          axis.ticks = element_line(size = 1),
                          axis.line = element_line(size = 0.3, linetype = "solid"),
                          axis.text = element_text (size = 10, colour = 'black'),
                          plot.margin = margin(1, 1, 1, 1, "cm"),
                          legend.position= "none") +
            #stat_boxplot(geom ='errorbar', lwd = 1, width = 0.15)+
  facet_wrap(vars(model),ncol=4)+
  labs( 
    x="Incongruent audiovisual input",
    y="Mutual Information (MI)",
    col="Comparation condition",
    title= "Descriptive plots from data conditions ",
    subtitle="MI~input*condition+(1|model)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))




myggp+ theme(axis.text.y = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 16))+ theme(axis.title.y = element_text(size=16))+
  theme(axis.title.x = element_text(size=16))
##############################
c=mysummary(correlationResult,measurevar="MI",groupvars = c("input","condition","model"))

factor1_n=length(levels(correlationResult$input))
factor2_n=length(levels(correlationResult$condition))
colors=brewer.pal(factor1_n*factor2_n,"Dark2")
rt_log10_boxplot.plot = ggplot(c, aes(x = input, y = MI, fill = condition)) +
  geom_boxplot() +
  geom_point(aes(color = viridis(factor1_n*factor2_n)))
################################
print("coding fromm dummy to categorical")
correlationResult$condition=revalue(correlationResult$condition, c("1"="Ba", "2"="Ga","3"="Da"))
correlationResult$input=revalue(correlationResult$input, c("4"="McGurk", "5"="BaGa",
                                   "6"="BaDa","7"="GaDa",
                                   "8"="DaBa","9"="DaGa"))
correlationResult$model=revalue(correlationResult$model, c("0"="SOIMA_1", "1"="SOIMA_2",
                                                           "2"="SOIMA_3","3"="SOIMA_4",
                                                           "4"="SOIMA_5","5"="SOIMA_6",
                                                           "6"="SOIMA_7","7"="SOIMA_8",
                                                           "8"="SOIMA_9","9"="SOIMA_10")) 
  

fixed1=filter(correlationResult,input==c("McGurk","BaGa"))
constante=1.5
levels(fixed1$input)=c("-.5",".5","0","0","0","0")
fixed1["num_input"]=as.numeric(fixed1$input)
d=fixed1$input
fixed1$num_input=fixed1$num_input-constante
fixed2=filter(correlationResult,input==c("BaDa","GaDa"))
levels(fixed2$input)=c("-.5",".5","0","0","0","0")
fixed2["num_input"]=as.numeric(fixed2$input)-3.5
fixed3=filter(correlationResult,input==c("DaGa","DaBa"))
#levels(fixed3$input)=as.numeric(fixed3$input)
d=as.numeric(fixed3$input)
fixed3$input=d-constante-4
fixed4=filter(correlationResult,condition==c("Ga","Ba"))
#levels(fixed4$condition)=as.numeric(fixed4$condition)
d=as.numeric(fixed4$condition)
fixed4$condition=d-constante
fixed5=filter(correlationResult,condition==c("Da","Ga"))
#levels(fixed5$condition)=as.numeric(fixed5$condition)
d=as.numeric(fixed5$condition)
fixed5$condition=d-constante-1
#plot1#
plot_fixed1=ggplot(fixed1,aes(x = num_input,y=MI))+
  geom_boxplot(aes(group=num_input),width=0.5,alpha=1, lwd = .4) +
  geom_rug(aes(col=condition),alpha=0.8)+
  geom_jitter(aes( col=condition),alpha=0.4,
              position = position_jitter(width = .3, height=-0.7),
              size=.3) +
  geom_abline(intercept = myt@beta[1], slope = myt@beta[2])+
  scale_colour_manual(values = brewer.pal(6,"Dark2"))+
  scale_x_continuous( breaks = c(-.5,.5),labels = c("McGurk", "BaGa"))+
  theme_linedraw()+
  labs( 
    x="Incongruent audiovisual input",
    y="Mutual Information (MI)",
    col="Comparation condition",
    title= "Mutual information for significant fixed effect",
    subtitle="MI~input*condition+(1|model)")
 plot_fixed1+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#plot2#
 plot_fixed2=ggplot(fixed2,aes(x = num_input,y=MI))+
   geom_boxplot(aes(group=num_input),width=0.5,alpha=1, lwd = .4) +
   geom_rug(aes(col=condition),alpha=0.8)+
   geom_jitter(aes( col=condition),alpha=0.4,
               position = position_jitter(width = .3, height=-0.7),
               size=.3) +
   geom_abline(intercept = myt@beta[1], slope = myt@beta[4])+
   scale_colour_manual(values = brewer.pal(6,"Dark2"))+
   scale_x_continuous(breaks = c(-.5,.5), labels = c("GaDa","BaDa"))+
   theme_linedraw()+
   labs( 
     x="Incongruent audiovisual input",
     y="Mutual Information (MI)",
     col="Comparation condition",
     title= "Mutual information for significant fixed effect",
     subtitle="MI~input*condition+(1|model)")
 plot_fixed2+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 
 #plot3#
 plot_fixed3=ggplot(fixed3,aes(x = input,y=MI))+
   geom_boxplot(aes(group=input),width=0.5,alpha=1, lwd = .4) +
   geom_rug(aes(col=condition),alpha=0.8)+
   geom_jitter(aes( col=condition),alpha=0.4,
               position = position_jitter(width = .3, height=-0.7),
               size=.3) +
   geom_abline(intercept = myt@beta[1], slope = myt@beta[6])+
   scale_colour_manual(values = brewer.pal(6,"Dark2"))+
   scale_x_continuous(breaks = c(-.5,.5), labels = c("DaBa","DaGa"))+
   theme_linedraw()+
   labs( 
     x="Incongruent audiovisual input",
     y="Mutual Information (MI)",
     col="Comparation condition",
     title= "Mutual information for significant fixed effect",
     subtitle="MI~input*condition+(1|model)")
 plot_fixed3+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 #plot4#
 plot_fixed4=ggplot(fixed4,aes(x = condition,y=MI))+
   geom_boxplot(aes(group=condition),width=0.5,alpha=1, lwd = .4) +
   geom_rug(aes(col=input),alpha=0.8)+
   geom_jitter(aes( col=input),alpha=0.4,
               position = position_jitter(width = .3, height=-0.7),
               size=.3) +
   geom_abline(intercept = myt@beta[1], slope = myt@beta[7])+
   scale_colour_manual(values = brewer.pal(6,"Paired"))+
   scale_x_continuous(breaks = c(-.5,.5), labels = c("Ba","Ga"))+
   theme_linedraw()+
   labs( 
     x="Comparation condition",
     y="Mutual Information (MI)",
     col="Input condition",
     title= "Mutual information for significant fixed effect",
     subtitle="MI~input*condition+(1|model)")
 plot_fixed4+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 #plot5#
 plot_fixed5=ggplot(fixed5,aes(x = condition,y=MI))+
   geom_boxplot(aes(group=condition),width=0.5,alpha=1, lwd = .4) +
   geom_rug(aes(col=input),alpha=0.8)+
   geom_jitter(aes( col=input),alpha=0.4,
               position = position_jitter(width = .3, height=-0.7),
               size=.3) +
   geom_abline(intercept = myt@beta[1], slope = myt@beta[8])+
   scale_colour_manual(values = brewer.pal(6,"Paired"))+
   scale_x_continuous(breaks = c(-.5,.5), labels = c("Ga","Da"))+
   theme_linedraw()+
   labs( 
     x="Input condition",
     y="Mutual Information (MI)",
     col="Incongruent audiovisual input",
     title= "Mutual information for significant fixed effect",
     subtitle="MI~input*condition+(1|model)")
 plot_fixed5+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
