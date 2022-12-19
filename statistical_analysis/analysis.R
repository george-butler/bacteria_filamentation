library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)
library(chisq.posthoc.test)
library(cowplot)
library(grid)
library(gridExtra)
library(lsmeans)
library(extrafont)
########################################################################################################################################################
align_legend <- function(p, hjust = 0.5)
{
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    legend$grobs[[gi]] <- guides
  }
  
  g$grobs[[legend_index]] <- legend
  g
}

holder<-read.csv("~/bacteria_filamentation/statistical_analysis/final_data_set.csv")
holder$treatment<-as.character(holder$treatment)
holder$MIC<-0
holder[holder$treatment == 30,]$MIC<-0.5
holder[holder$treatment == 60,]$MIC<-1
holder$MIC<-as.character(holder$MIC)

########################################################################################################################################################
#test proportion of cells in each treatment that have a longer nucleoid vs filament length 
holder_temp<-holder[holder$treatment != 0,]
holder_temp$diff<-holder_temp$max_chromosome_length/holder_temp$max_filament_length
holder_temp$nuc_ext<-"N"
holder_temp[holder_temp$diff > 1,]$nuc_ext<-"Y"
tab1<-table(holder_temp$treatment, holder_temp$nuc_ext)
chisq.test(tab1, correct=TRUE)
########################################################################################################################################################
#filament length
holder_temp<-holder

holder_temp$treatment<-as.factor(holder_temp$treatment)
filament_length_aov<-lmer(log(max_filament_length) ~ treatment + (1|uni_id), data=holder_temp, REML=TRUE)
summary(filament_length_aov)
anova(filament_length_aov)

mult_comp<-lsmeans(filament_length_aov,list(pairwise~treatment),adjust = "bonferroni")
mult_comp

holder_temp$mean_filament_length<-NA
holder_temp$se_filament_length<-NA
mod_coef<-summary(filament_length_aov)[["coefficients"]]
holder_temp[holder_temp$treatment == 0,]$mean_filament_length[1]<-mod_coef[1,1]
holder_temp[holder_temp$treatment == 0,]$se_filament_length[1]<-mod_coef[1,2]

holder_temp$treatment<-relevel(holder_temp$treatment,2)
filament_length_aov<-lmer(log(max_filament_length) ~ treatment + (1|uni_id), data=holder_temp, REML=TRUE)
mod_coef<-summary(filament_length_aov)[["coefficients"]]
holder_temp[holder_temp$treatment == 30,]$mean_filament_length[1]<-mod_coef[1,1]
holder_temp[holder_temp$treatment == 30,]$se_filament_length[1]<-mod_coef[1,2]

holder_temp$treatment<-relevel(holder_temp$treatment,3)
filament_length_aov<-lmer(log(max_filament_length) ~ treatment + (1|uni_id), data=holder_temp, REML=TRUE)
mod_coef<-summary(filament_length_aov)[["coefficients"]]
holder_temp[holder_temp$treatment == 60,]$mean_filament_length[1]<-mod_coef[1,1]
holder_temp[holder_temp$treatment == 60,]$se_filament_length[1]<-mod_coef[1,2]

holder_temp$MIC<-as.character(holder_temp$MIC)
holder_temp[holder_temp$MIC == 0,]$MIC<-"Control"
holder_temp[holder_temp$MIC == 0.5,]$MIC<-"0.5xMIC"
holder_temp[holder_temp$MIC == 1,]$MIC<-"1xMIC"
holder_temp$MIC<-as.factor(holder_temp$MIC)
holder_temp$MIC<-factor(holder_temp$MIC, levels=c("Control","0.5xMIC","1xMIC"))
#MIC Figure
ggplot(holder_temp, aes(x=factor(MIC), y=log(max_filament_length), color=factor(MIC)))+
  geom_point(data=holder_temp,alpha=0.1,color="black",size=4)+
  geom_errorbar(aes(ymax = mean_filament_length + 1.96*se_filament_length, ymin = mean_filament_length - 1.96*se_filament_length),width=0.2,size=1.25)+
  geom_point(aes(y=mean_filament_length), size = 5)+scale_colour_manual(values = c("black","#3182bd","#e6550d"))+
  labs(color = "Ciprofloxacin \nconcentration (ng)")+
  xlab(expression(paste("xMIC"))) +  ylab(expression(paste("Maximum filament length (",mu,"m)")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 30,color = "black",family="Ariel"),
        axis.title = element_text(size = 30,color = "black",family="Ariel"),
        legend.text = element_text(size = 30,color = "black",family="Ariel"),
        legend.title = element_text(size = 30,color = "black",family="Ariel"),
        strip.text = element_text(size= 30,color = "black",family="Ariel"),
        axis.title.x=element_blank(),
        legend.position = "none")+
  scale_y_continuous(breaks=c(0,1.098612,1.945910,2.995732,4.007333,5.010635,5.991465), labels=c(1,3,7,20,55,150,400), limits = c(0,6))

########################################################################################################################################################
#Nucleoid length
holder_temp<-holder

holder_temp$treatment<-as.factor(holder_temp$treatment)
chromo_length_aov<-lmer(log(max_chromosome_length) ~ treatment + (1|uni_id), data=holder_temp, REML=TRUE)
summary(chromo_length_aov)
anova(chromo_length_aov)

mult_comp<-lsmeans(chromo_length_aov,list(pairwise~treatment),adjust = "bonferroni")

holder_temp$mean_chromo_length<-NA
holder_temp$se_chromo_length<-NA
mod_coef<-summary(chromo_length_aov)[["coefficients"]]
holder_temp[holder_temp$treatment == 0,]$mean_chromo_length[1]<-mod_coef[1,1]
holder_temp[holder_temp$treatment == 0,]$se_chromo_length[1]<-mod_coef[1,2]

holder_temp$treatment<-relevel(holder_temp$treatment,2)
chromo_length_aov<-lmer(log(max_chromosome_length) ~ treatment + (1|uni_id), data=holder_temp, REML=TRUE)
mod_coef<-summary(chromo_length_aov)[["coefficients"]]
holder_temp[holder_temp$treatment == 30,]$mean_chromo_length[1]<-mod_coef[1,1]
holder_temp[holder_temp$treatment == 30,]$se_chromo_length[1]<-mod_coef[1,2]

holder_temp$treatment<-relevel(holder_temp$treatment,3)
chromo_length_aov<-lmer(log(max_chromosome_length) ~ treatment + (1|uni_id), data=holder_temp, REML=TRUE)
mod_coef<-summary(chromo_length_aov)[["coefficients"]]
holder_temp[holder_temp$treatment == 60,]$mean_chromo_length[1]<-mod_coef[1,1]
holder_temp[holder_temp$treatment == 60,]$se_chromo_length[1]<-mod_coef[1,2]

holder_temp$MIC<-as.character(holder_temp$MIC)
holder_temp[holder_temp$MIC == 0,]$MIC<-"Control"
holder_temp[holder_temp$MIC == 0.5,]$MIC<-"0.5xMIC"
holder_temp[holder_temp$MIC == 1,]$MIC<-"1xMIC"
holder_temp$MIC<-as.factor(holder_temp$MIC)
holder_temp$MIC<-factor(holder_temp$MIC, levels=c("Control","0.5xMIC","1xMIC"))

ggplot(holder_temp, aes(x=factor(MIC), y=log(max_chromosome_length), color=factor(MIC)))+
  geom_point(data=holder_temp,alpha=0.1,color="black",size=4)+
  geom_violin(data=holder_temp[(holder_temp$lqt_cutoff == F) & (holder_temp$uqt_cutoff == F),],trim = FALSE)+
  geom_errorbar(aes(ymax = mean_chromo_length + 1.96*se_chromo_length, ymin = mean_chromo_length - 1.96*se_chromo_length),width=0.2,size=1.25)+
  geom_point(aes(y=mean_chromo_length), size = 5)+scale_colour_manual(values = c("black","#3182bd","#e6550d"))+
  labs(color = "Ciprofloxacin \nconcentration (ng)")+
  xlab(expression(paste("xMIC"))) +  ylab(expression(paste("Maximum nucleoid length (",mu,"m)")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 30,color = "black",family="Ariel"),
        axis.title = element_text(size = 30,color = "black",family="Ariel"),
        legend.text = element_text(size = 30,color = "black",family="Ariel"),
        legend.title = element_text(size = 30,color = "black",family="Ariel"),
        strip.text = element_text(size= 30,color = "black",family="Ariel"),
        axis.title.x=element_blank(),
        legend.position = "none")+
  scale_y_continuous(breaks=c(0,1.098612,1.945910,2.995732,4.007333,5.010635,5.991465), labels=c(1,3,7,20,55,150,400), limits = c(0,6))

########################################################################################################################################################
#Filament length vs the chromosome length
holder_temp<-holder[holder$treatment != 0,]

m1<-lmer(log(max_chromosome_length) ~ treatment + (1|uni_id), data = holder_temp)
m2<-lmer(log(max_chromosome_length) ~ log(max_filament_length):treatment + (1|uni_id), data = holder_temp)
anova(m1,m2)
summary(m2)
r.squaredGLMM(m2)

f_model<-lmer(log(max_chromosome_length) ~ log(max_filament_length):treatment+(1|uni_id), data = holder_temp)
mod_sum<-summary(f_model)
r.squaredGLMM(f_model)
holder_temp$y_pred<-NA
holder_temp[holder_temp$treatment == 30,]$y_pred<-mod_sum[["coefficients"]][1,1]+(mod_sum[["coefficients"]][2,1]*log(holder_temp[holder_temp$treatment == 30,]$max_filament_length))
holder_temp[holder_temp$treatment == 60,]$y_pred<-(mod_sum[["coefficients"]][1,1])+((mod_sum[["coefficients"]][3,1])*log(holder_temp[holder_temp$treatment == 60,]$max_filament_length))

holder_temp[holder_temp$MIC == 0.5,]$MIC<-"0.5xMIC"
holder_temp[holder_temp$MIC == 1,]$MIC<-"1xMIC"
p<-ggplot(holder_temp, aes(x = log(max_filament_length),y = y_pred, color=MIC))+
  geom_point(data= holder_temp, aes(x=log(max_filament_length),y=log(max_chromosome_length), color=MIC),size=4)+
  geom_line(size=3)+scale_colour_manual(values = c("#3182bd","#e6550d"))+
  labs(color = "xMIC")+
  xlab(expression(paste("Maximum filament length (",mu,"m)"))) + ylab(expression(paste("Maximum nucleoid length (",mu,"m)")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 30,colour = "black",family="Ariel"),
        axis.title = element_text(size = 30,colour = "black",family="Ariel"),
        legend.text = element_text(size = 30,colour = "black",family="Ariel"),
        legend.title = element_blank(),
        strip.text = element_text(size= 30,colour = "black",family="Ariel"),
        legend.position = c(0.2,0.9),
        legend.background = element_rect(fill = "#ffffffaa", colour = NA),
        legend.direction = "horizontal")+
  guides(color=guide_legend(ncol=1))+
  scale_y_continuous(breaks=c(1.098612,1.945910,2.995732,4.007333,5.010635,5.991465), labels=c(3,7,20,55,150,400), limits = c(1,6))+
  scale_x_continuous(breaks=c(1.098612,1.945910,2.995732,4.007333,5.010635,5.991465), labels=c(3,7,20,55,150,400), limits = c(1,6))

ggdraw(align_legend(p))
########################################################################################################################################################
#max number of proteins in each cell 
holder_temp<-holder
holder_temp[holder_temp$max_protein_load > 4,]$max_protein_load<-5

holder_temp[holder_temp$treatment == "0",]$treatment<-"Control"
holder_temp[holder_temp$treatment == "30",]$treatment<-"30ng"
holder_temp[holder_temp$treatment == "60",]$treatment<-"60ng"

holder_temp$treatment<-as.factor(holder_temp$treatment)
holder_temp$treatment<-factor(holder_temp$treatment, levels=c("Control","30ng","60ng"))

m1 = glmer(max_protein_load ~ 1 + (1|uni_id), family="poisson",data=holder_temp)
m2 = glmer(max_protein_load ~ treatment + (1|uni_id), family="poisson", data=holder_temp)
anova(m1, m2, test="Chisq")
summary(m2)

x_axis<- c("0","1","2","3","4","\u2265 5")

ordered_level<-levels(holder_temp$treatment)
poss_sig = array(data = NA, c(0,(length(ordered_level) + 1)))
for ( i in 1:length(ordered_level)){
  j1<-holder_temp
  j1$treatment<-relevel(j1$treatment, i)
  m1 = glmer(max_protein_load ~ treatment + (1|uni_id), family="poisson", data=j1)
  model_sum<-summary(m1)
  model_coef<-model_sum[["coefficients"]]
  b1_significance<-model_coef[c(1:3),4]
  h<-c()
  for(j in 1:length(ordered_level)){
    if ((i+1) == j){
      h<-append(h,1)
    }
    h<-append(h,b1_significance[j])
  }
  if((i == length(ordered_level)) & (j == length(ordered_level))){
    h<-append(h,1)
  }
  poss_sig  = rbind(poss_sig ,h)
}
poss_sig<-as.data.frame(poss_sig)
colnames(poss_sig)<-c("0",ordered_level)
rownames(poss_sig)<-ordered_level
poss_sig<-t(poss_sig)


holder_temp[holder_temp$MIC == "0",]$MIC<-"Control"
holder_temp[holder_temp$MIC == "0.5",]$MIC<-"0.5xMIC"
holder_temp[holder_temp$MIC == "1",]$MIC<-"1xMIC"
holder_temp$MIC<-as.factor(holder_temp$MIC)
holder_temp$MIC<-factor(holder_temp$MIC, levels=c("Control","0.5xMIC","1xMIC"))

ggplot(holder_temp,aes(x=max_protein_load,y=row_id,colour=MIC,fill=MIC))+
  geom_bar(aes(y = ..prop..))+scale_colour_manual(values = c("black","#3182bd","#e6550d"))+scale_fill_manual(values = c("black","#3182bd","#e6550d"))+
  scale_y_continuous(breaks = seq(0,0.7,0.1), limits = c(0,0.7))+
  ylab(expression(paste("Population proportion"))) + xlab(expression(paste("Maximum misfolded protein load")))+
  scale_x_continuous(breaks = seq(0,5,1), limits = c(-0.5,5.5),labels=x_axis)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 30,color = "black",family="Ariel"),
        axis.title = element_text(size = 30,color = "black",family="Ariel"),
        legend.text = element_text(size = 30,color = "black",family="Ariel"),
        legend.title = element_text(size = 30,color = "black",family="Ariel"),
        strip.text = element_text(size= 30,color = "black",family="Ariel"),
        legend.position = "none")+
  facet_wrap( ~ MIC)

########################################################################################################################################################
#chromosome length vs max_protein_load
holder_temp<-holder[holder$treatment != 0,]
holder_temp[holder_temp$max_protein_load > 4,]$max_protein_load<-5
m1<-lmer(log(max_chromosome_length) ~ treatment + (1|uni_id), data = holder_temp)
m2<-lmer(log(max_chromosome_length) ~ max_protein_load*treatment + (1|uni_id), data = holder_temp)
anova(m1,m2)
summary(m2)

f_model<-lmer(log(max_chromosome_length) ~ max_protein_load*treatment+(1|uni_id), data = holder_temp)
mod_sum<-summary(f_model)
r.squaredGLMM(f_model)
holder_temp$y_pred<-NA
holder_temp[holder_temp$treatment == 30,]$y_pred<-mod_sum[["coefficients"]][1,1]+(mod_sum[["coefficients"]][2,1]*holder_temp[holder_temp$treatment == 30,]$max_protein_load)
holder_temp[holder_temp$treatment == 60,]$y_pred<-(mod_sum[["coefficients"]][1,1]+mod_sum[["coefficients"]][3,1])+((mod_sum[["coefficients"]][2,1]+mod_sum[["coefficients"]][4,1])*holder_temp[holder_temp$treatment == 60,]$max_protein_load)

holder_temp$mean_chromosome_length<-NA
m_chr<-tapply(log(holder_temp[holder$treatment == 30,]$max_chromosome_length),holder_temp[holder$treatment == 30,]$max_protein_load,mean)
for ( i in 1:length(m_chr)){
   holder_temp[(holder_temp$max_protein_load == row.names(m_chr)[i]) & (holder_temp$treatment == 30),]$mean_chromosome_length<-m_chr[i]
}
m_chr<-tapply(log(holder_temp[holder$treatment == 60,]$max_chromosome_length),holder_temp[holder$treatment == 60,]$max_protein_load,mean)
for ( i in 1:length(m_chr)){
   holder_temp[(holder_temp$max_protein_load == row.names(m_chr)[i])  & (holder_temp$treatment == 60),]$mean_chromosome_length<-m_chr[i]
}

x_axis<- c("0","1","2","3","4","\u2265 5")

holder_temp[holder_temp$MIC == "0.5",]$MIC<-"0.5xMIC"
holder_temp[holder_temp$MIC == "1",]$MIC<-"1xMIC"

ggplot(holder_temp, aes(x = factor(max_protein_load),y = log(max_chromosome_length), color=MIC))+
  geom_point(data=holder_temp,alpha=0.2,color="black",size=4)+
  geom_point(aes(y=mean_chromosome_length), size = 5)+
  geom_line(data=holder_temp,aes(x=max_protein_load+1,y=y_pred,color=MIC),size=3)+
  scale_colour_manual(values = c("#3182bd","#e6550d"))+
  labs(color = "Ciprofloxacin \nconcentration (ng)")+
  ylab(expression(paste("Maximum nucleoid length (",mu,"m)"))) + xlab(expression(paste("Maximum misfolded protein load")))+
  scale_x_discrete(labels=x_axis)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 30,color = "black",family="Ariel"),
        axis.title = element_text(size = 30,color = "black",family="Ariel"),
        legend.text = element_text(size = 30,color = "black",family="Ariel"),
        legend.title = element_text(size = 30,color = "black",family="Ariel"),
        strip.text = element_text(size= 30,color = "black",family="Ariel"),
        legend.position = "none")+
  facet_wrap(~MIC)+
  scale_y_continuous(breaks=c(1.098612,1.945910,2.995732,4.007333,5.010635,5.991465), labels=c(3,7,20,55,150,400), limits = c(1,6))
  
########################################################################################################################################################
holder_temp<-holder[holder$treatment != 0,]
holder_temp$prot_reduction<-"Yes"
holder_temp[holder_temp$post_max_protein_retention == 1,]$prot_reduction<-"No"
holder_temp$prot_reduction_numeric<-0
holder_temp[holder_temp$prot_reduction == "Yes",]$prot_reduction_numeric<-1

m1<-lmer(log(survival_time) ~ (log(max_chromosome_length):treatment+(1|uni_id)), data = holder_temp , REML = FALSE)
m2<-lmer(log(survival_time) ~ (log(max_chromosome_length):treatment+(1|uni_id)+prot_reduction), data = holder_temp , REML = FALSE)
m3<-lmer(log(survival_time) ~ (log(max_chromosome_length):treatment+(1|uni_id)+prot_reduction:treatment), data = holder_temp , REML = FALSE)

anova(m1,m2)
anova(m2,m3)
summary(m2)
r.squaredGLMM(m2)

final_model<-lmer(log(survival_time) ~ (log(max_chromosome_length):treatment+prot_reduction+(1|uni_id)), data = holder_temp , REML = TRUE)
model_sum<-summary(final_model)
r.squaredGLMM(final_model)

holder_temp$y_pred<-NA
holder_temp[(holder_temp$treatment == 30) & (holder_temp$prot_reduction_numeric == 0),]$y_pred<-model_sum[["coefficients"]][1,1] + 
  model_sum[["coefficients"]][3,1]*log(holder_temp[(holder_temp$treatment == 30) & (holder_temp$prot_reduction_numeric == 0),]$max_chromosome_length)+
  model_sum[["coefficients"]][2,1]*holder_temp[(holder_temp$treatment == 30) & (holder_temp$prot_reduction_numeric == 0),]$prot_reduction_numeric
holder_temp[(holder_temp$treatment == 30) & (holder_temp$prot_reduction_numeric == 1),]$y_pred<-model_sum[["coefficients"]][1,1] + 
  model_sum[["coefficients"]][3,1]*log(holder_temp[(holder_temp$treatment == 30) & (holder_temp$prot_reduction_numeric == 1),]$max_chromosome_length)+
  model_sum[["coefficients"]][2,1]*holder_temp[(holder_temp$treatment == 30) & (holder_temp$prot_reduction_numeric == 1),]$prot_reduction_numeric
holder_temp[(holder_temp$treatment == 60) & (holder_temp$prot_reduction_numeric == 0),]$y_pred<-model_sum[["coefficients"]][1,1] +
  model_sum[["coefficients"]][4,1]*log(holder_temp[(holder_temp$treatment == 60) & (holder_temp$prot_reduction_numeric == 0),]$max_chromosome_length)+
  model_sum[["coefficients"]][2,1]*holder_temp[(holder_temp$treatment == 60) & (holder_temp$prot_reduction_numeric == 0),]$prot_reduction_numeric
holder_temp[(holder_temp$treatment == 60) & (holder_temp$prot_reduction_numeric == 1),]$y_pred<-model_sum[["coefficients"]][1,1] +
  model_sum[["coefficients"]][4,1]*log(holder_temp[(holder_temp$treatment == 60) & (holder_temp$prot_reduction_numeric == 1),]$max_chromosome_length)+
    model_sum[["coefficients"]][2,1]*holder_temp[(holder_temp$treatment == 60) & (holder_temp$prot_reduction_numeric == 1),]$prot_reduction_numeric


1 - sum((log(holder_temp[holder_temp$treatment == 30,]$survival_time) - holder_temp[holder_temp$treatment == 30,]$y_pred)^2)/sum((log(holder_temp[holder_temp$treatment == 30,]$survival_time) - mean(log(holder_temp[holder_temp$treatment == 30,]$survival_time)))^2)
1 - sum((log(holder_temp[holder_temp$treatment == 60,]$survival_time) - holder_temp[holder_temp$treatment == 60,]$y_pred)^2)/sum((log(holder_temp[holder_temp$treatment == 60,]$survival_time) - mean(log(holder_temp[holder_temp$treatment == 60,]$survival_time)))^2)

holder_temp<-within(holder_temp, treat_prot_reduction<-factor(factor(prot_reduction):factor(treatment)))
holder_temp$treat_prot_reduction<-as.character(holder_temp$treat_prot_reduction)


holder_temp$prot_reduction<-as.factor(holder_temp$prot_reduction)
holder_temp$prot_reduction<-relevel(holder_temp$prot_reduction,2)

holder_temp[holder_temp$MIC == "0.5",]$MIC<-"0.5xMIC"
holder_temp[holder_temp$MIC == "1",]$MIC<-"1xMIC"

p<-ggplot(holder_temp, aes(x = log(max_chromosome_length),y = y_pred, color=MIC,alpha=prot_reduction))+
  geom_point(data= holder_temp, aes(x=log(max_chromosome_length),y=log(survival_time), color=MIC,alpha=prot_reduction),size=3)+
  geom_line(size=2)+scale_colour_manual(values = c("#3182bd","#e6550d"))+
  labs(alpha = "Misfolded protein reduction")+
  xlab(expression(paste("Maximum nucleoid length (",mu,"m)"))) + ylab(expression(paste("Survival time (h)")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 30,color = "black",family="Ariel"),
        axis.title = element_text(size = 30,color = "black",family="Ariel"),
        legend.text = element_text(size = 24,color = "black",family="Ariel"),
        legend.title = element_text(size = 24,color = "black",family="Ariel"),
        strip.text = element_text(size= 30,color = "black",family="Ariel"),
        legend.position = c(0.25,0.95),
        legend.background=element_blank(),
        legend.direction = "horizontal")+
  guides(alpha=guide_legend(ncol=1),color="none")+
  facet_wrap( ~ MIC)+
  scale_x_continuous(breaks=c(1.098612,1.945910,2.995732,4.007333,5.010635,5.991465), labels=c(3,7,20,55,150,400), limits = c(1,6))+
  scale_y_continuous(breaks=c(1.098612,1.504077,1.945910,2.484907,2.995732,3.401197), labels=c(3,4.5,7,12,20,30), limits = c(1,3.5))+
  scale_alpha_discrete(range = c(1,0.5))

ggdraw(align_legend(p))
########################################################################################################################################################
holder_temp<-holder[holder$treatment != 0,]
tab1<-table(holder_temp$bud_produced, holder_temp$treatment)
chisq.test(tab1, correct=TRUE)

sum_data = array(data=NA, c(0,5))
for(i in unique(holder_temp$treatment)){
  t<-holder_temp[holder_temp$treatment == i,]
  for (j in unique(t$bud_produced)){
    t1<-t[t$bud_produced == j,]
    c_val<-cbind(i,j,t1$MIC[1], nrow(t1),nrow(t1)/nrow(t))
    sum_data = rbind(sum_data, c_val)
  }
}
sum_data<-as.data.frame(sum_data)
colnames(sum_data)<-c("treatment","bud_produced","MIC","no_cells","prop")
sum_data$no_cells<-as.numeric(sum_data$no_cells)
sum_data$prop<-as.numeric(sum_data$prop)
sum_data$bud_produced<-as.factor(sum_data$bud_produced)
sum_data$bud_produced<-relevel(sum_data$bud_produced,2)
levels(sum_data$bud_produced)<-c("Yes","No")

sum_data[sum_data$MIC == "0.5",]$MIC<-"0.5xMIC"
sum_data[sum_data$MIC == "1",]$MIC<-"1xMIC"

ggplot(sum_data, aes(x=bud_produced,y=prop, color=bud_produced,fill = bud_produced)) + 
  geom_bar(stat="identity") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1))+
  scale_color_manual(values = c("black","black"))+scale_fill_manual(values=c("grey","black"))+
  guides(color = "none")+
  xlab(expression(paste("Bud produced"))) + ylab(expression(paste("Population proportion")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 30,color = "black",family="Ariel"),
        axis.title = element_text(size = 30,color = "black",family="Ariel"),
        legend.text = element_text(size = 30,color = "black",family="Ariel"),
        legend.title = element_text(size = 30,color = "black",family="Ariel"),
        strip.text = element_text(size= 30,color = "black",family="Ariel"),
        legend.position = "none")+
  facet_wrap( ~ MIC)
########################################################################################################################################################
holder_temp<-holder[holder$treatment != 0,]
holder_temp$prot_reduction<-"Yes"
holder_temp[holder_temp$post_max_protein_retention == 1,]$prot_reduction<-"No"
holder_temp$prot_reduction_numeric<-0
holder_temp[holder_temp$prot_reduction == "Yes",]$prot_reduction_numeric<-1

tab1<-table(holder_temp[holder_temp$treatment == 30,]$prot_reduction, holder_temp[holder_temp$treatment == 30,]$bud_produced)
chisq.test(tab1, correct=TRUE)

tab1<-table(holder_temp[holder_temp$treatment == 60,]$prot_reduction, holder_temp[holder_temp$treatment == 60,]$bud_produced)
chisq.test(tab1, correct=TRUE)

prop_temp = array(data=NA, c(0,5))
for (i in unique(holder_temp$treatment)){
  temp<-holder_temp[holder_temp$treatment == i,]
  for (j in unique(holder_temp$prot_reduction)){
    temp1<-temp[temp$prot_reduction == j,]
    for ( k in unique(temp1$bud_produced)){
      prop<-nrow(temp[(temp$prot_reduction == j) & (temp$bud_produced == k),])/nrow(temp1)
      c_tot<-cbind(i,j,temp1$MIC[1],k,prop)
      prop_temp = rbind(prop_temp,c_tot)
    }
  }
}
prop_temp<-as.data.frame(prop_temp)
colnames(prop_temp)<-c("treatment","prot_reduction","MIC","bud_produced","proportion")
prop_temp$proportion<-as.numeric(prop_temp$proportion)
prop_temp$bud_produced<-as.factor(prop_temp$bud_produced)
prop_temp$bud_produced<-relevel(prop_temp$bud_produced,2)
prop_temp$prot_reduction<-as.factor(prop_temp$prot_reduction)
prop_temp$prot_reduction<-relevel(prop_temp$prot_reduction,2)

prop_temp[prop_temp$MIC == "0.5",]$MIC<-"0.5xMIC"
prop_temp[prop_temp$MIC == "1",]$MIC<-"1xMIC"

ggplot(prop_temp, aes(x=bud_produced, y=proportion, color=MIC, fill=MIC,alpha=prot_reduction)) +
  geom_bar(stat="identity", position=position_dodge())+scale_color_manual(values = c("#3182bd","#e6550d"))+scale_fill_manual(values = c("#3182bd","#e6550d"))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1))+
  guides(color = "none")+
  scale_x_discrete(breaks=c("YES","NO"), labels=c("Yes","No"))+
  labs(alpha = "Misfolded protein reduction") + xlab(expression(paste("Bud produced"))) + ylab(expression(paste("Misfolded protein proportion")))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 30,color = "black",family="Ariel"),
        axis.title = element_text(size = 30,color = "black",family="Ariel"),
        legend.text = element_text(size = 24,color = "black",family="Ariel"),
        legend.title = element_text(size = 24,color = "black",family="Ariel"),
        strip.text = element_text(size= 30,color = "black",family="Ariel"),
        legend.position = c(0.20,0.9),
        legend.background=element_blank(),
        legend.direction = "vertical")+
  guides(alpha=guide_legend(ncol=1),fill="none")+
  facet_wrap( ~ MIC)+
  scale_alpha_discrete(range = c(1,0.5))
########################################################################################################################################################



