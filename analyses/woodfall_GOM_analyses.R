require(dplyr)
require(ggplot2)
require(vegan)
require(tidyr)
require(stringr)
library(ggrepel)
require(gridExtra)
library(tibble)

################palettes################
study_color_pallette <- scale_color_manual(values = c("Quercus rubra" = "indianred1",
                                                      "Quercus virginiana" = "navajowhite2",
                                                      "Pinus echinata" = "olivedrab3" , 
                                                      "Pinus elliottii" = "olivedrab4",
                                                      "Celtis laevigata" = "palegoldenrod", 
                                                      "1" = "grey30", 
                                                      "2" = "grey60"))
study_color_pallette2 <- scale_fill_manual(values = c("Quercus rubra" = "indianred1",
                                                      "Quercus virginiana" = "navajowhite2",
                                                      "Pinus echinata" = "olivedrab3" , 
                                                      "Pinus elliottii" = "olivedrab4",
                                                      "Celtis laevigata" = "palegoldenrod", 
                                                      "1" = "grey30", 
                                                      "2" = "grey60"))
                                             

################load data################
setwd("~/Dropbox/Return of the Woodfall/New Structure/Wood Falls/Data")
logbytaxa<- data.frame(read.csv("woodfall_logbytaxa.csv", header=TRUE, stringsAsFactors=FALSE))
glimpse(logbytaxa)
head(logbytaxa)

logsummary <-  data.frame(read.csv("woodfall_logs.csv", header=TRUE, stringsAsFactors=FALSE))

NCOL <- ncol(logbytaxa)

#removing stations and species with zero total abundance
logbytaxa2 <- logbytaxa %>%
  replace(is.na(.), 0) %>% #replace all NAs (missing data cells) with zeros
  mutate(rsum = rowSums(.[2:NCOL])) %>% #remove species (rows) with sums of zeros
  filter(rsum>0) %>%
  select(-rsum)

logbytaxa2 <- logbytaxa2[, colSums(logbytaxa2 != 0) > 0] #remove stations (cols) with sums of zero

#time to transpose for vegan
  # first remember the names
  n <- logbytaxa2$ID

  # transpose all but the first column (name)
  logbytaxa3 <- as.data.frame(t(logbytaxa2[,-1]))
  colnames(logbytaxa3) <- n
  logbytaxa3$log_id <- factor(row.names(logbytaxa3))

  str(logbytaxa3) # Check the column types


#################
  

#muthafuckin' diversity time
  logdiversity <- logbytaxa3 %>%
    mutate(Abundance = rowSums(logbytaxa3[,1:141]),                 
           S = specnumber(logbytaxa3[,1:141]),
           H = diversity(logbytaxa3[,1:141],index="shannon"),
           Simp = diversity(logbytaxa3[,1:141],index="simpson"),
           log10Abundance=log10(Abundance)) %>%
    select(Abundance, log10Abundance, S, H, Simp, log_id)
  
logdiversity$log_id<-str_replace_all(logdiversity$log_id, "[.]", "-")
  
logdiversity <- left_join(logdiversity, logsummary, by=c("log_id"="Log.ID"))  

logdiversity$log10Mass <- log10(logdiversity$wood.mass.final)

################

#now some plots

#########create my custon theme#########
theme_craig <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      # change stuff here
      axis.line = element_line(colour = "darkgrey"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position="none",
      strip.background = element_blank(),
      plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
}

p1 <- ggplot(data=logdiversity, aes(y=Abundance, x=log10(wood.mass.final), 
                                    color=wood.type, fill=wood.type))+
  geom_text(aes(label=log_id))+
  theme_craig()+
  study_color_pallette+
  study_color_pallette2+
  xlab("log10 Woodfall Mass (kg)")

p2 <- ggplot(data=logdiversity, aes(y=H, x=log10(wood.mass.final),color=wood.type, fill=wood.type))+
  geom_text(aes(label=log_id))+
  theme_craig()+
  study_color_pallette+
  study_color_pallette2+
  xlab("log10 Woodfall Mass (kg)")

p3 <- ggplot(data=logdiversity, aes(y=S, x=log10(wood.mass.final), 
                                    color=wood.type, fill=wood.type))+
  geom_text(aes(label=log_id))+
  theme_craig()+
  study_color_pallette+
  study_color_pallette2+
  xlab("log10 Woodfall Mass (kg)")

p4 <- ggplot(data=logdiversity, aes(y=Simp, x=log10(wood.mass.final), 
                                    color=wood.type, fill=wood.type))+
  geom_text(aes(label=log_id))+
  theme_craig()+
  study_color_pallette+
  study_color_pallette2+
  xlab("log10 Woodfall Mass (kg)")

grid.arrange(p1, p2, p3,p4, ncol=2)


###################MULTIVARIATE##########################
row.names(logbytaxa3) <- logbytaxa3$log_id
logbytaxa3 <-logbytaxa3 %>% select(-log_id)

#Hellinger pre-transformation of the species matrix
log.h <- decostand (logbytaxa3, "hellinger")
log.h.pca <- rda(log.h )
log.h.pca.summ <- summary(log.h.pca)

#pull coordinates and loadings
df1  <- data.frame(log.h.pca.summ$sites[,1:2])
logdiversity$PC1 <- df1[,1]
logdiversity$PC2 <- df1[,2]

#pull loadings
log.loadings  <- data.frame(log.h.pca.summ$species[,1:2])
log.loadings2 <- log.loadings %>%
  rownames_to_column('Species') %>%
  filter(PC1 > quantile(log.loadings$PC1, 0.95) | PC1 < quantile(log.loadings$PC1, 0.05) |
           PC2 > quantile(log.loadings$PC2, 0.95) | PC2 < quantile(log.loadings$PC2, 0.05) )

#figure
ggplot(data = logdiversity, aes(x=PC1, y=PC2, label=wood.type)) + 
  geom_point(aes(size=log10(wood.mass.final), color = wood.type))+
  geom_text(data=log.loadings2, aes(x=PC1, y=PC2, label=Species), color="grey50")+
  geom_text_repel(aes(label=log_id, color = wood.type),size=4, vjust=-2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_craig()+
  theme(legend.position="none")

##########################comparison to California woodfalls##########################





#callifornia wood falls
    setwd("~/Dropbox/Old Projects/Woodfall/Woodfall_data/")
    wood_orig_cali <- read.csv("woodfall_species.csv", header=TRUE) #reads in data
    glimpse(wood_orig_cali)
    
    #separate out just taxa
    wood_taxa_cali <- wood_orig_cali %>% 
      select(Xylophaga.zierenbergi:Limpet2)
    row.names(wood_taxa_cali) <- wood_orig_cali$Log
    
    #create a main log data holder
    wood_main_cali <- wood_orig_cali %>% 
      select(Log:Size.Order)
    
    wood_main_cali$log10mass <- log10(wood_main_cali$Weight..kg.)
    wood_main_cali$Set<-as.factor(wood_main_cali$Set)
    
    richness = data.frame(specnumber(wood_taxa_cali)); colnames(richness) = c("richness") 
    wood_main_cali = cbind(wood_main_cali, richness)
    
    ########plot of GoM and Calinfornia woodfalls#####
    ggplot(data=logdiversity, aes(y=S, x=log10Mass))+
      geom_point(cex=3, aes(color=wood.type, fill=wood.type))+
      geom_smooth(method=lm, se=FALSE, data=logdiversity, 
                  aes(y=S, x=log10Mass, group=wood.type, color=wood.type), linetype = "dotted")+
      geom_point(cex=3, data=wood_main_cali, aes(x=log10mass, y=richness, color=Set, fill=Set))+
      geom_smooth(method=lm, se=FALSE, data=wood_main_cali, 
                  aes(x=log10mass, y=richness, group=Set, color=Set), linetype = "dotted")+ 
      study_color_pallette+
      study_color_pallette2+
      theme_craig()+
      xlab("log10 Woodfall Mass (kg)")
##############plot of wood sizes#############
    logsummary <- logsummary %>%
      filter(experiment!="shallow") %>%
      mutate(log10Mass = log10(wood.mass.final))
    
    
    ggplot(data=logsummary, aes(x=log10Mass, y=Status, color=wood.type, fill=wood.type))+
      geom_point()+
      geom_text_repel(aes(label=tag.number))+
      facet_wrap(~experiment)+
      theme_craig()
    
    
    
###########ecosystem function###########
    
    logbytaxa_xylo <- logbytaxa3 %>% 
      select(contains("xylo"))
    
    logdiversity <- logdiversity %>%
      mutate(xylo_Abundance = rowSums(logbytaxa_xylo),                 
             xylo_S = specnumber(logbytaxa_xylo),
             xylo_H = diversity(logbytaxa_xylo,index="shannon"),
             xylo_Simp = diversity(logbytaxa_xylo,index="simpson"),
             xylo_log10Abundance=log10(xylo_Abundance)) 
    
    logdiversity$log10consumed <- log10(logdiversity$wood.mass.final-logdiversity$weight_final)
    logdiversity$consumed <- logdiversity$wood.mass.final-logdiversity$weight_final

    xylo_A <- ggplot(data=logdiversity, aes(x=xylo_Abundance, y=log10(consumed)))+
      geom_point(cex=4, aes(color=wood.type, fill=wood.type))+
      geom_smooth(method=lm, se=FALSE, color="grey70")+
      geom_text_repel(aes(label=log_id, color=wood.type))+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")
    
    xylo_S <- ggplot(data=logdiversity, aes(x=xylo_S, y=log10(consumed)))+
      geom_point(cex=4, aes(color=wood.type, fill=wood.type))+
      geom_smooth(method=lm, se=FALSE, color="grey70")+
      geom_text_repel(aes(label=log_id, color=wood.type))+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("S")
    
    xylo_H <- ggplot(data=logdiversity, aes(x=xylo_H,log10(consumed)))+
      geom_point(cex=4, aes(color=wood.type, fill=wood.type))+
      geom_smooth(method=lm, se=FALSE, color="grey70")+
      geom_text_repel(aes(label=log_id, color=wood.type))+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("H")
    
    xylo_J <- ggplot(data=logdiversity, aes(x=xylo_Simp,y=log10(consumed)))+
      geom_point(cex=4, aes(color=wood.type, fill=wood.type))+
      geom_smooth(method=lm, se=FALSE, color="grey70")+
      geom_text_repel(aes(label=log_id, color=wood.type))+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Simpson's J")
    
    grid.arrange(xylo_A, xylo_S, xylo_H, xylo_J)

    
    ####xylo plots####
    xylo1<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-1` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-1")
    
    
    xylo2<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-2` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-2")
    
    
    xylo3<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-3` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-3")
    
    
    xylo4<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-4` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-4")
    
    
    xylo5<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-5` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-5")
    
    
    xylo6<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-6` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-6")
    
    
    xylo7<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-7` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-7")
    
    
    xylo8<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-8` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-8")
    
    
    xylo9<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-9` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-9")
    
    
    xylo10<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-10` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-10")
    
    
    xylo11<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-11` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-11")
    
    
    xylo12<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-12` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-12")
    
    
    xylounid<- ggplot(data=logdiversity, aes(x=logbytaxa_xylo$`xylo-unid` , consumed, color=wood.type, fill=wood.type ))+
      geom_point(cex=4)+
      ggtitle("All Woodfalls")+
      theme_craig()+
      study_color_pallette+
      study_color_pallette2+
      ylab("Consumed Weight (kg)")+
      xlab("Abundance")+
      ggtitle("Xylo-Unid")
    
    grid.arrange(xylo1,xylo2,xylo3,xylo4,xylo7, xylo10,xylo11,xylo12,xylounid)  
    

#################################covariance/variance######################
    require(codyn)
    require(knitr)
    
    #create dataframe that works for species synchrony
      #create triplet from logbytaxa2
      wf_triplet <-pivot_longer(data=logbytaxa2, cols=L.Y.31:L.W.22, names_to = "woodfall", values_to = "count")
      #need to add in weights
      #damn log names don't match periods and dashes, seriously!
      wf_triplet$woodfall <- gsub(".","-",wf_triplet$woodfall,fixed=TRUE)
      #do a left join to get in weights and then keep only columns I want
      wf_triplet <- left_join(wf_triplet, logdiversity, by= c("woodfall" = "log_id")) %>%
        select(ID, woodfall, count, log10Mass, category) 
        
      
      wf_triplet$sizeclass <-cut(wf_triplet$log10Mass, breaks=c(-Inf, 0.5 , Inf), labels=c("small","big"))
      wf_triplet_small <- wf_triplet %>% filter (sizeclass=="small")
      wf_triplet_big <- wf_triplet %>% filter (sizeclass=="big")
    
      wf_synchrony_Gross <- synchrony(df = wf_triplet, 
                                     time.var = "log10Mass", 
                                     species.var = "ID",  
                                     abundance.var = "count", 
                                     metric = "Gross")
      
      wf_synchrony_Gross_sm <- synchrony(df = wf_triplet_small, 
                                      time.var = "log10Mass", 
                                      species.var = "ID",  
                                      abundance.var = "count", 
                                      metric = "Gross")
      
      wf_synchrony_Gross_lg <- synchrony(df = wf_triplet_big, 
                                      time.var = "log10Mass", 
                                      species.var = "ID",  
                                      abundance.var = "count", 
                                      metric = "Gross")
      wf_synchrony_Gross
      wf_synchrony_Gross_sm
      wf_synchrony_Gross_lg 
      

