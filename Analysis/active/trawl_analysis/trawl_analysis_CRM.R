setwd("~/Dropbox/Terrebonne/LUMCON_Terrebonne")
source("install-packages.R")#installs all necessary packages
source("datascript.R")#imports data with some minor cleanup

library(RColorBrewer)

####Theme Craig####
  theme_craig <- function () { 
    theme_bw(base_size=12) %+replace% 
      theme(
        # change stuff here
        axis.line = element_line(colour = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0))}



# create taxa by site matrix with common names
TB_trawl_commonsite <- TB_trawl_data %>%
  select(year, month, date_id, Common_name) %>%
  group_by(date_id, Common_name) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = Common_name, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()


##addressing outlier points
TB_trawl_outlier_rm <- TB_trawl_commonsite %>%
  filter(date_id %ni% c("2018-1","2007-7"))

#designate the color pallette


varCol = c("#8c510a", "#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#f5f5f5", "#c7eae5", "#c7eae5","#80cdc1","#35978f","#01665e", "#01665e")


TB_pallette <- scale_color_manual(values = c("2007" = "#8c510a",
                                             "2009" = "#8c510a",
                                             "2010" = "#bf812d",
                                             "2011" = "#dfc27d",
                                             "2012" = "#f6e8c3",
                                             "2013" = "#f5f5f5",
                                             "2014" = "#f5f5f5",
                                             "2015" = "#c7eae5",
                                             "2016" = "#c7eae5",
                                             "2017" = "#80cdc1",
                                             "2018" = "#35978f",
                                             "2019" = "#01665e",
                                             "2020" = "#01665e"))
TB_pallette2 <- scale_fill_manual(values = c("2007" = "#8c510a",
                                             "2009" = "#8c510a",
                                             "2010" = "#bf812d",
                                             "2011" = "#dfc27d",
                                             "2012" = "#f6e8c3",
                                             "2013" = "#f5f5f5",
                                             "2014" = "#f5f5f5",
                                             "2015" = "#c7eae5",
                                             "2016" = "#c7eae5",
                                             "2017" = "#80cdc1",
                                             "2018" = "#35978f",
                                             "2019" = "#01665e",
                                             "2020" = "#01665e"))




###Legendre PC with Hellinger

  #Hellinger pre-transformation of the species matrix
    tb.h <- decostand (TB_trawl_outlier_rm %>% select(-c(1:2)) %>% column_to_rownames("date_id"), "hellinger")
    tb.h.pca <- rda(tb.h)
    tb.h.pca.summ <- summary(tb.h.pca)

  
  #pull coordinates and loadings
  TB <- data.frame(matrix(ncol=2,nrow=105, dimnames=list(NULL, c("PC1", "PC2"))))
    
  df1  <- data.frame(tb.h.pca.summ$sites[,1:2])
  TB$PC1 <- df1[,1]
  TB$PC2 <- df1[,2]
  TB$date_id <- row.names(tb.h.pca.summ$sites)
  TB <- left_join(TB, TB_trawl_outlier_rm, by=c("date_id"="date_id"))
  
  Var = factor(TB$year, levels =c("2007","2009","2010","2011",
                                  "2012","2013","2014","2015",
                                  "2016","2017","2018","2019",
                                  "2020")) # factor variable for colours
  
  
  #pull loadings
  TB.loadings  <- data.frame(tb.h.pca.summ$species[,1:2])
  TB.loadings$Species <- row.names(tb.h.pca.summ$species)
  #filtering out top species contributions
  TB.loadings.top <- TB.loadings %>%
    filter(PC1>quantile(PC1, prob=.99) | 
             PC1<quantile(PC1, prob=.01)| 
             PC2>quantile(PC2, prob=.99)  | 
             PC2<quantile(PC2, prob=.01)) 
  #analyses
  tb.h.dl <- dist(tb.h)
  pc_test<- betadisper(tb.h.dl, TB$year)
  anova(pc_test)
  permutest(pc_test, pairwise = TRUE)
  
  TB.centroids <- TB %>%
    group_by(year) %>%
    summarise(
      meanPC1 = mean(PC1),
      meanPC2 = mean(PC2)
    )
  
  TB$year <- as.integer(TB$year)
  

  
  #figure
  ggplot(data = TB, aes(x=PC1, y=PC2, color=as.factor(year))) + 
    geom_point(aes(color=as.factor(year)), cex=3)+
    ggrepel::geom_text_repel(data=TB.loadings.top, aes(x=PC1, y=PC2, label=Species), color="grey30")+
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    TB_pallette+
    TB_pallette2+
    geom_path(data=TB.centroids, aes(x=meanPC1, y=meanPC2), arrow=arrow(), size=1, color="grey20")+
    geom_text(data=TB.centroids, aes(x=meanPC1, y=meanPC2, label=year), color="grey20")+
    theme_craig()+
    labs(color = "Year")
  
###Cluster Analysis
  spe <- TB_trawl_outlier_rm %>% select(-c(1:2)) %>% column_to_rownames("date_id")
  spe.norm <- decostand(spe, "normalize")
  spe.ch <- vegdist(spe.norm, "euc")  
  
  #single linkage
  spe.ch.single <- hclust(spe.ch, method="single")  
  plot(spe.ch.single, cex = 0.6)  
  
  #compelte linkage
  spe.ch.complete <- hclust(spe.ch, method="complete")
  plot(spe.ch.complete, cex = 0.6) 
  
  #UPGMA
  spe.ch.UPGMA <- hclust(spe.ch, method="average")
  plot(spe.ch.UPGMA, cex = 0.6)  
  
  #Centroid
  spe.ch.centroid <- hclust(spe.ch, method="centroid")
  plot(spe.ch.centroid, cex = 0.6)  

  #Ward's
  spe.ch.ward <- hclust(spe.ch, method="ward")
  plot(spe.ch.ward, cex = 0.6)    
  
  #final dendrogram
  spe.chwo <-reorder(spe.ch.ward, spe.ch)
  k=4
  
  plot(spe.chwo, hang=-1, xlab="4 groups", sub="", ylab="Height", 
     main="Chord-Ward (reordered)", labels=TB$year, cex=0.6)

  rect.hclust(spe.chwo, k=k)
  hcd = as.dendrogram(spe.ch.ward)

  require(ggtree)
  
  #trying to plot with ggtree, can't get colors right
  ggtree(hcd)+
    geom_tippoint(cex=3, color=as.factor(TB$year))+
    geom_tiplab(angle=90, hjust=1)+
    layout_dendrogram()+
    TB_pallette+
    TB_pallette2
    
  
# Heat map
# ********
  
      # Heat map of the dissimilarity matrix ordered with the dendrogram
      dend <- as.dendrogram(spe.chwo)
      heatmap(as.matrix(spe.ch), Rowv=dend, symm=TRUE, margin=c(3,3))
      
      # Ordered community table
      # Species are ordered by their weighted averages on site scores
      or <- vegemite(spe, spe.chwo)
      
      # Heat map of the doubly ordered community table, with dendrogram
      heatmap(t(spe[rev(or$species)]), Rowv=NA, Colv=dend,
              col=c("white", brewer.pal(5,"Greens")), scale="none", margin=c(4,4), 
              ylab="Species (weighted averages of sites)", xlab="Sites")
      
      
# k-means partitioning of the pre-transformed species data
# ********************************************************

      # With 4 groups
      spe.kmeans <- kmeans(spe.norm, centers=4, nstart=100)
      
      # Comparison with the 4-group classification derived from Ward clustering:
      spebc.ward.g <- cutree(spe.ch.ward, k)
      table(spe.kmeans$cluster, spebc.ward.g)
      
      # k-means partitioning, 2 to 10 groups
      spe.KM.cascade <- cascadeKM(spe.norm, inf.gr=2, sup.gr=10, iter=100, 
                                  criterion="ssi")
      summary(spe.KM.cascade)
      spe.KM.cascade$results
      spe.KM.cascade$partition
      plot(spe.KM.cascade, sortg=TRUE)
      
      # Reorder the sites according to the k-means result
      spe[order(spe.kmeans$cluster),]
      
      # Reorder sites and species using function vegemite()
      ord.KM <- vegemite(spe, spe.kmeans$cluster)
      spe[ord.KM$sites, ord.KM$species]

#Partitioning around medoids (PAM)
# Computed on the chord distance matrix
# *************************************
      
      require(cluster)
  
      # Choice of the number of clusters
      # Loop to compute average silhouette width for 2 to 28 clusters.
      asw <- numeric(nrow(spe))
      for (k in 2:(nrow(spe)-1)) 
        asw[k] <- pam(spe.ch, k, diss=TRUE)$silinfo$avg.width
      k.best <- which.max(asw) 
      cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
          "with an average silhouette width of", max(asw), "\n")
      dev.new(title="PAM")
      plot(1:nrow(spe), asw, type="h", main="Choice of the number of clusters", 
           xlab="k (number of clusters)", ylab="Average silhouette width")
      axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
           col.axis="red")
      points(k.best, max(asw), pch=16, col="red", cex=1.5)
      
      # PAM for k = 4 clusters
      spe.ch.pam <- pam(spe.ch, k=4, diss=TRUE)
      summary(spe.ch.pam)
      spe.ch.pam.g <- spe.ch.pam$clustering
      spe.ch.pam$silinfo$widths
      
      # Compare with classification from Ward clustering and from k-means
      table(spe.ch.pam.g, spebc.ward.g)
      table(spe.ch.pam.g, spe.kmeans$cluster)
      
      # Silhouette profile for k = 4 groups, k-means and PAM
      dev.new(title="Silhouettes - k-means and PAM", width=12, height=8)
      par(mfrow=c(1,2))
      k <- 4
      sil <- silhouette(spe.kmeans$cluster, spe.ch)
      rownames(sil) <- row.names(spe)
      plot(sil, main="Silhouette plot - k-means", 
           cex.names=0.8, col=2:(k+1))
      plot(silhouette(spe.ch.pam), main="Silhouette plot - PAM", cex.names=0.8, 
           col=2:(k+1))

#Relationships between TB clusters and year
# based on the k-means clustering results (four groups)
# *****************************************************************

    # Boxplots of quantitative environmental variables:
    boxplot(TB$year ~ spe.kmeans$cluster, main="Year", las=1, 
          ylab="Year", col=2:5, varwidth=TRUE)

    # Test of ANOVA assumptions
    # Normality of residuals
    shapiro.test(resid(aov(TB$year ~ as.factor(spe.kmeans$cluster))))
  
    
    # Homogeneity of variances
    bartlett.test(TB$year, as.factor(spe.kmeans$cluster))

    
    # ANOVA of the testable variables
    summary(aov(TB$year ~ as.factor(spe.kmeans$cluster)))

    # Kruskal-Wallis test of variable alt
    kruskal.test(TB$year ~ as.factor(spe.kmeans$cluster))

# Contingency table of two typologies
# ***********************************

  # Year-based typology 
  env2 <- TB$year
  env.de <- vegdist(scale(env2), "euc")
  env.kmeans <- kmeans(env.de, centers=4, nstart=100)
  env.KM.4 <- env.kmeans$cluster
  # Table crossing the k-means and environment 4-group typologies
  table(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster))
  
  # Test the relationship using a chi-square test 
  chisq.test(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster))
  # Change the testing procedure to a permutation test
  chisq.test(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster), 
             simulate.p.value=TRUE)

# Mean abundances on k-means site clusters
# ****************************************
  library(knitr)
  
    groups <- as.factor(spe.kmeans$cluster)
    spe.means <- matrix(0, ncol(spe), length(levels(groups)))
    row.names(spe.means) <- colnames(spe)
    for(i in 1:ncol(spe))
    {
      spe.means[i,] <- tapply(spe[,i], spe.kmeans$cluster, mean)
    }
    # Mean species abundances of the four groups
    group1 <- round(sort(spe.means[,1], decreasing=TRUE), 2)
    group2 <- round(sort(spe.means[,2], decreasing=TRUE), 2)
    group3 <- round(sort(spe.means[,3], decreasing=TRUE), 2)
    group4 <- round(sort(spe.means[,4], decreasing=TRUE), 2)
    # Species with abundances greater than group mean species abundance
    group1.domin <- which(group1 > mean(group1))
    group1
    kable(group1.domin, format="rst")
    write.csv(group1.domin, "group1.csv")
    
    #... same for other groups
        group2.domin <- which(group2 > mean(group2))
        kable(group2.domin, format="rst")
        write.csv(group2.domin, "group2.csv")
        
        group3.domin <- which(group3 > mean(group3))
        kable(group3.domin, format="rst")
        write.csv(group3.domin, "group3.csv")
        
        group4.domin <- which(group4 > mean(group4))
        kable(group4.domin, format="rst")
        write.csv(group4.domin, "group4.csv")
    

        group1.domin <-group1.domin %>%
          as.data.frame() %>%
          rownames_to_column("Species") %>%
          rename ("Rank1" = ".")
        
        group2.domin <-group2.domin %>%
          as.data.frame() %>%
          rownames_to_column("Species") %>%
          rename ("Rank2" = ".")
        
        group3.domin <-group3.domin %>%
          as.data.frame() %>%
          rownames_to_column("Species") %>%
          rename ("Rank3" = ".")
        
        group4.domin <-group4.domin %>%
          as.data.frame() %>%
          rownames_to_column("Species") %>%
          rename ("Rank4" = ".")
        
        dominant_species <- group1.domin %>% 
          full_join(group2.domin, by="Species") %>%
          full_join(group3.domin, by="Species") %>%
          full_join(group4.domin, by="Species")
        
        dominant_species <- dominant_species %>%
          filter(Rank1 <= 10 | Rank2 <= 10 | Rank3 <= 10 | Rank3 <= 10) %>%
          pivot_longer(-Species, names_to = "Group", values_to = "Rank") %>%
          mutate(Group = recode(Group,
                                "Rank1" = "2010",
                                "Rank2" = "2018",
                                "Rank3" = "2013",
                                "Rank4" = "2015"))
          

        library(directlabels)
        ggplot(dominant_species, aes(x=Group, y=Rank, group=Species, color=Species))+
          geom_line()+
          geom_label(aes(label=Rank),label.size = NA, fill = "white")+
          theme_craig()+
          ylab("Year Group")+
          scale_y_reverse()+
          geom_dl(aes(label = Species), method = list(dl.trans(x = x + 0.4), "last.points", cex = 0.8))+
          theme(legend.position = "none")+
          scale_x_discrete(expand=c(0.05, 0, 0.4, 0))+
          theme(axis.line=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank())
          
          
 
  
        

# Fuzzy c-means clustering of the species data
# *************************************************
    
    k <- 4		# Choose the number of clusters
    spe.fuz <- fanny(spe.ch, k=k, memb.exp=1.5)
    summary(spe.fuz)
    
    # Site fuzzy membership
    spe.fuz$membership
    # Nearest crisp clustering
    spe.fuz$clustering
    spefuz.g <- spe.fuz$clustering
    
    # Silhouette plot
    dev.new(title="Fuzzy clustering of fish data - Silhouette plot")
    plot(silhouette(spe.fuz), main="Silhouette plot - Fuzzy clustering", 
         cex.names=0.8, col=spe.fuz$silinfo$widths+1)
    
    # Ordination of fuzzy clusters (PCoA)
    # -----------------------------------
    # Step 1: ordination (PCoA) of the fish chord distance matrix
    dc.pcoa <- cmdscale(spe.ch)
    dc.scores <- scores(dc.pcoa, choices=c(1,2))
    
    # Step 2: ordination plot of fuzzy clustering result
    dev.new(title="Fuzzy clustering of fish data - Ordination plot")
    plot(scores(dc.pcoa), asp=1, type="n",
         main="Ordination of fuzzy clusters (PCoA)")
    abline(h=0, lty="dotted")
    abline(v=0, lty="dotted")
    
    # Step 3: representation of fuzzy clusters
    for (i in 1:k)
    {
      gg <- dc.scores[spefuz.g==i,]
      hpts <- chull(gg)
      hpts <- c(hpts, hpts[1])
      lines(gg[hpts,], col=i+1)
    }
    stars(spe.fuz$membership, location=scores(dc.pcoa), draw.segments=TRUE,
          add=TRUE, scale=FALSE, len=0.1, col.segments=2:(k+1))
    legend(locator(1), paste("Cluster", 1:k, sep=" "),
           pch=15, pt.cex=2, col=2:(k+1), bty="n")
    
