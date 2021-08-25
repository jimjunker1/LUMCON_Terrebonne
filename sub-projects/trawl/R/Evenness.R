here::i_am("sub-projects/trawl/R/Evenness.R")
#*******************************************************************************
#********************************************************************************
# 
## R scipts "Evenness" for Chao and Ricotta (2019) Ecology paper. 
## This R code is for computing Figures 2, 3 and 4 of Chao and Ricotta (2019) paper.
# NOTE: The packages "ggplot2", "dplyr", "ade4", "reshape2", "ggpubr", "phytools", "ape" must be 
# installed and loaded before running the scripts. 
# 
#
# The following R scripts include two parts:  
# (1). Script for computing the profiles for six classes of evenness measures (see Figure 2 in Chao and Ricotta's paper).
# (2). Script for computing the contribution of each species/node to taxonomic dissimarity and/or phylogenetic
#      dissimarity measures (see Figures 3 and 4 in Chao and Ricotta's paper)
# 
#
#*******************************************************************************
#*******************************************************************************

# library(ggplot2)
# library(dplyr)
# library(ade4)
# library(reshape2)
library(ggpubr)
# library(phytools)
# library(ape)

####################################################################################
#
# (1). Computing the profiles for six classes of evenness measures (Table 1 and Figure 2)
#
####################################################################################

qD <- function(p,q){
  p <- p[p>0]
  if(q!=1){
    (sum(p^q))^(1/(1-q))
  }else{
    exp(-sum(p*log(p)))
  }
}

#' new_fun computes all six classes of evenness measures.
#' @param x is an observed species abundance or frequency vector. 
#' @param q.order is a vector of diversity orders: user must specify a sequence (suggested range is from 0 to 2 in an increment of 0.05).
#' @return the profiles of all six classes of evenness indices listed in Table 1; see Figure 2 for output.
new_fun <- function(x,q.order,evenness.type,...){
  FUN <- qD
  n <- sum(x)
  p <- x/n
  q_profile_evenness <- function(q){
    qDest <- FUN(p,q)
    #S <- sum(x>0)
    S <- sum(x>0)
    E1 <- ifelse(q!=1, (1-qDest^(1-q))/(1-S^(1-q)), log(qDest)/log(S))
    E2 <- ifelse(q!=1, (1-qDest^(q-1))/(1-S^(q-1)), log(qDest)/log(S))
    E3 <- (qDest-1)/(S-1)
    E4 <- (1-1/qDest)/(1-1/S)
    E5 <- log(qDest)/log(S)
    if(q==0){
      p <- p[p>0]
      nu <- abs(p - (1/S))
      nu <- nu[nu > 0]
      sub1 <- (sum(log(abs(nu)))/sum(nu>0)-(log(1-1/S)+(1-S)*log(S))/S)
      E6 <- 1-exp(sub1)
    }else{
      p <- p[p>0]
      E6 <- 1-(sum(abs(p-1/S)^q)/((1-1/S)^q+(S-1)*S^(-q)))^(1/q)
    }
    
    #E6 <- ifelse(q=1, 1-sum(abs(p-1/S)^(1-q))/(abs(1-1/S)^(1-q)+)
    return(c(E1,E2,E3,E4,E5,E6))
  }
  out <- as.matrix(t(sapply(q.order, q_profile_evenness)))
  colnames(out) <- c("E1", "E2", "E3", "E4", "E5", "E6")
  out[,c(evenness.type)]
}


####################################################################################
#
# (2). Computing the contribution of each species/node to dissimilarity measures (Figures 3 and 4)
#      Jaccard-type (1-U_qN) and Sorensen-type (1-C_qN)
#
####################################################################################
#' dis1 computes the contribution of each species/node to the two types of dissimilarity measures.
#' @param x is the species-by-assemblages abundance matrix with species names as rownames.
#' @param q is value for the diversity order.
#' @param type is tax (taxonomic) or phy (phylogenetic).
#' @param type2 is "species" or "k"."species" means the contribution of each species/node to the two types of dissimilarity measures 
#' (Jaccard-type dissimilarity and Sorensen-type dissimilarity).
#' "k" means the contribution of each assemblage/location/site to the two types of dissimilarity measures 
#' (Jaccard-type dissimilarity and Sorensen-type dissimilarity). In the worked example, the contribution of each assemblage/stage is not computed.
#' @tree is the pylog object of the phylogenetic tree of all assemblages.
#' @return the contribution of each species/node to the two types of dissimilarity measures: Jaccard-type (1-U_qN) and Sorensen-type (1-C_qN)
dis1 <- function(x, q, type = "tax", type2 = "species", tree = NULL){
  if(type2 == "species"){
    FUN <- rowSums
  }else{
    FUN <- colSums
  }
  if(type == "tax"){
    x <- as.matrix(x)
    x <- x[rowSums(x)>0, ]
    N <- ncol(x)
    zbar <- rowSums(x)/N
    x1 <- x[zbar>0, ]
    zbar1 <- zbar[zbar>0]
    if(q==0){
      UqN <- FUN(x==0)/((N-1)*(sum(rowSums(x)>0)))
      CqN <- FUN(x==0)/((N-1)*(sum(apply(x, 2, function(i){sum(i>0)}))))
    }else if(q==2){
      UqN <- FUN((x1-zbar1)^2)/((N^q-N)*sum(zbar1^q))
      CqN <- FUN((x1-zbar1)^2)/((1-N^(1-q))*sum(x1^q))
    }else if(q!=1){
      UqN <- FUN((x1)^q-(zbar1)^q)/((N^q-N)*sum(zbar1^q))
      CqN <- FUN((x1)^q-(zbar1)^q)/((1-N^(1-q))*sum(x1^q))
    }else{
      x2 <- x1/zbar1
      UqN <- FUN(x1*log(x2), na.rm = T)/((sum(x)*log(N)))
      CqN <- UqN
    }
  }else{
    Li <- c(tree1$leaves, tree1$nodes)
    cumtree = function(a, tree){
      a <- a[names(tree$leaves)]
      for(i in 1:length(tree$parts)){
        a[1+length(a)] <- sum(a[tree$parts[[i]]])
        names(a)[length(a)] <- names(tree$parts)[i]
      }
      a
    }
    ai <- apply(x, 2, cumtree, tree1)
    wt <- apply(ai, 1, function(x1)(sum(x1))^q/sum(Li*rowSums(ai, na.rm = T)^q))
    N <- ncol(ai)
    zbar <- rowSums(ai)/N
    x1 <- ai[zbar>0, ]
    zbar1 <- zbar[zbar>0]
    Li <- Li[zbar>0]
    T1 <- sum(rowSums(x1)*Li)
    if(q==0){
      if(type2 == "species"){
        rn <- nrow(x1)
        UqN <- sapply(1:rn, function(i){(Li[i]*sum(x1[i, ]==0))})/((N-1)*sum(Li)) 
        CqN <- sapply(1:rn, function(i){(Li[i]*sum(x1[i, ]==0))/((N-1)*sum(Li*rowSums(x1!=0)))})
      }else{
        UqN <- apply(x1, 2, function(x){sum(Li[x==0])})/((N-1)*sum(Li)) 
        CqN <- apply(x1, 2, function(x){sum(Li[x==0])/((N-1)*sum(Li*colSums(x1!=0)))})
      }
      
    }else if(q==2){
      UqN <- FUN(Li*((x1-zbar1)^2), na.rm = T)/((N^q-N)*sum(Li*zbar1^q))
      CqN <- FUN(Li*((x1-zbar1)^2), na.rm = T)/((1-N^(1-q))*sum(Li*x1^q))
    }else if(q!=1){
      UqN <- FUN(Li*((x1)^q-(zbar1)^q), na.rm = T)/((N^q-N)*sum(Li*zbar1^q))
      CqN <- FUN(Li*((x1)^q-(zbar1)^q), na.rm = T)/((1-N^(1-q))*sum(Li*x1^q))
    }else{
      x2 <- x1/zbar1
      UqN <- FUN(Li*x1*log(x2), na.rm = T)/(T1*log(N))
      CqN <- UqN
    }
  }
  
  # c(sum(UqN), sum(CqN))
  rbind(UqN, CqN)
}

#' draw_dis_spe plots the contribution of each species/node to dissimilarity (Jaccard-type dissimilarity and Sorensen-type dissimilarity).
#' @param data is a merged table of output values with three columns corresponding to output for q = 0, 1, 2.
#' @param title_name is the title name of plot. 
#' @type indicates the type of contribution: "tax" for taxonomic and "phy" for phylogenetic       
#' @return the plot of the contribution of each species/node.
draw_dis_spe <- function(data, title_name, type = "tax"){
  # colnames(data) <- c("q = 0", "q = 1", "q = 2")
  # data <- melt(data)
  g <- ggplot(data, aes(x = as.factor(species), y = UqN, fill = q_order))+
    geom_col(width = 0.2)+
    facet_grid(q_order~., scales = "free_y")+
    theme_bw()+
    # ylim(c(0, max(data[, 3])))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3), 
          axis.title = element_text(size = 14), 
          plot.title = element_text(hjust = 0.5))+
    guides(fill=FALSE)+
    ggtitle(title_name)
  
  if(type == "tax"){
    g <- g +
      xlab("Species")+
      ylab("Species contribution")
  }else{
    g <-  g +
      xlab("Species/node")+
      ylab("Species/node contribution")
  }
  return(g)
}



###################################################################################
#
# Example for (2). Alpine species example (See Figures 3 and 4 for output)
#
####################################################################################
######arrange data######
 # data1 <- read.table("Alpine_relative_abundance_data.txt")
# tree <- read.table("Alpine_phylo_tree.txt", header = F)[1,1]
# tree1 <- newick2phylog(tree)
# plot.phylog(tree1,
#             draw.box = TRUE, labels.nodes = names(tree1$nodes), clabel.leaves = 1, clabel.nodes = 1)

# 
# ######caculate the contribution of each species to taxonomic dissimarity and then plot #####
# t01 <- t(dis1(data1, 0, type = "tax", type2 = "species"))
# t11 <- t(dis1(data1, 1, type = "tax", type2 = "species"))
# t21 <- t(dis1(data1, 2, type = "tax", type2 = "species"))
# 
# tax_UqN_r <- cbind(t01[, 1], t11[, 1], t21[, 1])
# tax_CqN_r <- cbind(t01[, 2], t11[, 2], t21[, 2])
# draw_dis_spe(tax_UqN_r, "Jaccard-type taxonomic dissimilarity")
# draw_dis_spe(tax_CqN_r, "Sorensen-type taxonomic dissimilarity")
# 
# ######caculate the contribution of each species/node to phylogenetic dissimarity and then plot #####
# 
# p01 <- t(dis1(data1, 0, type = "phy", type2 = "species", tree = tree1))
# p11 <- t(dis1(data1, 1, type = "phy", type2 = "species", tree = tree1))
# p21 <- t(dis1(data1, 2, type = "phy", type2 = "species", tree = tree1))
# 
# phy_UqN_r <- cbind(p01[, 1], p11[, 1], p21[, 1])
# phy_CqN_r <- cbind(p01[, 2], p11[, 2], p21[, 2])
# 
# draw_dis_spe(phy_UqN_r, "Jaccard-type phylogenetic dissimilarity", type = "phy")
# draw_dis_spe(phy_CqN_r, "Sorensen-type phylogenetic dissimilarity", type = "phy")

##################################################
