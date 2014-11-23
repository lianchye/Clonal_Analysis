# Model based Clustering: Finite Mixture Model of Multinomial for Single Cell Analysis -----

# Required Libraries----
require(flexmix)          # For model bases clustering
require(RColorBrewer)     # For color palettes in plotting
require(ggplot2)          # For plotting
require(FactoMineR)       # For performing MCA
require(vegan)            # For calculating jaccard distance
require(gplots)           # For plotting heatmaps
require(Heatplus)         # For plotting heatmaps with annotations
require(reshape)          # For melting dataframe
require(adegenet)         # For performing the Directed Minimum Spanning tree for clones
require(igraph)           # For plotting of trees 

# Clear the workspace-------
rm(list=ls())

# Get all the genotype data with the cluster annotation -----
pt1.cluster.dir <- "./pat_1_cluster_dat"
pt2.cluster.dir <- "./pat_2_cluster_dat"
pt3.cluster.dir <- "./pat_3_cluster_dat"
pt4.cluster.dir <- "./pat_4_cluster_dat"
pt5.cluster.dir <- "./pat_5_cluster_dat"
pt6.cluster.dir <- "./pat_6_cluster_dat"

pt1_cluster <- read.table(pt1.cluster.dir,header=TRUE)
pt2_cluster <- read.table(pt2.cluster.dir,header=TRUE)
pt3_cluster <- read.table(pt3.cluster.dir,header=TRUE)
pt4_cluster <- read.table(pt4.cluster.dir,header=TRUE)
pt5_cluster <- read.table(pt5.cluster.dir,header=TRUE)
pt6_cluster <- read.table(pt6.cluster.dir,header=TRUE)

#####################----
# Declare Functions #
#####################----

# Create the dataframe for input----
my_create_in_df <- function(pt1_cluster){
  bin.mat <- as.matrix(pt1_cluster[,-which(colnames(pt1_cluster)=="clusters")])
  test.df <- data.frame(Freq=as.integer(rep(1,dim(bin.mat)[1])))
  test.df$Incidence <- bin.mat
  return(test.df)
}

# Perform MDS on the binary data----
my_MDS          <- function(pt1_cluster, pt1_fmm_best){
  MDS_in      <- data.frame(apply(pt1_cluster[,-which(colnames(pt1_cluster)=="clusters")],2,function(x){factor(x)}))
  res.mca     <- MCA(MDS_in, graph = FALSE)
  MDS_vars_df <- data.frame(res.mca$var$coord)
  MDS_inds_df <- data.frame(res.mca$ind$coord)
  MDS_inds_df$cluster_assign <- factor(flexmix::clusters(pt1_fmm_best))
  return(list(MDS_vars_df,MDS_inds_df))
}

# Perform the information criterion for cluster analysis----
my_ic_fun       <- function(pt_fmm){
  bic.df <- data.frame(nos_clus=seq(1:length(BIC(pt_fmm))),IC=BIC(pt_fmm))
  aic.df <- data.frame(nos_clus=seq(1:length(AIC(pt_fmm))),IC=AIC(pt_fmm))
  bic.df$measure <- rep("bic",dim(bic.df)[1])
  aic.df$measure <- rep("aic",dim(aic.df)[1])
  ic.df  <- rbind(bic.df,aic.df)
  g <- ggplot(ic.df,aes(x=nos_clus,y=IC,group=measure,color=measure))+
       geom_point(data=subset(ic.df,IC%in%c(min(bic.df$IC),min(aic.df$IC))),color="black",size=7.5,alpha=0.2)+
       geom_point(size=3)+
       geom_line(alpha=0.7)+
       scale_colour_discrete(name  ="Measure",
                            breaks=c("aic", "bic"),
                            labels=c("Akaike", "Bayesian"))+
       ggtitle("Number of Clones Selection Using Bayesian/Akaike Information Criterion")+
       xlab("Number of clusters")+
       ylab("Information Criterion")+
       theme_bw()
  g
  return(g)
}

# Generate the parameters from the FMM model---- 
my_params_df    <- function(pt_fmm_best , threshold=0.6){
  pt_params      <- data.frame(parameters(pt_fmm_best))
  pt_params$cat  <- (apply(pt_params,1,function(x){sum((x>threshold)+0)}))
  pt_params$snps <- row.names(pt_params)
  levels(pt_params$snps)  <- pt_params$snps[order(-pt_params$cat)]
  pt_params_melt          <- melt.data.frame(pt_params,id=c("snps","cat"))
  pt_params_melt$variable <- factor(gsub("Comp.","Clone ", pt_params_melt$variable))
  pt_params_melt$cat[pt_params_melt$value<threshold] <- "Below Thres"
  pt_params_melt$cat      <- as.factor(pt_params_melt$cat)
  return(pt_params_melt)
}

# Generate the clonal consensus genotype from the FMM parameters
my_param_2_cl   <- function(pt_fmm_best,bin_thres=0.65){
  pt_cl_geno           <- data.frame (t(parameters(pt_fmm_best)))
  pt_cl_geno[pt_cl_geno >= bin_thres] <- 1
  pt_cl_geno[pt_cl_geno <  bin_thres] <- 0
  rownames  (pt_cl_geno) <- gsub     ("Comp."   , "Clone_" ,rownames(pt_cl_geno))
  colnames  (pt_cl_geno) <- gsub     ("center." , ""       ,colnames(pt_cl_geno))
  return    (pt_cl_geno)
}

# Add primodial genotype as an additional clone
my_add_pri      <- function(pt_cl_geno){
  clone_0       <- (colSums(pt_cl_geno)==dim(pt_cl_geno)[1])+0
  clone_out     <- rbind(clone_0,pt_cl_geno)
  row.names(clone_out)[1] <- "clone_0"
  return(clone_out)
}

# Generate the Clonal timepoints based on mutational frequency
# clus_size     <- table     (flexmix::clusters(pt1_fmm_best))
my_pre_tree     <- function(pt_cl_geno,clus_size,primary=TRUE){
  mut_freq      <- apply     (pt_cl_geno,2,function(x){
                   return((sum(clus_size*as.vector(x))/sum(clus_size)))})
  log_cl_freq   <- apply     (pt_cl_geno,1,function(x){
                   return(ceiling(-log10(prod(mut_freq[which(x==1)]))))})
  clone_time    <- c(Sys.Date()+log_cl_freq)
  
  if(primary)   {pt_cl_geno           <- my_add_pri(pt_cl_geno)
                 clone_time           <- c(Sys.Date(),Sys.Date()+log_cl_freq)
                 names(clone_time)[1] <- "Clone_0"}
  out           <- list(pt_cl_geno,clone_time)
  return(out)
}

# Perform the tree generation algorithm
my_run_tree     <- function(clone_gen_dis,clone_time){
  clone_gen_dis <- as.matrix (vegdist(clone_gen_dis,method="jaccard"))
  sqtk.res.add  <- seqTrack  (clone_gen_dis, 
                              x.names = row.names(clone_gen_dis), 
                              x.dates = clone_time)
  g_pri         <- plot(sqtk.res.add,vertex.size=4)
  return(list(g_pri,sqtk.res.add))
}

# Tree plotting layout function ----
my_layout<-function(g_pri,root_in=1){
  layout.mat  <- layout.reingold.tilford(g_pri,root=root_in)     # Get the layout matrix using the general reingold method
  scaling_fac <- E(g_pri)$weight/E(g_pri)$weight[1] # Generate the scaling factor for each of the factor
  layout.el   <- get.edgelist(g_pri)                # Get the edge list matrix
  edges_tot   <- length(layout.el[,1])              # Get the total number of edges
  for (i in seq(1,edges_tot)){
    message(paste0("Scaling for edge: ",i))
    xy_parent_node  <- layout.mat[as.numeric(layout.el[i,1]),]
    xy_daught_node  <- layout.mat[as.numeric(layout.el[i,2]),]
    update_x        <- xy_daught_node[1]
    update_y        <- xy_parent_node[2]-scaling_fac[i]
    layout.mat[as.numeric(layout.el[i,2]),1] <- update_x
    layout.mat[as.numeric(layout.el[i,2]),2] <- update_y
  }
  return(layout.mat)
}

# Generate the heatmaps----
# Returns: Heatmap[[1]] -> EM clustering heatmap
# Returns: Heatmap[[2]] -> EM contrast with hclust
# Returns: vector       -> hclust results
# ---------------------
my_heatmap      <- function(pt_cluster , pt_fmm_best , nos_clust,nos_mut){
  bin.mat      <- as.matrix(pt_cluster[,-which(colnames(pt_cluster)=="clusters")])
  EM_cluster   <- factor(flexmix::clusters           (pt_fmm_best))
  nos_cluster  <- max   (as.numeric(flexmix::clusters(pt_fmm_best)))
  row_order    <- order (as.numeric(flexmix::clusters(pt_fmm_best)))
  data.dist    <- vegdist  (bin.mat     , method = "jaccard")
  data.dist.g  <- vegdist  (t(bin.mat)  , method = "jaccard")
  row.clus     <- hclust   (data.dist   , "ward.D2")
  col.clus     <- hclust   (data.dist.g , "ward.D2")
  hclust_ass   <- cutree   (row.clus    , nos_clust)
  color_scheme <- colorRampPalette(c("white", "#660000"), space = "rgb")(2)
  p            <- annHeatmap2(bin.mat[row_order,],
                  scale  = "none", 
                  col    = color_scheme, breaks = 2,
                  legend = 3,
                  dendrogram = list(Col = list(dendro = as.dendrogram(col.clus)),
                                    Row = list(status = "no"                   ) ),
                  cluster    = list(Col = list(cuth   = col.clus$height[length(col.clus$height)-nos_mut+1]     ),
                                    Row = list(grp    = EM_cluster[row_order],
                                               col    = brewer.pal(nos_cluster,"Set2")[seq(1,nos_cluster,by=1)]) ),
                  ann        = list(Row = list(data   = data.frame(EM_cluster=EM_cluster[row_order])))
                  )
  p_row        <- annHeatmap2(bin.mat,
                  scale  = "none", 
                  col    = color_scheme, breaks = 2,
                  legend = 3,
                  dendrogram = list(Col = list(dendro = as.dendrogram(col.clus)),
                                    Row = list(dendro = as.dendrogram(row.clus)) ),
                  cluster    = list(Col = list(cuth   = col.clus$height[length(col.clus$height)-nos_mut+1]    ),
                                    Row = list(cuth   = row.clus$height[length(row.clus$height)-nos_clust+1]) ),
                  ann        = list(Row = list(data   = data.frame(EM_cluster=EM_cluster)))
                  )
   return(list(p,p_row,hclust_ass))
}

# Generates P.values of mutations in each clone given estimated Allele dropout Rate ------
# --------------------
my_pvalue_mut   <- function(pt_cluster , pt_fmm_best, ADO){
  bin.df        <- pt_cluster[,-which(colnames(pt_cluster)=="clusters")]
  pvalue.df     <- do.call("rbind", as.list(
                   by(bin.df, list(clusters=paste0("Clone_",flexmix::clusters(pt_fmm_best))),
                      function (x){
                         clu_size <- dim(x)[1]
                         p.value  <- apply(x,2,function(y){
                         return((binom.test(sum(y), clu_size, p=1-ADO, alternative = "greater"))$p.value)})
                      return   (p.value)}
                   )))
  return(pvalue.df)
}

# Generates P.values of mutations between clone ------
# --------------------
my_pvalue_clone <- function(pt_cluster , pt_fmm_best){
  bin.df        <- pt_cluster[,-which (colnames(pt_cluster)=="clusters")]
  bin.df$clus   <- factor    (flexmix::clusters(pt_fmm_best))
  nos_clus      <- length    (levels  (bin.df$clus))
  output.list   <- list()
  output.vec    <- vector    ()
  for(i in levels(bin.df$clus)){
     for(j in levels(bin.df$clus)){
        clus_df1      <- subset    (bin.df  ,clus==i)
        clus_df2      <- subset    (bin.df  ,clus==j)
        clus_df1_cs   <- colSums   (clus_df1  [,-which(colnames(clus_df1)=="clus")])
        clus_df2_cs   <- colSums   (clus_df2  [,-which(colnames(clus_df2)=="clus")])
        clus_df1.size <- dim       (clus_df1) [1]
        clus_df2.size <- dim       (clus_df2) [1]
        clus_prop.df  <- as.matrix (cbind (clus_df1_cs   , clus_df2_cs  ))
        clus_in.df    <- clus_prop.df [clus_prop.df[,1] > 0 & 
                                       clus_prop.df[,2] > 0 ,]
        prob          <- rep       (clus_df1.size/(clus_df1.size + clus_df2.size),
                                    dim(clus_in.df)[1])
        test.out      <- prop.test (clus_in.df, 
                                    alternative = "two.sided",
                                    p           = prob , 
                                    conf.level  = 0.95 ,
                                    correct     = TRUE)
        output.vec    <- rbind     (output.vec,test.out$p.value)
     }
  }
  output        <- matrix(output.vec , nrow = nos_clus)
  return(output)
}

# Generates P.values of clone specifying mutations -----
# -----
my_pvalue_cl_mu <- function(pt_cluster , pt_fmm_best){
  bin.df        <- pt_cluster[,-which (colnames(pt_cluster)=="clusters")]
  bin.df$clus   <- factor    (flexmix::clusters(pt_fmm_best))
  nos_clus      <- length    (levels  (bin.df$clus))
  output.vec    <- vector    ()
  for(i in levels(bin.df$clus)){
    for(j in levels(bin.df$clus)){
      clus_df1      <- subset    (bin.df  ,clus==i)
      clus_df2      <- subset    (bin.df  ,clus==j)
      clus_df1_cs   <- colSums   (clus_df1  [,-which(colnames(clus_df1)=="clus")])
      clus_df2_cs   <- colSums   (clus_df2  [,-which(colnames(clus_df2)=="clus")])
      clus_df1.size <- dim       (clus_df1) [1]
      clus_df2.size <- dim       (clus_df2) [1]
      clus_prop.df  <- as.matrix (cbind (clus_df1_cs   , clus_df2_cs  ))
      clus_in.df    <- clus_prop.df [clus_prop.df[,1] > 0 & 
                                     clus_prop.df[,2] > 0 ,]
      prob          <- rep       (clus_df1.size/(clus_df1.size + clus_df2.size),
                                  dim(clus_in.df)[1])
      test.out      <- prop.test (clus_in.df, 
                                  alternative = "two.sided",
                                  p           = prob , 
                                  conf.level  = 0.95 ,
                                  correct     = TRUE)
      output.vec    <- rbind     (output.vec,test.out$p.value)
    }
  }
  output        <- matrix(output.vec , nrow = nos_clus)
  return(output)
}

# Generates a random MONOCLONAL cell dataset -----
# -----
my_rand_mon_bin <- function(nos_cell=20, nos_mut=10, ADO=0.3) {
  random.vec    <- vector()
  counter_mut   <- 0
  repeat{
    per_mut     <- rbinom(nos_cell,1,1-ADO)
    random.vec  <- c(random.vec,per_mut)
    counter_mut <- counter_mut + 1
    if(counter_mut == nos_mut){break}
  }
  out.df           <- data.frame(matrix(random.vec,nrow=nos_cell))
  colnames(out.df) <- paste0("mut_",seq(1,nos_mut,by=1))
  out.df$clusters  <- as.integer(rep(1,nos_cell))
  return(out.df)
}

# Generates a random clonal Genotype -----
# -----
my_rand_clo_bin <- function(nos_clone=4, nos_mut=10 ,prob=0.5) {
  random.vec    <- vector()
  counter_mut   <- 0
  repeat{
    per_mut     <- rbinom(nos_clone,1,prob)
    random.vec  <- c(random.vec,per_mut)
    counter_mut <- counter_mut + 1
    if(counter_mut == nos_mut){break}
  }
  out.df           <- data.frame(matrix(random.vec,nrow=nos_clone))
  colnames(out.df) <- paste0("mut_",seq(1,nos_mut,by=1))
  return(out.df)
}

# Generates random cell data from clonal genotype

my_rand_clu     <- function(clonal_genotype,nos_cell_vec,ADO,noise){
  out_all.df    <- NULL
  for (i in seq(1,dim(clonal_genotype)[1],by=1)){
      random.vec       <- vector()
      for (j in seq(1,dim(clonal_genotype)[2],by=1)){
          if(clonal_genotype[i,j]==1){
            per_mut    <- rbinom(nos_cell_vec[i], 1 , 1-ADO)
            random.vec <- c(random.vec,per_mut)
          }else{
            per_mut    <- rbinom(nos_cell_vec[i], 1 , noise)
            random.vec <- c(random.vec,per_mut)
          }
      }
      out.df           <- data.frame (matrix(random.vec,nrow=nos_cell_vec[i]))
      colnames(out.df) <- paste0     ("mut_",seq(1,dim(clonal_genotype)[2],by=1))
      out.df$clusters  <- as.integer (rep   (i,nos_cell_vec[i]))  
      out_all.df       <- rbind      (out_all.df,out.df)
  }
  return(out_all.df)
}

# Runs the random function to get the clustering simulation----
# ----
my_rand_iter    <- function(iter_clones=5, iter_mutations=60, iter_ADO=0.4, repli=3){
  out.df <- NULL
  for (i in seq(1    , iter_clones    , by=1   )) {
  for (j in seq(10   , iter_mutations , by=10   )) {
  for (k in seq(0.2  , iter_ADO       , by=0.05)) {
  count_rep=1
  repeat{
        message(paste0("clone: ",i," Mut: ",j," ADO: ",k," Repli: ",count_rep))
        pti_random_cl   <- my_rand_clo_bin (nos_clone    = i  ,
                                            nos_mut      = j  ,
                                            prob         = 0.5)       # Generate random clone genotypes
        pti_cluster     <- my_rand_clu     (pti_random_cl                 ,
                                            nos_cell_vec = sample(1:50,i) ,
                                            ADO          = k              ,
                                            noise        = 0.00001)   # Generate random patient data
        pti_mb_clus.df  <- my_create_in_df (pti_cluster)              # Creates the input for EM
        pti_fmm         <- stepFlexmix     (Incidence    ~ 1              ,
                                            weights      = ~ Freq         , 
                                            data         = pti_mb_clus.df ,
                                            model        = FLXMCmvbinary(truncated = TRUE ),
                                            control      = list         (minprior  = 0.005), 
                                            k            = 1:7, 
                                            nrep         = 5)         # Perform the model base clustering
        pti_ic_plot     <- my_ic_fun       (pti_fmm)                  # Generate Plots using Information Criterion
        bic_k           <- dim   (posterior(getModel(pti_fmm,"BIC")))[2]
        aic_k           <- dim   (posterior(getModel(pti_fmm,"BIC")))[2]
        out_vec         <- c     (i,j,k,count_rep,bic_k,aic_k)
        message(paste0("sim_clone: ",out_vec[1]," sim_mut: ",out_vec[2]," sim_ADO: ",out_vec[3],
                       " Repli: "   ,out_vec[4]," bic_k: "  ,out_vec[5]," aic_k: "  ,out_vec[6]))
        out.df          <- rbind (out.df,out_vec)
        count_rep       <- count_rep + 1
        if(count_rep == repli+1){break}
  }
  }
  }
  }
  colnames(out.df)<-c("sim_clone","sim_mut","sim_ADO","Repli","bic_k","aic_k")
  return(out.df)
}

#Sensitivity simulation results
#-----

my_sens_iter    <- function(iter_sens=c(0.01,0.02,0.05,0.10,0.25), iter_mutations=60, iter_ADO=0.4, repli=3, num_cells=100){
  out.df <- NULL
  for (i in iter_sens) {
  for (j in seq(60   , iter_mutations , by=10   )) {
  for (k in seq(0.2 , iter_ADO       , by=0.05)) {
        count_rep=1
        repeat{
          message(paste0("clone: ",i," Mut: ",j," ADO: ",k," Repli: ",count_rep))
          pti_random_cl   <- my_rand_clo_bin (nos_clone    = 2  ,
                                              nos_mut      = j  ,
                                              prob         = 0.5)       # Generate random clone genotypes
          pti_cluster     <- my_rand_clu     (pti_random_cl                 ,
                                              nos_cell_vec = ceiling(c(i,1-i)*num_cells) ,
                                              ADO          = k              ,
                                              noise        = 0.00001)   # Generate random patient data
          pti_mb_clus.df  <- my_create_in_df (pti_cluster)              # Creates the input for EM
          pti_fmm         <- stepFlexmix     (Incidence    ~ 1              ,
                                              weights      = ~ Freq         , 
                                              data         = pti_mb_clus.df ,
                                              model        = FLXMCmvbinary(truncated = TRUE ),
                                              control      = list         (minprior  = 0.001), 
                                              k            = 1:4, 
                                              nrep         = 5)         # Perform the model base clustering
          pti_ic_plot     <- my_ic_fun       (pti_fmm)                  # Generate Plots using Information Criterion
          bic_k           <- dim   (posterior(getModel(pti_fmm,"BIC")))[2]
          aic_k           <- dim   (posterior(getModel(pti_fmm,"BIC")))[2]
          out_vec         <- c     (i,j,k,count_rep,bic_k,aic_k)
          message(paste0(" sensitivity: ",out_vec[1]," sim_mut: ",out_vec[2]," sim_ADO: ",out_vec[3],
                         " Repli: "      ,out_vec[4]," bic_k: "  ,out_vec[5]," aic_k: "  ,out_vec[6]))
          out.df          <- rbind (out.df,out_vec)
          count_rep       <- count_rep + 1
          if(count_rep == repli+1){break}
        }
      }
    }
  }
  colnames(out.df)<-c("sensitivity","sim_mut","sim_ADO","Repli","bic_k","aic_k")
  return(out.df)
}

sens_simu_data<-my_sens_iter(iter_sens=c(seq(0.01,0.10,by=0.01),seq(0.2,0.5,by=0.1)), iter_mutations=50, iter_ADO=0.35, repli=3, num_cells=100)
save(sens_simu_data,file="./sens_simu_data.rdata_100")

for(nos_cell_sim in seq(25,75,by=25)){
  assign      (paste0("sens_simu_data_num_cell_",nos_cell_sim),
               my_sens_iter(iter_sens      = c(seq(0.01,0.10,by=0.01),seq(0.2,0.5,by=0.1)),
                            iter_mutations = 50, 
                            iter_ADO       = 0.35,
                            repli          = 3,
                            num_cells      = nos_cell_sim))
  write.out <- get(paste0("sens_simu_data_num_cell_",nos_cell_sim))
  write.csv   (write.out,file=paste0("./sens_simu_data.rdata_",nos_cell_sim))
  save        (write.out,file=paste0("./sens_simu_data.rdata_",nos_cell_sim))
}

for(nos_cell_sim in seq(125,200,by=25)){
  assign      (paste0("sens_simu_data_num_cell_",nos_cell_sim),
               my_sens_iter(iter_sens      = c(seq(0.01,0.10,by=0.01),seq(0.2,0.5,by=0.1)),
                            iter_mutations = 50, 
                            iter_ADO       = 0.35,
                            repli          = 3,
                            num_cells      = nos_cell_sim))
  write.out <- get(paste0("sens_simu_data_num_cell_",nos_cell_sim))
  write.csv   (write.out,file=paste0("./sens_simu_data.csv.rdata_" ,nos_cell_sim))
  save        (write.out,file=paste0("./sens_simu_data.rdata_"     ,nos_cell_sim))
}

for(nos_cell_sim in seq(225,300,by=25)){
  assign      (paste0("sens_simu_data_num_cell_",nos_cell_sim),
               my_sens_iter(iter_sens      = seq(0.01,0.10,by=0.01),
                            iter_mutations = 100, 
                            iter_ADO       = 0.2,
                            repli          = 3,
                            num_cells      = nos_cell_sim))
  write.out <- get(paste0("sens_simu_data_num_cell_",nos_cell_sim))
  write.csv   (write.out,file=paste0("./sens_simu_data.csv.rdata_" ,nos_cell_sim))
  save        (write.out,file=paste0("./sens_simu_data.rdata_"     ,nos_cell_sim))
}

for(nos_cell_sim in seq(25,200,by=25)){
  assign      (paste0("sens_simu_data_num_cell_",nos_cell_sim,"_2"),
               my_sens_iter(iter_sens      = seq(0.01,0.10,by=0.01),
                            iter_mutations = 100, 
                            iter_ADO       = 0.2,
                            repli          = 3,
                            num_cells      = nos_cell_sim))
  write.out <- get(paste0("sens_simu_data_num_cell_",nos_cell_sim,"_2"))
  write.csv   (write.out,file=paste0("./sens_simu_data.rdata_",nos_cell_sim,"_2"))
  save        (write.out,file=paste0("./sens_simu_data.rdata_",nos_cell_sim,"_2"))
}


all.sim.out<-
  rbind(
  cbind(sens_simu_data_num_cell_25 ,rep(25 ,dim(sens_simu_data_num_cell_25 )[1])),
  cbind(sens_simu_data_num_cell_50 ,rep(50 ,dim(sens_simu_data_num_cell_50 )[1])),
  cbind(sens_simu_data_num_cell_75 ,rep(75 ,dim(sens_simu_data_num_cell_75 )[1])),
  cbind(sens_simu_data,rep(100,dim(sens_simu_data)[1])),
  cbind(sens_simu_data_num_cell_125,rep(125,dim(sens_simu_data_num_cell_125)[1])),
  cbind(sens_simu_data_num_cell_150,rep(150,dim(sens_simu_data_num_cell_150)[1])),
  cbind(sens_simu_data_num_cell_175,rep(175,dim(sens_simu_data_num_cell_175)[1])),
  cbind(sens_simu_data_num_cell_200,rep(200,dim(sens_simu_data_num_cell_200)[1]))
  )
all.sim.df<-as.data.frame(all.sim.out)

all.sim2.out<-
  rbind(
    cbind(sens_simu_data_num_cell_25_2 ,rep(25 ,dim(sens_simu_data_num_cell_25_2 )[1])),
    cbind(sens_simu_data_num_cell_50_2 ,rep(50 ,dim(sens_simu_data_num_cell_50_2 )[1])),
    cbind(sens_simu_data_num_cell_75_2 ,rep(75 ,dim(sens_simu_data_num_cell_75_2 )[1])),
    cbind(sens_simu_data_num_cell_100_2,rep(100,dim(sens_simu_data_num_cell_100_2)[1])),
    cbind(sens_simu_data_num_cell_125_2,rep(125,dim(sens_simu_data_num_cell_125_2)[1])),
    cbind(sens_simu_data_num_cell_150_2,rep(150,dim(sens_simu_data_num_cell_150_2)[1])),
    cbind(sens_simu_data_num_cell_175_2,rep(175,dim(sens_simu_data_num_cell_175_2)[1])),
    cbind(sens_simu_data_num_cell_200_2,rep(200,dim(sens_simu_data_num_cell_200_2)[1])),
    cbind(sens_simu_data_num_cell_225,rep(225,dim(sens_simu_data_num_cell_225)[1])),
    cbind(sens_simu_data_num_cell_250,rep(250,dim(sens_simu_data_num_cell_250)[1])),
    cbind(sens_simu_data_num_cell_275,rep(275,dim(sens_simu_data_num_cell_275)[1])),
    cbind(sens_simu_data_num_cell_300,rep(300,dim(sens_simu_data_num_cell_300)[1]))
  )
all.sim2.df<-as.data.frame(all.sim2.out)
save        (all.sim2.out,file="./all_sim2.rdata")

all.sim_oct_4<-rbind(all.sim.df,all.sim2.df)
save        (all.sim_oct_4,file="./all.sim_oct_4")

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="sim_mut") { 
    value=paste0("Mutations: ",value)
  }
  if (var=="sensitivity") { 
    value=paste0("sensitivity: ",as.numeric(value)*100,"%")
  }
  return(value)
}

ggplot(subset(all.sim.df,sensitivity<0.08 & sim_ADO==0.2 & Repli==2),
       aes   (x=V7,y=bic_k,group=factor(sim_ADO),color=factor(bic_k),shape=factor(bic_k)))+
  geom_point(size=3)+
  geom_line (color="grey80")+
  facet_grid(sim_mut~sensitivity,labeller=mf_labeller)+
  scale_y_continuous(breaks=c(1,2))+
  xlab("Number of cells")+
  ylab("Number of detected Clones")+
  scale_color_discrete(name="Nos of Clones",
                       breaks=c("1", "2"),
                       labels=c("1", "2"))+
  scale_shape_discrete(name="Nos of Clones",
                       breaks=c("1", "2"),
                       labels=c("1", "2"))+
  ggtitle("Sensitivity Detection Limit Simulation For ADO=0.2")+
  theme_bw  ()


save        (all.sim.out,file="./all_sim.rdata")

# Generates sub clones for sub-ancestors
my_sub_geno <- function(pt_cluster,nos_grp,grp_sel){
  bin.mat         <- as.matrix(pt_cluster[,-which(colnames(pt_cluster)=="clusters")])
  data.dist.g     <- vegdist  (t(bin.mat)  , method = "jaccard")
  col.clus        <- hclust   (data.dist.g , "ward.D2")
  groups          <- cutree(col.clus, k=nos_grp)
  sub_geno        <- (groups %in% grp_sel)+0
  names(sub_geno) <- names(groups)
  return(sub_geno)
}
#########################################################------
# Iterative Random multiclonal Patient Data for testing #
#########################################################------

rand4.out     <- my_rand_iter()
rand4.df      <- data.frame(rand4.out)
randall.df    <- rbind(randall.df,rand4.df)

simu_plot<- ggplot        (subset(randall.df,sim_mut<55),
               aes   (x=sim_clone, y=aic_k  , group=factor(sim_mut)))+
  geom_point  ()+
  ylab("Clusters inferred using Akaike Information Criterion")+
  xlab("Known Simulated number of Clusters ")+
  facet_grid  (sim_ADO ~ sim_mut,scales="free_y")+
  geom_abline (intercept=0      , slope=1, alpha=0.5,color="brown")+
  geom_smooth (method=lm)+
  theme_bw    ()+
  coord_equal ()
simu_plot
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/simu_plot.pdf",simu_plot,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/simu_plot.svg",simu_plot,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/simu_plot.eps",simu_plot,width=10,height=10)

write.csv(randall.df,"bayes_simu_4.csv")

###############################################------
# Random multiclonal Patient Data for testing #
###############################################------
ptM_random_cl   <- my_rand_clo_bin(nos_clone    = 5   ,
                                   nos_mut      = 10  ,
                                   prob         = 0.5)       # Generate random clone genotypes
ptM_cluster     <- my_rand_clu    (ptM_random_cl                 ,
                                   nos_cell_vec = c(20,25,10,5,15) ,
                                   ADO          = 0.3            ,
                                   noise        = 0.00001)   # Generate random patient data
ptM_mb_clus.df  <- my_create_in_df(ptM_cluster)              # Creates the input for EM
ptM_fmm         <- stepFlexmix    (Incidence ~ 1,
                                   weights = ~ Freq, 
                                   data    = ptM_mb_clus.df,
                                   model   = FLXMCmvbinary(truncated = TRUE ),
                                   control = list         (minprior  = 0.005), 
                                   k       = 1:7, 
                                   nrep    = 5)              # Perform the model base clustering
ptM_ic_plot     <- my_ic_fun      (ptM_fmm)                  # Generate Plots using Information Criterion
dev.off()                                                    # Clear all the plots
ptM_ic_plot                                                  # Do the actual plot
nos_clones      <- 1                                         # IMPT! -> Select the required number of clones from criterion plots
ptM_fmm_best    <- getModel       (ptM_fmm, nos_clones)      # Select the required number of clones
ptM_params_melt <- my_params_df   (ptM_fmm_best)             # Get the Parameters of Best Model
ptM_MDS         <- my_MDS         (ptM_cluster,ptM_fmm_best) # Perform MCA on the binary data
ptM_pos.df      <- posterior      (ptM_fmm_best)             # Get the posterior probability for each cell
nos_row_hc      <- 3                                         # Number of hclust row clusters 
nos_mut         <- 3                                         # Number of mutations categories
ptM_heatmaps    <- my_heatmap     (ptM_cluster,ptM_fmm_best,
                                   nos_row_hc,nos_mut)       # Performs the heatmaps + jaccard dist clustering

ggplot(ptM_params_melt,aes(x=snps,y=value,fill=cat))+
  geom_bar          (stat="identity")+
  geom_hline        (yintercept=0.6  , color="brown",linetype="dashed")+
  facet_wrap        ( ~ variable     , ncol=1)+
  ggtitle           ("Inferred Clone profile")+
  scale_fill_manual (values = c(brewer.pal (length(levels(ptM_params_melt$cat))-1,"Set2"),"grey80"))+
  theme(axis.text.x      = element_text (angle = 90, hjust = 1),
        panel.background = element_blank(),
        axis.line        = element_line (colour = "black"))

ggplot(subset(ptM_params_melt,value>0.5 & value <1),aes(x=1-value))+
  geom_histogram(binwidth = 0.02)+
  ggtitle       ("Inferred Experimental Dropout Rate")+
  theme(axis.text.x       = element_text (angle  = 90, hjust = 1),
        panel.background  = element_blank(),
        axis.line         = element_line (colour = "black"))

ggplot(data = ptM_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(ptM_MDS[[2]]))) + 
  geom_hline    (yintercept = 0, colour = "gray70") + 
  geom_vline    (xintercept = 0, colour = "gray70") +
  geom_point    (aes(color  = factor(cluster_assign)), alpha = 0.7,size=3) +
  geom_density2d(aes(color  = factor(cluster_assign)), alpha = 0.3) +
  geom_text     (data = ptM_MDS[[1]],aes(x = Dim.1, y = Dim.2, label = rownames(ptM_MDS[[1]])),
                 size = 4,alpha = 0.4) +
  scale_colour_manual ( name   ="Clones", 
                        values = brewer.pal(nos_clones, "Set2"),
                        breaks = levels    (ptM_MDS[[2]]$cluster_assign),
                        labels = paste0    ("clone ",levels(ptM_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggplot(data = ptM_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(ptM_MDS[[2]]))) + 
  geom_hline    (yintercept  = 0, colour = "gray70") + 
  geom_vline    (xintercept  = 0, colour = "gray70") +
  geom_point    (aes(color   = factor(cluster_assign)), alpha = 0.7, size = 3) +
  geom_density2d(aes(color   = factor(cluster_assign)), alpha = 0.3) +
  scale_colour_manual(name   ="Clones",
                      values = brewer.pal(nos_clones, "Set2"),
                      breaks = levels    (ptM_MDS[[2]]$cluster_assign),
                      labels = paste0    ("clone ",levels(ptM_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

dev.off()
plot(ptM_heatmaps[[1]])
dev.off()
plot(ptM_heatmaps[[2]])

###################################------
# Random Patient Data for testing #
###################################------

ptr_cluster     <- my_rand_mon_bin(nos_cell= 50 ,
                                   nos_mut = 20 ,
                                   ADO     = 0.3)            # Generate random patient data
ptr_mb_clus.df  <- my_create_in_df(ptr_cluster)              # Creates the input for EM
ptr_fmm         <- stepFlexmix    (Incidence ~ 1,
                                   weights = ~ Freq, 
                                   data    = ptr_mb_clus.df,
                                   model   = FLXMCmvbinary(truncated = TRUE ),
                                   control = list         (minprior  = 0.005), 
                                   k       = 1:7, 
                                   nrep    = 3)              # Perform the model base clustering
ptr_ic_plot     <- my_ic_fun      (ptr_fmm)                  # Generate Plots using Information Criterion
dev.off()                                                    # Clear all the plots
ptr_ic_plot                                                  # Do the actual plot
nos_clones      <- 1                                         # IMPT! -> Select the required number of clones from criterion plots
ptr_fmm_best    <- getModel       (ptr_fmm, nos_clones)      # Select the required number of clones
ptr_params_melt <- my_params_df   (ptr_fmm_best)             # Get the Parameters of Best Model
ptr_MDS         <- my_MDS         (ptr_cluster,ptr_fmm_best) # Perform MCA on the binary data
ptr_pos.df      <- posterior      (ptr_fmm_best)             # Get the posterior probability for each cell
nos_row_hc      <- 3                                         # Number of hclust row clusters 
nos_mut         <- 3                                         # Number of mutations categories
ptr_heatmaps    <- my_heatmap     (ptr_cluster,ptr_fmm_best,
                                   nos_row_hc,nos_mut)       # Performs the heatmaps + jaccard dist clustering

ggplot(ptr_params_melt,aes(x=snps,y=value,fill=cat))+
  geom_bar          (stat="identity")+
  geom_hline        (yintercept=0.6  , color="brown",linetype="dashed")+
  facet_wrap        ( ~ variable     , ncol=1)+
  ggtitle           ("Inferred Clone profile")+
  scale_fill_manual (values = c(brewer.pal (length(levels(ptr_params_melt$cat))-1,"Set2"),"grey80"))+
  theme(axis.text.x      = element_text (angle = 90, hjust = 1),
        panel.background = element_blank(),
        axis.line        = element_line (colour = "black"))

ggplot(subset(ptr_params_melt,value>0.5 & value <1),aes(x=1-value))+
  geom_histogram(binwidth = 0.02)+
  ggtitle       ("Inferred Experimental Dropout Rate")+
  theme(axis.text.x       = element_text (angle  = 90, hjust = 1),
        panel.background  = element_blank(),
        axis.line         = element_line (colour = "black"))

ggplot(data = ptr_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(ptr_MDS[[2]]))) + 
  geom_hline    (yintercept = 0, colour = "gray70") + 
  geom_vline    (xintercept = 0, colour = "gray70") +
  geom_point    (aes(color  = factor(cluster_assign)), alpha = 0.7,size=3) +
  geom_density2d(aes(color  = factor(cluster_assign)), alpha = 0.3) +
  geom_text     (data = ptr_MDS[[1]],aes(x = Dim.1, y = Dim.2, label = rownames(ptr_MDS[[1]])),
                 size = 4,alpha = 0.4) +
  scale_colour_manual ( name   ="Clones", 
                        values = brewer.pal(nos_clones, "Set2"),
                        breaks = levels    (ptr_MDS[[2]]$cluster_assign),
                        labels = paste0    ("clone ",levels(ptr_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggplot(data = ptr_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(ptr_MDS[[2]]))) + 
  geom_hline    (yintercept  = 0, colour = "gray70") + 
  geom_vline    (xintercept  = 0, colour = "gray70") +
  geom_point    (aes(color   = factor(cluster_assign)), alpha = 0.7, size = 3) +
  geom_density2d(aes(color   = factor(cluster_assign)), alpha = 0.3) +
  scale_colour_manual(name   ="Clones",
                      values = brewer.pal(nos_clones, "Set2"),
                      breaks = levels    (ptr_MDS[[2]]$cluster_assign),
                      labels = paste0    ("clone ",levels(ptr_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

dev.off()
plot(ptr_heatmaps[[1]])
dev.off()
plot(ptr_heatmaps[[2]])

#############################----
# Run functions on Patients #
#############################----

###Patient 1--------
pt1_mb_clus.df  <- my_create_in_df(pt1_cluster)              # Creates the input for EM
pt1_fmm         <- stepFlexmix    (Incidence ~ 1,
                                   weights = ~ Freq, 
                                   data    = pt1_mb_clus.df,
                                   model   = FLXMCmvbinary(truncated = TRUE ),
                                   control = list         (minprior  = 0.005), 
                                   k       = 1:7, 
                                   nrep    = 3)              # Perform the model base clustering
pt1_ic_plot     <- my_ic_fun      (pt1_fmm)                  # Generate Plots using Information Criterion
pt1_ic_plot                                                  # Do the actual plot
nos_clones      <- 4                                         # IMPT! -> Select the required number of clones from criterion plots
pt1_fmm_best    <- getModel       (pt1_fmm, nos_clones)      # Select the required number of clones
pt1_params_melt <- my_params_df   (pt1_fmm_best)             # Get the Parameters of Best Model
pt1_MDS         <- my_MDS         (pt1_cluster,pt1_fmm_best) # Perform MCA on the binary data
pt1_pos.df      <- posterior      (pt1_fmm_best)             # Get the posterior probability for each cell
nos_row_hc      <- 4                                         # Number of hclust row clusters 
nos_mut         <- 3                                         # Number of mutations categories
pt1_heatmaps    <- my_heatmap     (pt1_cluster,pt1_fmm_best,
                                   nos_row_hc,nos_mut)       # Performs the heatmaps + jaccard dist clustering
pt1_est_ado     <- median         (1-pt1_params_melt$value[  pt1_params_melt$value >  0.5 
                                                           & pt1_params_melt$value <= 1  ])


pt1_gg_clone_profile <- ggplot(pt1_params_melt,aes(x=snps,y=value,fill=cat))+
  geom_bar          (stat="identity")+
  geom_hline        (yintercept=0.6  , color="brown",linetype="dashed")+
  facet_wrap        ( ~ variable     , ncol=1)+
  ggtitle           ("Inferred Clone profile")+
  scale_fill_manual (values = c(brewer.pal (length(levels(pt1_params_melt$cat))-1,"Set1"),"grey80"))+
  theme(axis.text.x      = element_text (angle = 90, hjust = 1),
        panel.background = element_blank(),
        axis.line        = element_line (colour = "black"))
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_clone_profile.pdf",pt1_gg_clone_profile,width=7.5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_clone_profile.svg",pt1_gg_clone_profile,width=7.5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_clone_profile.eps",pt1_gg_clone_profile,width=7.5,height=10)

pt1_gg_ADO <- ggplot(subset(pt1_params_melt,value>0.5 & value <1),aes(x=1-value))+
  geom_histogram(binwidth = 0.02)+
  ggtitle       ("Inferred Experimental Dropout Rate")+
  xlab          ("Allele Dropout Rate") +
  theme(axis.text.x       = element_text (angle  = 90, hjust = 1),
        panel.background  = element_blank(),
        axis.line         = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_ADO.pdf",pt1_gg_ADO)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_ADO.svg",pt1_gg_ADO)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_ADO.eps",pt1_gg_ADO)


pt1_gg_MCA_1 <- ggplot(data = pt1_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt1_MDS[[2]]))) + 
  geom_hline    (yintercept = 0, colour = "gray70") + 
  geom_vline    (xintercept = 0, colour = "gray70") +
  geom_point    (aes(color  = factor(cluster_assign)), alpha = 0.7,size=6) +
  geom_density2d(aes(color  = factor(cluster_assign)), alpha = 0.3) +
  geom_text     (data = pt1_MDS[[1]],aes(x = Dim.1, y = Dim.2, label = rownames(pt1_MDS[[1]])),
                 size = 5,alpha = 0.4) +
  scale_colour_manual ( name   ="Clones", 
                        values = brewer.pal(nos_clones, "Set2"),
                        breaks = levels    (pt1_MDS[[2]]$cluster_assign),
                        labels = paste0    ("clone ",levels(pt1_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_MCA_1.pdf",pt1_gg_MCA_1,scale=2.5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_MCA_1.svg",pt1_gg_MCA_1,scale=2.5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_MCA_1.eps",pt1_gg_MCA_1,scale=2.5)

pt1_gg_MCA <- ggplot(data = pt1_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt1_MDS[[2]]))) + 
  geom_hline    (yintercept  = 0, colour = "gray70") + 
  geom_vline    (xintercept  = 0, colour = "gray70") +
  geom_point    (aes(color   = factor(cluster_assign)), alpha = 0.7, size = 6) +
  geom_density2d(aes(color   = factor(cluster_assign)), alpha = 0.3) +
  scale_colour_manual(name   ="Clones",
                      values = brewer.pal(nos_clones, "Set2"),
                      breaks = levels    (pt1_MDS[[2]]$cluster_assign),
                      labels = paste0    ("clone ",levels(pt1_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_MCA.pdf",pt1_gg_MCA,scale=2.5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_MCA.svg",pt1_gg_MCA,scale=2.5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_MCA.eps",pt1_gg_MCA,scale=2.5)

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_em_clus_heatmap.pdf",width=5,height=7)
plot(pt1_heatmaps[[1]])
dev.off()

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_em_compare_heatmap.pdf",width=5,height=7)
plot(pt1_heatmaps[[2]])
dev.off()

pt1_cl_geno   <- my_param_2_cl (pt1_fmm_best)
pt1_cl_sel    <- c(1,2,3)
pt1_pre_tree  <- my_pre_tree   (pt1_cl_geno[pt1_cl_sel,],
                                table(flexmix::clusters(pt1_fmm_best))[pt1_cl_sel],
                                primary=TRUE)
pt1_tree      <- my_run_tree   (pt1_pre_tree[[1]]            ,
                                pt1_pre_tree[[2]])

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt1_tree.pdf",width=7,height=7)
plot(pt1_tree[[1]],
     layout          = my_layout,
     vertex.size     = c(0,as.vector(table(flexmix::clusters(pt1_fmm_best))[pt1_cl_sel]))+8,
     vertex.label    = rownames(pt1_tree[[2]]),
     vertex.color    = c("grey80",brewer.pal(length(grep("_0",rownames(pt1_pre_tree[[1]]),invert=TRUE)),"Set2")),
     edge.arrow.size = 0.3,
     edge.color      = "black")
dev.off()

###Patient 2--------
pt2_mb_clus.df  <- my_create_in_df(pt2_cluster)              # Creates the input for EM
pt2_fmm         <- stepFlexmix    (Incidence ~ 1,
                                   weights = ~ Freq, 
                                   data    = pt2_mb_clus.df,
                                   model   = FLXMCmvbinary(truncated = TRUE ),
                                   control = list         (minprior  = 0.005), 
                                   k       = 1:7, 
                                   nrep    = 3)              # Perform the model base clustering
pt2_ic_plot     <- my_ic_fun      (pt2_fmm)                  # Generate Plots using Information Criterion
pt2_ic_plot                                                  # Do the actual plot
nos_clones      <- 5                                         # IMPT! -> Select the required number of clones from criterion plots
pt2_fmm_best    <- getModel       (pt2_fmm, nos_clones)      # Select the required number of clones
pt2_params_melt <- my_params_df   (pt2_fmm_best)             # Get the Parameters of Best Model
pt2_MDS         <- my_MDS         (pt2_cluster,pt2_fmm_best) # Perform MCA on the binary data
pt2_pos.df      <- posterior      (pt2_fmm_best)             # Get the posterior probability for each cell
nos_row_hc      <- 6                                         # Number of hclust row clusters 
nos_mut         <- 4                                         # Number of mutations categories
pt2_heatmaps    <- my_heatmap     (pt2_cluster,pt2_fmm_best,
                                   nos_row_hc,nos_mut)       # Performs the heatmaps + jaccard dist clustering
pt2_est_ado     <- median         (1-pt2_params_melt$value[  pt2_params_melt$value >  0.5 
                                                           & pt2_params_melt$value <= 1  ])


pt2_gg_clone_profile <- ggplot(pt2_params_melt,aes(x=snps,y=value,fill=cat))+
  geom_bar          (stat="identity")+
  geom_hline        (yintercept=0.6  , color="brown",linetype="dashed")+
  facet_wrap        ( ~ variable     , ncol=1)+
  ggtitle           ("Inferred Clone profile")+
  scale_fill_manual (values = c(brewer.pal (length(levels(pt2_params_melt$cat))-1,"Set1"),"grey80"))+
  theme(axis.text.x         = element_text (angle = 90, hjust = 1),
        panel.background    = element_blank(),
        axis.line           = element_line (colour = "black"))
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_clone_profile.pdf",pt2_gg_clone_profile,width=5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_clone_profile.svg",pt2_gg_clone_profile,width=5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_clone_profile.eps",pt2_gg_clone_profile,width=5,height=10)

pt2_gg_ADO <- ggplot(subset(pt2_params_melt,value>0.5 & value <1),aes(x=1-value))+
  geom_histogram(binwidth = 0.02)+
  ggtitle       ("Inferred Experimental Dropout Rate")+
  xlab          ("Allele Dropout Rate") +
  theme(axis.text.x       = element_text (angle  = 90, hjust = 1),
        panel.background  = element_blank(),
        axis.line         = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_ADO.pdf",pt2_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_ADO.svg",pt2_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_ADO.eps",pt2_gg_ADO,width=5,height=5)

pt2_gg_MCA_1 <- ggplot(data = pt2_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt2_MDS[[2]]))) + 
  geom_hline    (yintercept = 0, colour = "gray70") + 
  geom_vline    (xintercept = 0, colour = "gray70") +
  geom_point    (aes(color  = factor(cluster_assign)), alpha = 0.7,size=6) +
  geom_density2d(aes(color  = factor(cluster_assign)), alpha = 0.3) +
  geom_text     (data = pt2_MDS[[1]],aes(x = Dim.1, y = Dim.2, label = rownames(pt2_MDS[[1]])),
                 size = 4,alpha = 0.4) +
  scale_colour_manual ( name   ="Clones", 
                        values = brewer.pal(nos_clones, "Set2"),
                        breaks = levels    (pt2_MDS[[2]]$cluster_assign),
                        labels = paste0    ("clone ",levels(pt2_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_MCA_1.pdf",pt2_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_MCA_1.svg",pt2_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_MCA_1.eps",pt2_gg_MCA_1,width=10,height=10)

pt2_gg_MCA<-ggplot(data = pt2_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt2_MDS[[2]]))) + 
  geom_hline    (yintercept  = 0, colour = "gray70") + 
  geom_vline    (xintercept  = 0, colour = "gray70") +
  geom_point    (aes(color   = factor(cluster_assign)), alpha = 0.7, size = 5) +
  geom_density2d(aes(color   = factor(cluster_assign)), alpha = 0.3) +
  scale_colour_manual(name   ="Clones",
                      values = brewer.pal(nos_clones, "Set2"),
                      breaks = levels    (pt2_MDS[[2]]$cluster_assign),
                      labels = paste0    ("clone ",levels(pt2_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_MCA.pdf",pt2_gg_MCA,width=18,height=12)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_MCA.svg",pt2_gg_MCA,width=18,height=12)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_MCA.eps",pt2_gg_MCA,width=18,height=12)

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_em_clus_heatmap.pdf",width=5,height=7)
plot(pt2_heatmaps[[1]])
dev.off()

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_em_compare_heatmap.pdf",width=5,height=7)
plot(pt2_heatmaps[[2]])
dev.off()

pt2_cl_geno   <- my_param_2_cl (pt2_fmm_best)
pt2_cl_sel    <- c(1,2,4,5)
pt2_pre_tree  <- my_pre_tree   (pt2_cl_geno[pt2_cl_sel,],
                                table(flexmix::clusters(pt2_fmm_best))[pt2_cl_sel],
                                primary=FALSE)
pt2_tree      <- my_run_tree   (pt2_pre_tree[[1]]            ,
                                pt2_pre_tree[[2]])
dev.off()
pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_tree.pdf",width=7,height=7)
plot(pt2_tree[[1]],
     layout          = my_layout(pt2_tree[[1]], root=2),
     vertex.size     = c(as.vector(table(flexmix::clusters(pt2_fmm_best)))[pt2_cl_sel])*0.75,
     vertex.label    = rownames(pt2_tree[[2]]),
     vertex.color    = c(brewer.pal(length(rownames(pt2_cl_geno)),"Set2")[pt2_cl_sel]),
     edge.arrow.size = 0.3,
     edge.color      = "black")
dev.off()

#for patient 2 we add in a sub clone ancestor

sub           <- my_sub_geno   (pt2_cluster ,6 , c(4,2))
pt2_cl_geno   <- my_param_2_cl (pt2_fmm_best)
pt2_cl_sel    <- c(1,2,4,5)
pt2_cl_geno_p <- rbind         (pt2_cl_geno[pt2_cl_sel,],sub)
row.names(pt2_cl_geno_p)[length(pt2_cl_geno_p[,1])]<-"Clone_sub"
pt2_cl_geno_s <- c(table(flexmix::clusters(pt2_fmm_best))[pt2_cl_sel],0)
pt2_pre_tr_s  <- my_pre_tree   (pt2_cl_geno_p,
                                pt2_cl_geno_s,
                                primary=FALSE)
pt2_tr_s      <- my_run_tree   (pt2_pre_tr_s[[1]]            ,
                                pt2_pre_tr_s[[2]])

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt2_tree_subClone.pdf",width=7,height=7)
plot(pt2_tr_s[[1]],
     layout          = my_layout(pt2_tr_s[[1]], root=2),
     vertex.size     = pt2_cl_geno_s,
     vertex.label    = rownames(pt2_tr_s[[2]]),
     vertex.color    = c(brewer.pal(length(rownames(pt2_cl_geno)),"Set2")[pt2_cl_sel],"grey30"),
     edge.arrow.size = 0.3,
     edge.color      = "black")
dev.off()

###Patient 3--------
pt3_mb_clus.df  <- my_create_in_df(pt3_cluster)              # Creates the input for EM
pt3_fmm         <- stepFlexmix    (Incidence ~ 1,
                                   weights = ~ Freq, 
                                   data    = pt3_mb_clus.df,
                                   model   = FLXMCmvbinary(truncated = TRUE ),
                                   control = list         (minprior  = 0.005), 
                                   k       = 1:7, 
                                   nrep    = 3)              # Perform the model base clustering
pt3_ic_plot     <- my_ic_fun      (pt3_fmm)                  # Generate Plots using Information Criterion
pt3_ic_plot                                                  # Do the actual plot
nos_clones      <- 5                                         # IMPT! -> Select the required number of clones from criterion plots
pt3_fmm_best    <- getModel       (pt3_fmm, nos_clones)      # Select the required number of clones
pt3_params_melt <- my_params_df   (pt3_fmm_best)             # Get the Parameters of Best Model
pt3_MDS         <- my_MDS         (pt3_cluster,pt3_fmm_best) # Perform MCA on the binary data
pt3_pos.df      <- posterior      (pt3_fmm_best)             # Get the posterior probability for each cell
nos_row_hc      <- 5                                         # Number of hclust row clusters 
nos_mut         <- 5                                         # Number of mutations categories
pt3_heatmaps    <- my_heatmap     (pt3_cluster,pt3_fmm_best,
                                   nos_row_hc,nos_mut)       # Performs the heatmaps + jaccard dist clustering
pt3_est_ado     <- median         (1-pt3_params_melt$value[  pt3_params_melt$value >  0.5 
                                                           & pt3_params_melt$value <= 1  ])

pt3_gg_clone_profile<-ggplot(pt3_params_melt,aes(x=snps,y=value,fill=cat))+
  geom_bar          (stat="identity")+
  geom_hline        (yintercept=0.6  , color="brown",linetype="dashed")+
  facet_wrap        ( ~ variable     , ncol=1)+
  ggtitle           ("Inferred Clone profile")+
  scale_fill_manual (values = c(brewer.pal (length(levels(pt3_params_melt$cat))-1,"Set1"),"grey80"))+
  theme(axis.text.x         = element_text (angle = 90, hjust = 1),
        panel.background    = element_blank(),
        axis.line           = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_clone_profile.pdf",pt3_gg_clone_profile,width=7.5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_clone_profile.svg",pt3_gg_clone_profile,width=7.5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_clone_profile.eps",pt3_gg_clone_profile,width=7.5,height=10)

pt3_gg_ADO <- ggplot(subset(pt3_params_melt,value>0.5 & value <1),aes(x=1-value))+
  geom_histogram(binwidth = 0.02)+
  ggtitle       ("Inferred Experimental Dropout Rate")+
  xlab          ("Allele Dropout Rate") +
  theme(axis.text.x       = element_text (angle  = 90, hjust = 1),
        panel.background  = element_blank(),
        axis.line         = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_ADO.pdf",pt3_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_ADO.svg",pt3_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_ADO.eps",pt3_gg_ADO,width=5,height=5)

pt3_gg_MCA_1 <- ggplot(data = pt3_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt3_MDS[[2]]))) + 
  geom_hline    (yintercept = 0, colour = "gray70") + 
  geom_vline    (xintercept = 0, colour = "gray70") +
  geom_point    (aes(color  = factor(cluster_assign)), alpha = 0.7,size=6) +
  geom_density2d(aes(color  = factor(cluster_assign)), alpha = 0.3) +
  geom_text     (data = pt3_MDS[[1]],aes(x = Dim.1, y = Dim.2, label = rownames(pt3_MDS[[1]])),
                 size = 4,alpha = 0.4) +
  scale_colour_manual ( name   ="Clones", 
                        values = brewer.pal(nos_clones, "Set2"),
                        breaks = levels    (pt3_MDS[[2]]$cluster_assign),
                        labels = paste0    ("clone ",levels(pt3_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_MCA_1.pdf",pt3_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_MCA_1.svg",pt3_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_MCA_1.eps",pt3_gg_MCA_1,width=10,height=10)

pt3_gg_MCA <- ggplot(data = pt3_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt3_MDS[[2]]))) + 
  geom_hline    (yintercept  = 0, colour = "gray70") + 
  geom_vline    (xintercept  = 0, colour = "gray70") +
  geom_point    (aes(color   = factor(cluster_assign)), alpha = 0.7, size = 6) +
  geom_density2d(aes(color   = factor(cluster_assign)), alpha = 0.3) +
  scale_colour_manual(name   ="Clones",
                      values = brewer.pal(nos_clones, "Set2"),
                      breaks = levels    (pt3_MDS[[2]]$cluster_assign),
                      labels = paste0    ("clone ",levels(pt3_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_MCA.pdf",pt3_gg_MCA,width=18,height=12)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_MCA.svg",pt3_gg_MCA,width=18,height=12)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_MCA.eps",pt3_gg_MCA,width=18,height=12)

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_em_clus_heatmap.pdf",width=5,height=7)
plot(pt3_heatmaps[[1]])
dev.off()

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_em_compare_heatmap.pdf",width=5,height=7)
plot(pt3_heatmaps[[2]])
dev.off()

pt3_cl_geno   <- my_param_2_cl (pt3_fmm_best)
pt3_cl_sel    <- c(1,2,3,4,5)
pt3_pre_tree  <- my_pre_tree   (pt3_cl_geno[pt3_cl_sel,],
                                table(flexmix::clusters(pt3_fmm_best))[pt3_cl_sel],
                                primary=TRUE)
pt3_tree      <- my_run_tree   (pt3_pre_tree[[1]],
                                pt3_pre_tree[[2]])

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt3_tree.pdf",width=7,height=7)
plot(pt3_tree[[1]],
     layout          = my_layout(pt3_tree[[1]], root=1),
     vertex.size     = c(1,as.vector(table(flexmix::clusters(pt3_fmm_best)))[pt3_cl_sel])*0.75,
     vertex.label    = rownames(pt3_tree[[2]]),
     vertex.color    = c("grey30",brewer.pal(length(rownames(pt3_cl_geno)),"Set2")[pt3_cl_sel]),
     edge.arrow.size = 0.3,
     edge.color      = "black")
dev.off()

###Patient 4--------

# Perform model based clustering
pt4_mb_clus.df  <- my_create_in_df(pt4_cluster)              # Creates the input for EM
pt4_fmm         <- stepFlexmix    (Incidence ~ 1,
                                   weights = ~ Freq, 
                                   data    = pt4_mb_clus.df,
                                   model   = FLXMCmvbinary(truncated = TRUE ),
                                   control = list         (minprior  = 0.005), 
                                   k       = 1:8, 
                                   nrep    = 3)              # Perform the model base clustering

# Perform information criterion
pt4_ic_plot     <- my_ic_fun      (pt4_fmm)                  # Generate Plots using Information Criterion
pt4_ic_plot                                                  # Do the actual plot

# Model Selection
nos_clones      <- 5                                         # IMPT! -> Select the required number of clones from criterion plots
pt4_fmm_best    <- getModel       (pt4_fmm, nos_clones)      # Select the required number of clones
pt4_params_melt <- my_params_df   (pt4_fmm_best)             # Get the Parameters of Best Model
pt4_pos.df      <- posterior      (pt4_fmm_best)             # Get the posterior probability for each cell

# MDS analysis
pt4_MDS         <- my_MDS         (pt4_cluster,pt4_fmm_best) # Perform MCA on the binary data

# Hclust Analysis
nos_row_hc      <- 4                                         # Number of hclust row clusters 
nos_mut         <- 4                                         # Number of mutations categories
pt4_heatmaps    <- my_heatmap     (pt4_cluster,pt4_fmm_best,
                                   nos_row_hc,nos_mut)       # Performs the heatmaps + jaccard dist clustering
View(my_pvalue_clone(pt4_cluster,pt4_fmm_best))

pt4_est_ado     <- median         (1-pt4_params_melt$value[  pt4_params_melt$value >  0.5 
                                                          &  pt4_params_melt$value <= 1  ])

pt4_gg_clone_profile <- ggplot(pt4_params_melt,aes(x=snps,y=value,fill=cat))+
  geom_bar          (stat="identity")+
  geom_hline        (yintercept=0.6  , color="brown",linetype="dashed")+
  facet_wrap        ( ~ variable     , ncol=1)+
  ggtitle           ("Inferred Clone profile")+
  scale_fill_manual (values = c(brewer.pal (length(levels(pt4_params_melt$cat))-1,"Set1"),"grey80"))+
  theme(axis.text.x         = element_text (angle = 90, hjust = 1),
        panel.background    = element_blank(),
        axis.line           = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_clone_profile.pdf",pt4_gg_clone_profile,width=11.5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_clone_profile.svg",pt4_gg_clone_profile,width=11.5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_clone_profile.eps",pt4_gg_clone_profile,width=11.5,height=10)


pt4_gg_ADO<-ggplot(subset(pt4_params_melt,value>0.5 & value <1),aes(x=1-value))+
  geom_histogram(binwidth = 0.02)+
  ggtitle       ("Inferred Experimental Dropout Rate")+
  xlab          ("Allele Dropout Rate") +
  theme(axis.text.x       = element_text (angle  = 90, hjust = 1),
        panel.background  = element_blank(),
        axis.line         = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_ADO.pdf",pt4_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_ADO.svg",pt4_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_ADO.eps",pt4_gg_ADO,width=5,height=5)

pt4_gg_MCA_1 <- ggplot(data = pt4_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt4_MDS[[2]]))) + 
  geom_hline    (yintercept = 0, colour = "gray70") + 
  geom_vline    (xintercept = 0, colour = "gray70") +
  geom_point    (aes(color  = factor(cluster_assign)), alpha = 0.7,size=5) +
  geom_density2d(aes(color  = factor(cluster_assign)), alpha = 0.3) +
  geom_text     (data = pt4_MDS[[1]],aes(x = Dim.1, y = Dim.2, label = rownames(pt4_MDS[[1]])),
                 size = 4,alpha = 0.4) +
  scale_colour_manual ( name   ="Clones", 
                        values = brewer.pal(nos_clones, "Set2"),
                        breaks = levels    (pt4_MDS[[2]]$cluster_assign),
                        labels = paste0    ("clone ",levels(pt4_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_MCA_1.pdf",pt4_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_MCA_1.svg",pt4_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_MCA_1.eps",pt4_gg_MCA_1,width=10,height=10)

pt4_gg_MCA <- ggplot(data = pt4_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt4_MDS[[2]]))) + 
  geom_hline    (yintercept  = 0, colour = "gray70") + 
  geom_vline    (xintercept  = 0, colour = "gray70") +
  geom_point    (aes(color   = factor(cluster_assign)), alpha = 0.7, size = 5) +
  geom_density2d(aes(color   = factor(cluster_assign)), alpha = 0.3) +
  scale_colour_manual(name   ="Clones",
                      values = brewer.pal(nos_clones, "Set2"),
                      breaks = levels    (pt4_MDS[[2]]$cluster_assign),
                      labels = paste0    ("clone ",levels(pt4_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_MCA.pdf",pt4_gg_MCA,width=18,height=12)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_MCA.svg",pt4_gg_MCA,width=18,height=12)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_MCA.eps",pt4_gg_MCA,width=18,height=12)

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_em_clus_heatmap.pdf",width=11.5,height=14)
plot(pt4_heatmaps[[1]])
dev.off()

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_em_compare_heatmap.pdf",width=8.5,height=10)
plot(pt4_heatmaps[[2]])
dev.off()

pt4_cl_geno   <- my_param_2_cl (pt4_fmm_best)
pt4_cl_sel    <- c(1,2,3,5)
pt4_pre_tree  <- my_pre_tree   (pt4_cl_geno[pt4_cl_sel,],
                                table(flexmix::clusters(pt4_fmm_best))[pt4_cl_sel],
                                primary=TRUE)
pt4_pre_tree[[2]][5]<-pt4_pre_tree[[2]][5]+2
pt4_tree      <- my_run_tree   (pt4_pre_tree[[1]],
                                pt4_pre_tree[[2]])

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt4_tree.pdf",width=7,height=7)
plot(pt4_tree[[1]],
     layout          = my_layout(pt4_tree[[1]], root=1),
     vertex.size     = c(1,as.vector(table(flexmix::clusters(pt4_fmm_best)))[pt4_cl_sel])*0.75,
     vertex.label    = rownames(pt4_tree[[2]]),
     vertex.color    = c("grey30",brewer.pal(length(rownames(pt4_cl_geno)),"Set2")[pt4_cl_sel]),
     edge.arrow.size = 0.3,
     edge.color      = "black")
dev.off()

###Patient 5--------

# Perform model based clustering
pt5_mb_clus.df  <- my_create_in_df(pt5_cluster)              # Creates the input for EM
pt5_fmm         <- stepFlexmix    (Incidence ~ 1,
                                   weights = ~ Freq, 
                                   data    = pt5_mb_clus.df,
                                   model   = FLXMCmvbinary(truncated = TRUE ),
                                   control = list         (minprior  = 0.005), 
                                   k       = 1:8, 
                                   nrep    = 5)              # Perform the model base clustering

# Perform information criterion
pt5_ic_plot     <- my_ic_fun      (pt5_fmm)                  # Generate Plots using Information Criterion
pt5_ic_plot                                                  # Do the actual plot

# Model Selection
nos_clones      <- 4                                         # IMPT! -> Select the required number of clones from criterion plots
pt5_fmm_best    <- getModel       (pt5_fmm, nos_clones)      # Select the required number of clones
pt5_params_melt <- my_params_df   (pt5_fmm_best)             # Get the Parameters of Best Model
pt5_pos.df      <- posterior      (pt5_fmm_best)             # Get the posterior probability for each cell

# MDS analysis
pt5_MDS         <- my_MDS         (pt5_cluster,pt5_fmm_best) # Perform MCA on the binary data

# Hclust Analysis
nos_row_hc      <- 5                                         # Number of hclust row clusters 
nos_mut         <- 4                                         # Number of mutations categories
pt5_heatmaps    <- my_heatmap     (pt5_cluster,pt5_fmm_best,
                                   nos_row_hc,nos_mut)       # Performs the heatmaps + jaccard dist clustering
pt5_est_ado     <- median         (1-pt5_params_melt$value[  pt5_params_melt$value >  0.5 
                                                           & pt5_params_melt$value <= 1  ])
View(my_pvalue_clone(pt5_cluster,pt5_fmm_best))

pt5_gg_clone_profile <- ggplot(pt5_params_melt,aes(x=snps,y=value,fill=cat))+
  geom_bar          (stat="identity")+
  geom_hline        (yintercept=0.6  , color="brown",linetype="dashed")+
  facet_wrap        ( ~ variable     , ncol=1)+
  ggtitle           ("Inferred Clone profile")+
  scale_fill_manual (values = c(brewer.pal (length(levels(pt5_params_melt$cat))-1,"Set1"),"grey80"))+
  theme(axis.text.x         = element_text (angle = 90, hjust = 1),
        panel.background    = element_blank(),
        axis.line           = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_clone_profile.pdf",pt5_gg_clone_profile,width=14,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_clone_profile.svg",pt5_gg_clone_profile,width=14,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_clone_profile.eps",pt5_gg_clone_profile,width=14,height=10)


pt5_gg_ADO<-ggplot(subset(pt5_params_melt,value>0.5 & value <1),aes(x=1-value))+
  geom_histogram(binwidth = 0.02)+
  ggtitle       ("Inferred Experimental Dropout Rate")+
  xlab          ("Allele Dropout Rate") +
  theme(axis.text.x       = element_text (angle  = 90, hjust = 1),
        panel.background  = element_blank(),
        axis.line         = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_ADO.pdf",pt5_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_ADO.svg",pt5_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_ADO.eps",pt5_gg_ADO,width=5,height=5)

pt5_gg_MCA_1<-ggplot(data = pt5_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt5_MDS[[2]]))) + 
  geom_hline    (yintercept = 0, colour = "gray70") + 
  geom_vline    (xintercept = 0, colour = "gray70") +
  geom_point    (aes(color  = factor(cluster_assign)), alpha = 0.7,size=6) +
  geom_density2d(aes(color  = factor(cluster_assign)), alpha = 0.3) +
  geom_text     (data = pt5_MDS[[1]],aes(x = Dim.1, y = Dim.2, label = rownames(pt5_MDS[[1]])),
                 size = 4,alpha = 0.4) +
  scale_colour_manual ( name   ="Clones", 
                        values = brewer.pal(nos_clones, "Set2"),
                        breaks = levels    (pt5_MDS[[2]]$cluster_assign),
                        labels = paste0    ("clone ",levels(pt5_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_MCA_1.pdf",pt5_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_MCA_1.svg",pt5_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_MCA_1.eps",pt5_gg_MCA_1,width=10,height=10)

pt5_gg_MCA<-ggplot(data = pt5_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt5_MDS[[2]]))) + 
  geom_hline    (yintercept  = 0, colour = "gray70") + 
  geom_vline    (xintercept  = 0, colour = "gray70") +
  geom_point    (aes(color   = factor(cluster_assign)), alpha = 0.7, size = 6) +
  geom_density2d(aes(color   = factor(cluster_assign)), alpha = 0.3) +
  scale_colour_manual(name   ="Clones",
                      values = brewer.pal(nos_clones, "Set2"),
                      breaks = levels    (pt5_MDS[[2]]$cluster_assign),
                      labels = paste0    ("clone ",levels(pt5_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_MCA.pdf",pt5_gg_MCA,width=18,height=15)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_MCA.svg",pt5_gg_MCA,width=18,height=15)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_MCA.eps",pt5_gg_MCA,width=18,height=15)

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_em_clus_heatmap.pdf",width=8,height=11)
plot(pt5_heatmaps[[1]])
dev.off()

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_em_compare_heatmap.pdf",width=8,height=11)
plot(pt5_heatmaps[[2]])
dev.off()

pt5_cl_geno   <- my_param_2_cl (pt5_fmm_best)
pt5_cl_sel    <- c(2,3,4)
pt5_pre_tree  <- my_pre_tree   (pt5_cl_geno[pt5_cl_sel,],
                                table(flexmix::clusters(pt5_fmm_best))[pt5_cl_sel],
                                primary=TRUE)

pt5_tree      <- my_run_tree   (pt5_pre_tree[[1]],
                                pt5_pre_tree[[2]])

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt5_tree.pdf",width=7,height=7)
plot(pt5_tree[[1]],
     layout          = my_layout(pt5_tree[[1]], root=1),
     vertex.size     = c(1,as.vector(table(flexmix::clusters(pt5_fmm_best)))[pt5_cl_sel])*0.75,
     vertex.label    = rownames(pt5_tree[[2]]),
     vertex.color    = c("grey30",brewer.pal(length(rownames(pt5_cl_geno)),"Set2")[pt5_cl_sel]),
     edge.arrow.size = 0.3,
     edge.color      = "black")
dev.off()

###Patient 6--------

# Perform model based clustering
pt6_mb_clus.df  <- my_create_in_df(pt6_cluster)              # Creates the input for EM
pt6_fmm         <- stepFlexmix    (Incidence ~ 1,
                                   weights = ~ Freq, 
                                   data    = pt6_mb_clus.df,
                                   model   = FLXMCmvbinary(truncated = TRUE ),
                                   control = list         (minprior  = 0.005), 
                                   k       = 1:7, 
                                   nrep    = 5)              # Perform the model base clustering

# Perform information criterion
pt6_ic_plot     <- my_ic_fun      (pt6_fmm)                  # Generate Plots using Information Criterion
pt6_ic_plot                                                  # Do the actual plot

# Model Selection
nos_clones      <- 2                                         # IMPT! -> Select the required number of clones from criterion plots
pt6_fmm_best    <- getModel       (pt6_fmm, nos_clones)      # Select the required number of clones
pt6_params_melt <- my_params_df   (pt6_fmm_best)             # Get the Parameters of Best Model
pt6_pos.df      <- posterior      (pt6_fmm_best)             # Get the posterior probability for each cell

# MDS analysis
pt6_MDS         <- my_MDS         (pt6_cluster,pt6_fmm_best) # Perform MCA on the binary data

# Hclust Analysis
nos_row_hc      <- 2                                         # Number of hclust row clusters 
nos_mut         <- 2                                         # Number of mutations categories
pt6_heatmaps    <- my_heatmap     (pt6_cluster,pt6_fmm_best,
                                   nos_row_hc,nos_mut)       # Performs the heatmaps + jaccard dist clustering
pt6_est_ado     <- median         (1-pt6_params_melt$value[  pt6_params_melt$value >  0.5 
                                                             & pt6_params_melt$value <= 1  ])

pt6_gg_clone_profile<-ggplot(pt6_params_melt,aes(x=snps,y=value,fill=cat))+
  geom_bar          (stat="identity")+
  geom_hline        (yintercept=0.6  , color="brown",linetype="dashed")+
  facet_wrap        ( ~ variable     , ncol=1)+
  ggtitle           ("Inferred Clone profile")+
  scale_fill_manual (values = c(brewer.pal (length(levels(pt6_params_melt$cat))-1,"Set1"),"grey80"))+
  theme(axis.text.x         = element_text (angle = 90, hjust = 1),
        panel.background    = element_blank(),
        axis.line           = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_clone_profile.pdf",pt6_gg_clone_profile,width=5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_clone_profile.svg",pt6_gg_clone_profile,width=5,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_clone_profile.eps",pt6_gg_clone_profile,width=5,height=10)

pt6_gg_ADO <- ggplot(subset(pt6_params_melt,value>0.5 & value <1),aes(x=1-value))+
  geom_histogram(binwidth = 0.02)+
  ggtitle       ("Inferred Experimental Dropout Rate")+
  theme(axis.text.x       = element_text (angle  = 90, hjust = 1),
        panel.background  = element_blank(),
        axis.line         = element_line (colour = "black"))

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_ADO.pdf",pt6_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_ADO.svg",pt6_gg_ADO,width=5,height=5)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_ADO.eps",pt6_gg_ADO,width=5,height=5)

pt6_gg_MCA_1<-ggplot(data = pt6_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt6_MDS[[2]]))) + 
  geom_hline    (yintercept = 0, colour = "gray70") + 
  geom_vline    (xintercept = 0, colour = "gray70") +
  geom_point    (aes(color  = factor(cluster_assign)), alpha = 0.7,size=6) +
  geom_density2d(aes(color  = factor(cluster_assign)), alpha = 0.3) +
  geom_text     (data = pt6_MDS[[1]],aes(x = Dim.1, y = Dim.2, label = rownames(pt6_MDS[[1]])),
                 size = 4,alpha = 0.4) +
  scale_colour_manual ( name   ="Clones", 
                        values = brewer.pal(nos_clones, "Set2"),
                        breaks = levels    (pt6_MDS[[2]]$cluster_assign),
                        labels = paste0    ("clone ",levels(pt6_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_MCA_1.pdf",pt6_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_MCA_1.svg",pt6_gg_MCA_1,width=10,height=10)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_MCA_1.eps",pt6_gg_MCA_1,width=10,height=10)

pt6_gg_MCA<-ggplot(data = pt6_MDS[[2]], aes(x = Dim.1, y = Dim.2, label = rownames(pt6_MDS[[2]]))) + 
  geom_hline    (yintercept  = 0, colour = "gray70") + 
  geom_vline    (xintercept  = 0, colour = "gray70") +
  geom_point    (aes(color   = factor(cluster_assign)), alpha = 0.7, size = 6) +
  geom_density2d(aes(color   = factor(cluster_assign)), alpha = 0.3) +
  scale_colour_manual(name   ="Clones",
                      values = brewer.pal(nos_clones, "Set2"),
                      breaks = levels    (pt6_MDS[[2]]$cluster_assign),
                      labels = paste0    ("clone ",levels(pt6_MDS[[2]]$cluster_assign)))+
  ggtitle       ("MCA plot Single Cell Data With Finite Mixed Model Based Clustering") +
  theme_bw()

ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_MCA.pdf",pt6_gg_MCA,width=18,height=12)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_MCA.svg",pt6_gg_MCA,width=18,height=12)
ggsave("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_MCA.eps",pt6_gg_MCA,width=18,height=12)

pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_em_clus_heatmap.pdf",width=15,height=19)
plot(pt6_heatmaps[[1]])
dev.off()
pdf("/datastore/winstonk/chuck_clonal_plots_EM_cluster/pt6_em_compare_heatmap.pdf",width=10,height=12)
plot(pt6_heatmaps[[2]])
dev.off()