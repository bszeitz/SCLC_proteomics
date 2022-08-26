##############################
#
# Utility functions
#
# Written by Beáta Szeitz
# Project: SCLC proteomics
#
##############################



##############################
# Read in txt files that are in Perseus-like format
# 
# Arguments for this function:
# - raw: the txt file to be read in (directly exported from Perseus)
# - inj.order: TRUE or FALSE. If samples are ordered based on the 
#   injection order, then you can set this to TRUE. Then it will create a new  
#   column called "MSorder" in the sample annotation table.
#
# This function returns a list of tables (protein expression,  
# sample annotation and protein info).
#
# Packages required: none
#
read_in_Perseus <- function (raw, inj.order) {
  raw = read.delim(raw, 
                   stringsAsFactors = FALSE, header=TRUE); 
  if("Contaminant" %in% colnames(raw)){
    raw <- subset(raw, raw[,"Contaminant"] != "TRUE")
  }
  raw2 <- cbind(as.data.frame(matrix(nrow=nrow(raw),ncol=1)),raw)
  colnames(raw2)[1] <- "Annot"
  for (i in 1:nrow(raw2)){
    if (grepl("#", raw[i,1])) {
      x1 = strsplit(raw2[i,2],split= "}")
      x2 = strsplit(x1[[1]][1],split= "[\\{\\]")
      raw2[i,1] <- x2[[1]][2]
      raw2[i,2] <- x1[[1]][2]
    }
  }
  for (i in 1:nrow(raw2)){
    if (grepl("C:", raw2[i,1])) {
      x1 = strsplit(raw2[i,1],split= "C:")
      raw2[i,1] <- x1[[1]][2]
    }
  }
  
  annotations <- subset(raw2, nchar(raw2[,1]) > 0)
  ann <- annotations[,-1]
  rownames(ann) <- annotations[,1]
  ann <- as.data.frame(t(ann))
  ann <- subset(ann, ann$Type=="E")
  ann[,1] <- row.names(ann)
  colnames(ann)[1] <- "MainID"
  
  expressi <- subset(raw2, is.na(nchar(raw2[,1])))
  expressi <- expressi[,-1]
  traw = as.data.frame(t(raw2))
  infos <- subset(traw[-1,], traw[-1,1] != "E")
  infos <- as.data.frame(t(infos))
  infos <- infos[infos$Accession !="" & infos$Accession !="T",]
  row.names(infos) <- infos$Accession
  
  main_table <- expressi[,!names(expressi) %in% c(colnames(infos))]
  row.names(main_table) <- expressi$Accession
  main_table <- as.data.frame(apply(main_table, c(1,2), as.numeric))
  
  options(contrasts = c("contr.sum", "contr.poly"))
  ann <- ann[colnames(main_table),]
  if (inj.order){
    ann$MSorder <- seq(1, nrow(ann), 1)
  }
  
  return(created <- list(main_table,ann,infos))
}



##############################
# Prepare table in long format for ggplot
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
#
# This function returns the protein expression table in a long
# format with columns Sample, Protein and Intensity.
#
# Packages required: reshape2, varhandle
#
prepare_for_ggplot <- function (df) {
  df <- t(as.matrix(df))
  df <- cbind(row.names(df),df)
  df <- reshape2::melt(df, id.vars="V1")
  df <- varhandle::unfactor(df)
  df$value <- suppressWarnings(as.numeric(as.character(df$value)))
  df <- df[df[,2] !="",]
  colnames(df) <- c("Sample", "Protein","Intensity")
  return(df)
}



##############################
# Visualize protein intensity distributions across samples
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
# - project and description: the title and subtitle of the plot
#
# This function returns a series of boxplots to visualize the
# protein intensity distributions in each sample.
#
# Packages required: ggplot2
# 
plot_sample_histograms <- function(df, project, description){
  sampleorder <- colnames(df)
  df <- prepare_for_ggplot(df)
  df.expression.gg_noNA <- na.omit(df)
  df.expression.gg_noNA$Sample <- factor(df.expression.gg_noNA$Sample, levels = sampleorder)
  median.value <- median(df.expression.gg_noNA$Intensity)
  ggplot(df.expression.gg_noNA, aes(x=Sample, y=Intensity))+
    geom_violin() + geom_hline(yintercept = median.value)+
    geom_boxplot(width=0.1) +
    stat_summary(fun.y=median, geom="point", size=2, color="red") +
    labs(title=project, subtitle= description, 
         y="Protein Intensities", x="Samples")+ 
    theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(text= element_text(size=12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.0),
          legend.position = "none")
  
  
}



##############################
# Perform median normalization (centering around the global median)
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows

# This function returns the a new protein expression table where 
# samples are median-normalized (centered around the global median)
#
# Packages required: none
#
normalize_median <- function (df) {
  df2 <- df
  med <- median(as.matrix(df), na.rm=T)
  for (i in 1:ncol(df)){
    df2[,i] <- (df[,i] - median(df[,i], na.rm = T) + med)
  }
  return(df2)
}



##############################
# Short function to create the column annotations for heatmaps
# 
# Arguments for this function:
# - Annotation: the sample annotation table
# - columnnames: the annotation table's columns to be included
# - colorlist: the list of colors
#
# This function returns the column heatmap annotation.
#
# Packages required: ComplexHeatmap
# 
create_heatmapannot <- function(Annotation, columnnames, colorlist){
  column.ha = HeatmapAnnotation(df =Annotation[,columnnames],
                                which="col",
                                col=colorlist,
                                annotation_name_side = "right",
                                gp = gpar(col = "grey"),
                                show_legend = TRUE,
                                show_annotation_name = TRUE)
  return(column.ha)
}



##############################
# Take the median of repeated measurements
# 
# Arguments for this function:
# - df: the protein expression table
# - slist: the list of summarizing names for repeated measurements (e.g. for 
#   SCLC01_DIA_R1 and SCLC01_DIA_R2, the summarizing name is "SCLC01")
# - sepa: the separator with which repeated measurements are
#   distinguished from each other (e.g. for SCLC01_DIA_R1 and SCLC01_DIA_R2,
#   the separator is "_")
#
# This function returns a new dataframe with the median values of
# repeated measurements (in the project, it is used for summarizing measurements
# in the same batches).
#
# Packages required: stats
# 
median_of_repeated_measurements <- function (df, slist, sepa) {
  df2 <- as.data.frame(matrix(nrow=nrow(df), ncol=length(slist)))
  row.names(df2) <- row.names(df)
  for (i in 1:length(slist)){
    name <- paste(slist[i],sepa, sep="")
    collist <- grep(name, colnames(df))
    for (j in 1:nrow(df)) {
      df2[j,i] <- stats::median(as.matrix(df[j,c(collist)]), na.rm = T)
      colnames(df2)[i] <- slist[i]
    }
  }
  return(df2)
}



##############################
# Take the mean of repeated measurements
# 
# Arguments for this function:
# - df: the protein expression table
# - slist: the list of summarizing names for repeated measurements (e.g. for 
#   SCLC01_Batch1 and SCLC01_Batch2, the summarizing name is "SCLC01")
# - sepa: the separator with which repeated measurements are
#   distinguished from each other (e.g. for SCLC01_Batch1 and SCLC01_Batch2,
#   the separator is "_")
#
# This function returns a new dataframe with the mean values of
# repeated measurements (in the project, it is used for summarizing measurements 
# in different batches).
#
# Packages required: none
# 
mean_of_repeated_measurements <- function (df, slist, sepa) {
  df2 <- as.data.frame(matrix(nrow=nrow(df), ncol=length(slist)))
  row.names(df2) <- row.names(df)
  for (i in 1:length(slist)){
    name <- paste(slist[i],sepa, sep="")
    collist <- grep(name, colnames(df))
    for (j in 1:nrow(df)) {
      df2[j,i] <- base::mean(as.matrix(df[j,c(collist)]), na.rm = T)
      colnames(df2)[i] <- slist[i]
    }
  }
  return(df2)
}



##############################
# Filter for valid values based on percentage
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
# - percentage: the minimum percentage of valid values required for 
#   each protein
#
# This function returns a protein expression table that is filtered
# for valid values (> or = the given percentage).
#
# Packages required: none
# 
filter_missingvalues <- function (df, percentage) {
  perc <- ceiling(ncol(df)*(percentage/100))
  vv <- rowSums(!is.na(df))
  to.be.kept <- unlist(lapply(vv, function(x){ 
    if (x < perc) { y <- F}
    else y <- T}))
  return(df[to.be.kept,])
}



##############################
# Batch correction by linear regression
# 
# Arguments for this function:
# - df: a protein expression table where samples are in columns and
#   proteins are in rows
# - annot: a sample annotation table with 2 columns: MainID (= sample names) 
#   and MSbatch (= batches)
#
# This function returns a protein expression table that is corrected
# for batch effect (in our case, there were 2 batches).
#
# Packages required: stats
#
batch_correction_LM <- function(df, annot){
  df <- as.data.frame(df)
  Epxr_batchcorr <- df
  for (i in 1:nrow(df)) {
    prot.data <- as.data.frame(t(df[i,]))
    prot.data$MainID <- row.names(prot.data)
    data <- merge(prot.data, annot,by="MainID")
    row.names(data) <- data$MainID
    data <- data[,-grep("MainID",colnames(data))] # remove the "sample" column
    #VVs.batch1 <- sum(!is.na(data[data[,2] == 1,1]))
    #VVs.batch2 <- sum(!is.na(data[data[,2] == 2,1]))
    #if (VVs.batch1 == 0 | VVs.batch2 ==0){
    #  next
    #}
    data[,2] <- ifelse(data[,2] ==1, 0, 1)
    lm_res = lm(data[,1] ~ data[,2],data=data, na.action=na.exclude)
    newval <- as.data.frame(lm_res[["residuals"]] + 
                              lm_res[["coefficients"]][["(Intercept)"]])
    newdata <- merge(data, newval, by='row.names', all.x=T)
    Epxr_batchcorr[i,] <- t(newdata[,4])
    row.names(Epxr_batchcorr)[i] <- names(data)[1]
  }
  df <- df[,order(colnames(df))]
  colnames(Epxr_batchcorr) <- colnames(df)
  return(Epxr_batchcorr)
}



##############################
# Perform Principal Component Variance Analysis (PVCA) and plot results
# 
# Arguments for this function:
# - m: a protein expression table where samples are in columns and
#   proteins are in rows
# - annotation: a sample annotation table
# - ID.column: the column name that entails the sample names
# - columns.to.include: the annotation names to be included in the PVCA
# - pct.threshold: argument of pvca() (the % of the minimum amount of the 
#   variabilities that the selected principal components need to explain), 
#   default is 0.6
# - color: the color of the barplots
# - title: the figure title
# - figure: draw plot? (TRUE or FALSE). If FALSE, then only PVCA results are
#   returned.
#
# This function returns PVCA results and optionally makes a barplot figure of
# the results.
#
# Packages required: PVCA, BioBase
#
plot_PVCA <- function(m, annotation, ID.column, columns.to.include, pct.threshold = 0.6, color="blue", title, figure = T){
  m <- t(scale(t(m)))
  row.names(annotation) <- annotation[,ID.column]
  phenoData.original <- annotation[colnames(m),c(ID.column, columns.to.include)]
  colnames(phenoData.original) <- sapply(colnames(phenoData.original), function(x){gsub(" ", ".", x,fixed=T)})
  phenoData <- new("AnnotatedDataFrame", data=phenoData.original)
  m.expressionSet <-  ExpressionSet(assayData=as.matrix(na.omit(m)),
                                     phenoData=phenoData)
  
  batch.factors <- sapply(columns.to.include, function(x){gsub(" ", ".", x,fixed=T)})
  pvcaObj <- pvcaBatchAssess (m.expressionSet, batch.factors, pct.threshold)
  pvcaObj.DF <- data.frame(label = pvcaObj$label, dat = as.numeric(pvcaObj$dat),
                           row.names = pvcaObj$label)
  effect.order <- pvcaObj.DF[order(pvcaObj.DF$dat, decreasing = T),"label"]
  pvcaObj.DF <- pvcaObj.DF[c(effect.order[effect.order!="resid"], "resid"),]
  if (figure){
    return(ggplot(pvcaObj.DF, aes(x=label, y=dat))+ ggtitle(title)+ 
             geom_bar(stat="identity", fill=color, colour=color) + 
             geom_text(aes(label=round(dat,3), y=dat+0.04), size=5) +
             scale_x_discrete(limits=pvcaObj.DF$label) + 
             scale_y_continuous(limits = c(0,1)) +
             labs(x= "Effects", y= "Weighted avr. proportion var.") +theme_bw() +
             theme(plot.background = element_blank() ,panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank() ,panel.border = element_blank(), panel.background = element_blank())+
             theme(axis.line = element_line(color = 'black'),
                   axis.title.x = element_text(size = 15, vjust=-0.5),
                   axis.title.y = element_text(size = 15, vjust= 1.0),
                   axis.text = element_text(size = 12),
                   axis.text.x = element_text(angle = 90, vjust= 0.5, hjust=1)))
  } else {
    return(pvcaObj.DF)
  }
  
}



##############################
# Export txt files where "Accession" column shows the row names
# 
# Arguments for this function:
# - df: the dataframe to be exported
# - filename: the name of the new file (no .txt needed at the end)
# - row.name.needed: TRUE or FALSE. If TRUE, then a new column "Accession" 
#   will be created, which entails the row names of the original dataframe
#
# This function exports a dataframe with the specified name.
#
# Packages required: none
#
export_with_rowname <- function(df, filename, row.name.needed=T){
  if (row.name.needed){
    df <- cbind(row.names(df), df)
    colnames(df)[1] <- "Accession"
  }
  write.table(df, paste0(filename,".txt"), row.names = F, na = "", quote = F, sep="\t")
}



##############################
# Read in previously generated txt files (with export_with_rowname function)
# where "Accession" column shows the row names
# 
# Arguments for this function:
# - filename: the txt file to be read in (no .txt needed at the end)
# - include.rowname: TRUE or FALSE. If TRUE, then the column "Accession" will
#   be kept in the data frame, otherwise it is deleted
#
# This function returns a dataframe where the "Accession" column is used for  
# row names.
#
# Packages required: none
#
read_in_with_rowname <- function(filename, include.rowname){
  raw <- read.delim(paste0(filename,".txt"), sep="\t", check.names = F)
  row.names(raw) <- raw$Accession
  if (!include.rowname){
    raw$Accession <- NULL
  }
  return(raw)
}



##############################
# Export expression table in gct format
# 
# Arguments for this function:
# - df: the dataframe to be exported
# - ID.col: the column with row names
# - exportname: the name of the new gct file
#
# This function exports an expression table in gct format (required for ssGSEA).
#
# Packages required: none
#
save_as_gct <- function(df, ID.col,exportname){
  decoy <- matrix(nrow=3, ncol=ncol(df)+1)
  decoy[1,1] <- "#1.2"
  decoy[2,1:2] <- c(nrow(df), ncol(df)-1)
  decoy[3,] <- c("NAME","Description", colnames(df)[-grep(ID.col, colnames(df))])
  
  df <- cbind(df[,ID.col], df[,ID.col], df[,-grep(ID.col, colnames(df))])
  colnames(df)[1:2] <- c("NAME","Description")
  colnames(decoy) <- colnames(df)
  
  decoy.df <- rbind(decoy, df)
  
  write.table(decoy.df, exportname, row.names = F, col.names = F, sep = "\t", quote = F, na = "")
}



##############################
# Perform ANOVA and pairwise Tukey HSD test
# 
# Arguments for this function:
# - m: a protein expression table where samples are in columns and
#   proteins are in rows
# - annotation: the annotation table of the samples
# - colname_for_factor: the column name in the annotation table
#   that should be used as the independent variable
# - levelOrder: how should the levels of the independent variable
#   be ordered
# - analysisName: if the columns in the result table should be 
#   appended by a certain string
# - stars: whether the result table should contain columns that
#   indicate significance with stars (*** p<0.005, ** p<0.01,
#   * p<0.05, . p<0.10)
#
# This function returns a result table that contains:
# - Protein accession numbers
# - ANOVA test F, p-value, Benjamini-Hochberg adjusted p-value
# - Pairwise Tukey HSD test p-values
# - Pairwise Log2FC values
# - 95% confidence intervals for the pairwise log2FC values
# - Protein ranks for all pairwise comparisons (rank = -log10p * log2FC)
#
# Packages required: varhandle, DTK, reshape2, car
#
ANOVAandTK_padj_rank <- function(m, annotation, colname_for_factor, levelOrder, analysisName="", stars=F) {
  annotation = as.matrix(cbind(row.names(annotation), annotation[,colname_for_factor]))
  colnames(annotation) <- c("Sample", colname_for_factor)
  annotation <- as.data.frame(annotation)
  if (is.factor(annotation[,1])){
    annotation[,1] <- varhandle::unfactor(annotation[,1])
  }
  if (!is.factor(annotation[,2])){
    annotation[,2] <- factor(annotation[,2], levels=levelOrder)
  } else {
    annotation[,2] <- factor(varhandle::unfactor(annotation[,2]), levels=levelOrder) 
  }
  annotation <- na.omit(annotation)
  
  m <- m[,annotation$Sample]
  
  p.v <- matrix(nrow=nrow(m), ncol=1)
  Fvalues <- matrix(nrow=nrow(m), ncol=1)
  n <- as.numeric(length(levels(annotation[,2])))
  k <- 2 # number of comparations
  comb <- factorial(n)/(factorial(k)*factorial(n-k))
  
  #names
  m.noNA <- na.omit(m)
  prot.data <- as.data.frame(t(m.noNA[1,]))
  prot.data$'Sample' <- row.names(prot.data)
  data <- merge(prot.data,annotation,by="Sample")
  row.names(data) <- data$Sample
  data <- data[,which(colnames(data) !="Sample")] # remove the "sample" column
  data <- na.omit(data)
  tukey <- DTK::TK.test(x=data[,1],f=data[,2], a= 0.05)
  tukey <- reshape2::melt(tukey[["f"]])
  tukey.pv0 <- subset(tukey, Var2 =="p adj")
  comparisons <- tukey.pv0[,1]
  
  tukey.pv <- matrix(nrow=nrow(m), ncol=comb)
  colnames(tukey.pv) <- paste("tukey(",tukey.pv0[,1],")",sep="")
  tukey.diff <- matrix(nrow=nrow(m), ncol=comb)
  colnames(tukey.diff) <- paste("Log2FC(",tukey.pv0[,1],")",sep="")
  tukey.confint <- matrix(nrow=nrow(m), ncol=comb)
  colnames(tukey.confint) <- paste("Log2FC.95%.CI(",comparisons,")",sep="")
  ranking <- matrix(nrow=nrow(m), ncol=comb)
  colnames(ranking) <- paste("Rank(",tukey.pv0[,1],")",sep="")
  
  for (i in 1:nrow(m)) {
    prot.data <- as.data.frame(t(m[i,]))
    prot.data$'Sample' <- row.names(prot.data)
    data <- merge(prot.data,annotation,by="Sample")
    row.names(data) <- data$Sample
    data <- data[,which(colnames(data) !="Sample")] # remove the "sample" column
    data <- na.omit(data)
    
    lm_res = lm(data[,1] ~ data[,2],data=data, na.action=na.exclude)
    ANOVAres <- car::Anova(lm_res, type=2)
    Fvalues[i,1] <- paste0(round(ANOVAres$`F value`[1], 4), " (df=",ANOVAres$Df[1],")" )
    p.v[i,1] <- ANOVAres$`Pr(>F)`[1]
    tukey <- DTK::TK.test(x=data[,1],f=data[,2], a= 0.05)
    tukey <- reshape2::melt(tukey[["f"]])
    
    # tukey.PV
    tukey.pv0 <- subset(tukey, Var2 =="p adj")
    tukey.pv1 <- as.data.frame(t(tukey.pv0[,3]))
    colnames(tukey.pv1) <- paste("tukey(",tukey.pv0[,1],")",sep="")
    for (K in 1:ncol(tukey.pv1)){
      tukey.pv[i,colnames(tukey.pv1)[K]] <- tukey.pv1[1,K]
    }
    
    # tukey.diff
    tukey.diff0 <- subset(tukey, Var2 =="diff")
    tukey.diff1 <- as.data.frame(t(tukey.diff0[,3]))
    colnames(tukey.diff1) <- paste("Log2FC(",tukey.diff0[,1],")",sep="")
    for (K in 1:ncol(tukey.diff1)){
      tukey.diff[i,colnames(tukey.diff1)[K]] <- tukey.diff1[1,K]
    }
    
    # tukey.CI
    tukey.CI0 <- subset(tukey, Var2 =="lwr" | Var2 =="upr")
    tukey.CI0$Var1 <- varhandle::unfactor(tukey.CI0$Var1)
    comps <- unique(tukey.CI0$Var1)
    tukey.CI <- tukey.CI0[1:length(comps),]
    
    for (K in 1:length(comps)){
      uprCI <- round(tukey.CI0[tukey.CI0$Var1 == comps[K] & tukey.CI0$Var2 =="upr","value"],4)
      lwrCI <- round(tukey.CI0[tukey.CI0$Var1 == comps[K] & tukey.CI0$Var2 =="lwr","value"],4)
      tukey.CI[K,3] <- paste0("[",lwrCI," , ", uprCI,"]")
    }
    tukey.CI2 <- as.data.frame(t(tukey.CI[,3]))
    colnames(tukey.CI2) <- paste("Log2FC.95%.CI(",tukey.CI[,1],")",sep="")
    for (K in 1:ncol(tukey.CI2)){
      tukey.confint[i,colnames(tukey.CI2)[K]] <- tukey.CI2[1,K]
    }
    
    # ranking
    ranking.1 <- -log10(tukey.pv1)*tukey.diff1
    colnames(ranking.1) <- unlist(lapply(colnames(ranking.1), function(x){gsub("tukey","Rank", x)}))
    for (K in 1:ncol(ranking.1)){
      ranking[i,colnames(ranking.1)[K]] <- ranking.1[1,K]
    }
    
  }
  row.names(p.v) <- row.names(m)
  p.v <- as.data.frame(p.v)
  colnames(p.v) <- c(paste("p.v", colnames(annotation)[2]))
  row.names(tukey.pv) <- row.names(m)
  row.names(tukey.diff) <- row.names(m)
  row.names(ranking) <- row.names(m)
  adjp <- p.adjust(as.vector(p.v[,1]), "fdr")
  adjp <- as.data.frame(as.numeric(as.character(adjp)))
  p.v_ANOVA <- as.data.frame(cbind(p.v,adjp,tukey.pv))
  colnames(p.v_ANOVA)[2] <- c(paste("p.adj", colnames(annotation)[2])) 
  if (stars){
    p.v_star <- p.v_ANOVA
    for (C in 1:ncol(p.v_star)){
      p.v_star[,C] <- cut(as.numeric(p.v_ANOVA[,C]), breaks=c(-Inf,0.005, 0.01, 0.05, 0.1, Inf), 
                          label=c("***", "**", "*", ".",""))
    }
    colnames(p.v_star) <- paste("*",colnames(p.v_star) , sep="_")
    table_results <- cbind(Fvalues, p.v_ANOVA,tukey.diff,ranking, tukey.confint, p.v_star)
  } else {
    table_results <- cbind(Fvalues, p.v_ANOVA,tukey.diff,ranking, tukey.confint)
  }
  if (analysisName !=""){
    colnames(table_results) <- paste(colnames(table_results) ,analysisName, sep="_")
  }
  table_results <- as.data.frame(cbind(row.names(table_results),table_results))
  colnames(table_results)[1] <- "Accession"
  return(table_results)
}






##############################
# Perform Kruskal-Wallis and pairwise Wilcox tests
# 
# Arguments for this function:
# - m: an expression table where samples are in columns and
#   features are in rows
# - annotation: the annotation table of the samples
# - colname_for_factor: the column name in the annotation table
#   that should be used as the independent variable
# - levelOrder: how should the levels of the independent variable
#   be ordered
# - analysisName: if the columns in the result table should be 
#   appended by a certain string
#
# This function returns a result table that contains:
# - Feature names
# - Kruskal-Wallis test p-value, Benjamini-Hochberg adjusted p-value
# - Pairwise Wilcox test p-values
# - Pairwise Log2FC values
#
# Packages required: varhandle, DTK, reshape2
#
Kruskal_posthoc_adj <- function(m, annotation, colname_for_factor, levelOrder, analysisName="", min.perc.in.group) {
  annotation = as.matrix(cbind(row.names(annotation), annotation[,colname_for_factor]))
  colnames(annotation) <- c("Sample", colname_for_factor)
  annotation <- as.data.frame(annotation)
  if (is.factor(annotation[,1])){
    annotation[,1] <- varhandle::unfactor(annotation[,1])
  }
  if (!is.factor(annotation[,2])){
    annotation[,2] <- factor(annotation[,2], levels=levelOrder)
  } else {
    annotation[,2] <- factor(varhandle::unfactor(annotation[,2]), levels=levelOrder) 
  }
  annotation <- na.omit(annotation)
  
  m <- m[,annotation$Sample]
  
  
  min.sample.sizes <- vector(length=length(levels(annotation[,2])))
  names(min.sample.sizes) <- levels(annotation[,2])
  for (i in 1:length(min.sample.sizes)){
    min.sample.sizes[i] <- nrow(annotation[annotation[,2]==names(min.sample.sizes)[i],]) * min.perc.in.group
  }
  
  p.v <- matrix(nrow=nrow(m), ncol=1)
  n <- as.numeric(length(levels(annotation[,2])))
  k <- 2 
  comb <- factorial(n)/(factorial(k)*factorial(n-k))
  
  
  #names
  m.noNA <- as.matrix(na.omit(m))
  prot.data <- as.data.frame(m.noNA[1,])
  prot.data$'Sample' <- row.names(prot.data)
  data <- merge(prot.data,annotation,by="Sample")
  row.names(data) <- data$Sample
  data <- data[,which(colnames(data) !="Sample")] # remove the "sample" column
  data <- na.omit(data)
  tukey <- DTK::TK.test(x=data[,1],f=data[,2], a= 0.05)
  tukey <- reshape2::melt(tukey[["f"]])
  tukey.pv0 <- subset(tukey, Var2 =="p adj")
  
  comparisons <- tukey.pv0[,1]
  
  tukey.pv <- matrix(nrow=nrow(m), ncol=length(comparisons))
  colnames(tukey.pv) <- paste("wilcox(",tukey.pv0[,1],")",sep="")
  tukey.diff <- matrix(nrow=nrow(m), ncol=length(comparisons))
  colnames(tukey.diff) <- paste("Log2FC(",tukey.pv0[,1],")",sep="")
  
  m <- as.matrix(m)
  
  for (i in 1:nrow(m)) {
    #print(i)
    m.naomit <- na.omit(m[i,])
    if (length(m.naomit) < 5){ next}
    if (sd(m.naomit) <0.01){ next}
    prot.data <- as.data.frame(m[i,])
    prot.data$'Sample' <- row.names(prot.data)
    data <- merge(prot.data,annotation,by="Sample")
    row.names(data) <- data$Sample
    data <- data[,which(colnames(data) !="Sample")] # remove the "sample" column
    data <- na.omit(data)
    groups.to.include <- vector()
    
    for (N in 1:length(min.sample.sizes)){
      groups.to.include <- c(groups.to.include, ifelse(nrow(data[data[,2]==names(min.sample.sizes)[N],]) > min.sample.sizes[N], 
                                                       names(min.sample.sizes)[N], NA))
    }
    groups.to.include <- na.omit(groups.to.include)
    
    if (length(groups.to.include)==0 | length(groups.to.include)==1){ next}
    data <- data[which(data[,2] %in% groups.to.include),]
    
    p.v[i,1] <- kruskal.test(data[,1] ~ data[,2],data=data, na.action=na.exclude)[["p.value"]]
    tukey <- pairwise.wilcox.test(data[,1], data[,2],
                                  p.adjust.method = "BH", exact=F)
    tukey <- na.omit(reshape2::melt(tukey[["p.value"]]))
    tukey$Var1 <- paste(tukey$Var1, tukey$Var2, sep="-")
    tukey$Var2 <- NULL
    
    # tukey.PV
    tukey.pv1 <- as.data.frame(t(tukey[,2]))
    colnames(tukey.pv1) <- paste("wilcox(",tukey[,1],")",sep="")
    for (K in 1:ncol(tukey.pv1)){
      tukey.pv[i,colnames(tukey.pv1)[K]] <- tukey.pv1[1,K]
    }
    # tukey.diff
    tukey <- TukeyHSD(aov(data[,1] ~ data[,2],data=data, na.action=na.exclude))$`data[, 2]`
    tukey <- reshape2::melt(tukey)
    tukey.diff0 <- subset(tukey, Var2 =="diff")
    tukey.diff1 <- as.data.frame(t(tukey.diff0[,3]))
    colnames(tukey.diff1) <- paste("Log2FC(",tukey.diff0[,1],")",sep="")
    for (K in 1:ncol(tukey.diff1)){
      tukey.diff[i,colnames(tukey.diff1)[K]] <- tukey.diff1[1,K]
    }
  }
  row.names(p.v) <- row.names(m)
  p.v <- as.data.frame(p.v)
  colnames(p.v) <- c(paste("p.v", colnames(annotation)[2]))
  row.names(tukey.pv) <- row.names(m)
  row.names(tukey.diff) <- row.names(m)
  adjp <- p.adjust(as.vector(p.v[,1]), "fdr")
  adjp <- as.data.frame(as.numeric(as.character(adjp)))
  p.v_ANOVA <- as.data.frame(cbind(p.v,adjp,tukey.pv))
  colnames(p.v_ANOVA)[2] <- c(paste("p.adj", colnames(annotation)[2])) 
  p.v_ANOVA.star <- p.v_ANOVA
  for (C in 1:ncol(p.v_ANOVA.star)){
    p.v_ANOVA.star[,C] <- cut(p.v_ANOVA[,C], breaks=c(-Inf,0.005, 0.01, 0.05, 0.1, Inf), 
                              label=c("***", "**", "*", ".",""))
  }
  colnames(p.v_ANOVA.star) <- paste("*",colnames(p.v_ANOVA.star) , sep="_")
  table_results <- cbind(p.v_ANOVA,tukey.diff,p.v_ANOVA.star)
  colnames(table_results) <- paste(colnames(table_results) ,analysisName, sep="_")
  table_results <- as.data.frame(table_results)
  table_results$Accession <- row.names(m)
  return(table_results)
}


##############################
# Perform Shapiro-Wilk Normality Test
# 
# Arguments for this function:
# - m: an expression table where samples are in columns and
#   features are in rows
#
# This function returns a result table (dataframe) that contains:
# - Accession (i.e., the row names)
# - Shapiro-Wilk normality test p-value
#
# Packages required: stats
#
SW_normality_test <- function(m){
  #x <- GDSC1.SCLC.Resp.F[1,]
  SW.pvalues <- apply(m, 1, function(x){
    res <- shapiro.test(as.numeric(x))
    return(res$p.value)
  })
  return(data.frame(Accession = row.names(m),
                    Shapiro.Wilk.pv = SW.pvalues))
}

##############################
# Perform min-max scaling
# 
# Arguments for this function:
# - values: a vector with numbers that should be scaled
#
# This function returns a vector of scaled values.
#
# Packages required: none
#
min_max_scaling <- function(values){
  values <- as.numeric(values)
  min.value <- min(values, na.rm=T)
  max.value <- max(values, na.rm=T)
  scaled.values <- (values - min.value) / (max.value - min.value)
  return(scaled.values)
}
