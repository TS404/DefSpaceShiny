
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

# Import necessary files ----------------------------------------------------------------------
library(shiny)
library("rgl")
library("ggplot2")
library("Biostrings")
library("DECIPHER")
library("mclust")

folder    <- "C:\\Users\\T\\OneDrive\\1-Scripts\\GitHub\\DefSpace"
setwd(folder)

SAPCA.cis <- readRDS("data\\CisDef.reference.PCA.RDS")
SAPCA.tra <- readRDS("data\\TransDef.reference.PCA.RDS")

view.cis <- readRDS("data\\CisDef.viewangle.RDS")
view.tra <- readRDS("data\\TransDef.viewangle.RDS")

BLOSUM40 <- readRDS("data\\BLOSUM.RDS")

clusters.cis=c("'extreme' plant antimicrobial defensins",  #1
               "mostly plant antimicrobial defensins",     #2
               "proteins with a mixture of functions from across the eukarya", #3
               "plant signalling proteins",                #4
               "plant histidine-rich defensins",           #5
               "arthropod antimicrobial defensins",        #6
               "arthropod alpha neurotoxins")              #7

clusters.tra=c("theta defensins", #1
               "alpha defensins", #2
               "beta defensins",  #3
               "big defensins")   #4

# Outputs -------------------------------------------------------------------------

shinyServer(function(input, output) {
  
  # Calculations when button pressed
  DATA <- eventReactive(input$button,{
    
    newseq.cis <- seq.MSA.add(SAPCA.cis,input$query_sequence,"cis-Defensins")
    newseq.tra <- seq.MSA.add(SAPCA.tra,input$query_sequence,"trans-Defensins")
    
    if(input$query_type=="cis-Defensin" |
       input$query_type=="unknown" & newseq.cis$aln.hit.score >= newseq.tra$aln.hit.score){
        match        <- "cis-Defensin"
        SAPCA.match  <- SAPCA.cis
        newseq.match <- newseq.cis
        view         <- view.cis
        plotPCs      <- c(1,2,3)
    }else if(input$query_type=="trans-Defensin" |
             input$query_type=="unknown" & newseq.cis$aln.hit.score <= newseq.tra$aln.hit.score){
        match        <- "trans-Defensin"
        SAPCA.match  <- SAPCA.tra
        newseq.match <- newseq.tra
        view         <- view.tra
        plotPCs      <- c(1,2,4)
    }

    temp.seq  <- gsub("subject: |[[]1]| ","",capture.output(newseq.match$seq.unalignable)[3])
    newseq.r  <- seq.rotate    (SAPCA.match, newseq.match)
    newseq.c  <- seq.clust.add (SAPCA.match, newseq.r)
    SAPCA.add <- seq.SAPCA.add (SAPCA.match, newseq.r, newseq.c)
    
    if(match=="cis-Defensin"){
       clust.name   <- clusters.cis[SAPCA.add$seq.space.clusters$classification[1]]
    }else if(match=="trans-Defensin"){
       clust.name   <- clusters.tra[SAPCA.add$seq.space.clusters$classification[1]]
    }
    
    # data variables list
    DATA <- list(match        = match,
                 SAPCA.match  = SAPCA.match,
                 newseq.match = newseq.match,
                 view         = view,
                 plotPCs      = plotPCs,
                 clust.name   = clust.name,
                 temp.seq     = temp.seq,
                 SAPCA.add    = SAPCA.add)
    DATA
  })
  
  # Main seqspace plot
  output$mainplot <- renderPlot({
  
    # plot(DATA()$plotPCs)
    if(DATA()$match=="cis-Defensin"){
      colours<-palette(c("blue",            #1 Plant extreme
                         "darkolivegreen4", #2 Plant main
                         "grey",            #3 Intermed
                         "purple1",         #4 Plant sex
                         "orange",          #5 Plant his
                         "maroon",          #6 Invert
                         "red"))            #7 Tox
    }
    if(DATA()$match=="trans-Defensin"){
      colours<-palette(c("blue",     #1 Theta
                         "red",      #2 Aalpha
                         "orange",   #3 Beta
                         "purple"))  #4
    }

    plot_3Dclusters(DATA()$SAPCA.add,
                    plotPCs = DATA()$plotPCs,
                    labels = "query",
                    radius = c(2,rep(0.3,nrow(DATA()$SAPCA.add$numerical.alignment$MSA)-1)))
 
  })
  
  #####TEMP
  output$temp <- renderText({print(   DATA()$temp.seq  )})
  ####
  
  # Report on match
  # Superfamily type and best hit
  output$report <- renderText({
   
    if(input$query_type=="unknown"){
      if(quantile(DATA()$newseq.match$aln.all.score, 0.95)>=0.15){
        print(paste0("The submitted query sequence appears to be a ",
                     DATA()$match,
                     ". Its similarity to the nearest sequence is ",
                     percent(DATA()$newseq.match$aln.hit.score),
                     "."))
      }else{
        print(paste0("Sequence may not be a defensin. Its similarity to any known defensin is only ",
                     percent(DATA()$newseq.match$aln.hit.score),
                     "."))
      }
    }
  })
  
  # Matching residues
  output$report2 <- renderText({
    
    exceptions <- if(sum(strsplit(DATA()$temp.seq,"")[[1]]=="-")!=0){
                     paste0(", except for ",
                            sum(strsplit(DATA()$temp.seq,"")[[1]]=="-"),
                            ". The residues that were taken into account in calculating its sequence space position were therefore: '",
                            DATA()$temp.seq,
                            "'")
                   }
    
    paste0("All residues of the query sequence were alignable to the existing ",
           DATA()$match,
           " MSA",
           exceptions,
           ".")
  })
  
  # Cluster
  output$report3 <- renderText({
    
    paste0("The sequence falls within cluster ",
           DATA()$SAPCA.add$seq.space.clusters$classification[1],
           ", which contains ",
           DATA()$clust.name)
  })
  
  # Nearest neighbours fasta alignment
  output$fasta <- renderUI({
    
    HTML(
      as.fasta(DATA()$SAPCA.add$numerical.alignment$MSA[rownames(closest(DATA()$SAPCA.add,n = input$return_nearest,"query")),],
               decolgap = TRUE,
               print    = TRUE)
    )
  })
  
})


# Functions ------------------------------------------------------------------------------

numericise_MSA <- function(MSA,
                           res.prop,
                           cys){
  seq.names <- rownames(MSA)
  aln.len   <- ncol(MSA)
  res.props <- colnames(res.prop)
  res.avail <- row.names(res.prop)
  # Numericise MSA based on res.prop
  MSA.num.tall <- res.prop[t(MSA),]
  # Name data types
  rownames(MSA.num.tall) <- NULL
  sequence     <- rep(x = seq.names, each  = aln.len)
  residue      <- rep(x = 1:aln.len, times = length(seq.names))
  MSA.num.tall <- cbind(sequence, residue, MSA.num.tall)
  # Stack data into list of matrices
  MSA.num.stack <- NULL
  for (x in 1:length(res.props)) {
    col.names <- paste(1:aln.len,
                       rep(res.props[x],aln.len),
                       sep = ".")
    MSA.num.stack[[res.props[x]]] <- matrix(MSA.num.tall[,x+2],
                                            ncol     = aln.len,
                                            byrow    = TRUE,
                                            dimnames = list(seq.names,
                                                            col.names))
  }
  # Also reflow into single wide matrix
  MSA.num.wide <- MSA.num.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.num.wide <- cbind(MSA.num.wide, MSA.num.stack[[res.props[x]]])
  }
  
  
  ############################.
  # Scaling by property type #
  ############################.
  
  # Take means and variances of each property type
  prop.means <- NULL
  prop.vars  <- NULL
  for (x in 1:length(res.props)) {
    prop.means[x] <- mean(MSA.num.stack[[x]],na.rm=1)
    prop.vars[x]  <- var(tidyr::gather(data.frame(MSA.num.stack[[x]]))[2],na.rm=1)
  }
  names(prop.means) <- res.props
  names(prop.vars)  <- res.props
  
  # Scale numericised MSA to prop.means and prop.vars
  MSA.scale.stack <- NULL
  for (x in 1:length(res.props)) {
    MSA.scale.stack[[res.props[x]]] <- (MSA.num.stack[[res.props[x]]]- prop.means[x]) /
      sqrt(prop.vars[x])
  }
  
  # Replace gaps (currently "NA") with column average
  # Create na.colmean function
  na.colmean<-function(x){
    x[is.na(x)] <- mean(as.matrix(x),na.rm = 1)
    x
  }
  # For each property of MSA.num.stack, apply na.colmean function to each matrix comlumn
  for (x in 1:length(res.props)) {
    MSA.scale.stack[[x]] <- apply(MSA.scale.stack[[x]],2,na.colmean)
  }
  
  # Also reflow into singe wide matrix for PCA
  MSA.scale.wide <- MSA.scale.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.scale.wide <- cbind(MSA.scale.wide, MSA.scale.stack[[x]])
  }
  
  
  ##################.
  # Alignment list #
  ##################.
  
  numerical.alignment <- list(MSA             = MSA,
                              res.prop        = res.prop,
                              MSA.num.stack   = MSA.num.stack,
                              MSA.num.wide    = MSA.num.wide,
                              MSA.scale.stack = MSA.scale.stack,
                              MSA.scale.wide  = MSA.scale.wide,
                              prop.means      = prop.means,
                              prop.vars       = prop.vars,
                              seq.names       = seq.names,
                              aln.len         = aln.len)
  numerical.alignment
}



closest <- function (SAPCA,
                     sequence,
                     PC = 1:3,
                     n  = 10){
  
  coords     <- SAPCA$seq.space.PCA$coordinates
  centre     <- coords[sequence,PC]
  centre.m   <- matrix(rep(centre,nrow(coords)),
                       nrow  = nrow(coords),
                       byrow = TRUE)
  distances  <- SAPCA$seq.space.PCA$coordinates[,PC]-centre.m
  rootsquare <- sqrt(rowSums(distances^2))
  
  sorted     <- as.matrix(rootsquare[order(rootsquare)])
  colnames(sorted) <- "distance"
  
  head(sorted,n)
}



as.fasta <- function(matrix,degap=FALSE,decolgap=FALSE,write=FALSE,print=FALSE,name=NULL){
  
  # Remove empty columns
  if(decolgap){
    matrix<-matrix[,colMeans(matrix=="-")!=1]
  }
  
  # Convert alignment matrix to list of strings
  names <- paste(">",row.names(matrix),sep="")
  seqs  <- do.call("paste",c(data.frame(matrix),sep=""))
  
  # If just one sequence, this is how to name it
  if(is.null(dim(matrix))){
    names <- ">sequence"
    if(!is.null(name)){
      names <- paste(">",name,sep="")
    }
    seqs  <- paste(matrix,collapse="")
  }
  
  # Degap sequences
  if (degap){
    seqs <- gsub("-","",seqs)
  }
  
  # Interleave names and sequences
  ord1 <- 2*(1:length(names))-1
  ord2 <- 2*(1:length(seqs))
  
  # Output
  if (print==TRUE){
    paste0(c(names,seqs)[order(c(ord1,ord2))], sep = "<br/>")
  }else{
    if (write==FALSE){
      cat(c(names,seqs)[order(c(ord1,ord2))], sep = "\n")
    }
    else{
      if (!grepl(".fa",write,ignore.case=TRUE)){
        write<-paste(write,".fa",sep="")
      }
      cat(c(names,seqs)[order(c(ord1,ord2))], sep = "\n", file = write)
    }
  }
}



as.AAstring<-function(string, degap=FALSE){
  string <- paste(string,collapse="")
  if(degap==TRUE){
    string<-gsub("-","",string)
  }
  output <- Biostrings::AAString(string)
  output
}


as.AAstringSet<-function(MSA, degap=FALSE){
  MSA <- apply(MSA,1,paste,collapse="")
  if(degap==TRUE){
    MSA<-gsub("-","",MSA)
  }
  output <- Biostrings::AAStringSet(MSA)
  output
}



seq.MSA.add <- function(SAPCA,sequence,SAPCAname=NULL){
  MSA   <- SAPCA$numerical.alignment$MSA
  MSA2  <- as.AAstringSet(MSA,degap = TRUE)
  seqs  <- nrow(MSA)
  seq   <- as.AAstring(sequence, degap=FALSE)
  seq.d <- as.AAstring(sequence, degap=TRUE)
  BLOSUM40 <- blosum()
  
  aln.all <- Biostrings::pairwiseAlignment(MSA2,
                                           seq.d,
                                           substitutionMatrix = BLOSUM40,
                                           gapOpening   = 0,
                                           gapExtension = 8,
                                           scoreOnly    = TRUE)
  
  # Max possible similarity score
  aln.limit <- Biostrings::pairwiseAlignment(seq.d,
                                             seq.d,
                                             substitutionMatrix = BLOSUM40,
                                             gapOpening   = 0,
                                             gapExtension = 8,
                                             scoreOnly    = TRUE)
  
  # Similarity score as percentage of max
  aln.hit.score <- max(aln.all)/aln.limit
  
  # The sequence of the best matching member of the database
  aln.hit.num  <- which(aln.all==max(aln.all))[1]
  aln.hit.name <- SAPCA$numerical.alignment$seq.names[aln.hit.num]
  aln.hit.seq  <- paste(as.AAstring(MSA[aln.hit.num,],degap = 1))
  
  # Use "*" to indicate gaps in the best reference sequence (aln.hit)
  aln.hit <- gsub("-","*",as.AAstring(MSA[aln.hit.num,]))
  aln.add <- Biostrings::pairwiseAlignment(aln.hit,
                                           seq.d,
                                           substitutionMatrix = BLOSUM40,
                                           gapOpening   = 0,
                                           gapExtension = 8)
  
  # Has the new sequence introduced exrta gaps into the hit sequence alignement?
  aln.hit.orig         <- as.AAstring(MSA[aln.hit.num,])
  aln.hit.new          <- Biostrings::pattern(aln.add)
  aln.hit.seq.aln.orig <- unlist(strsplit(as.character(aln.hit.orig),""))
  aln.hit.seq.aln.new  <- unlist(strsplit(as.character(aln.hit.new),""))
  
  if(as.character(aln.hit.orig)!=as.character(aln.hit.new)){
    print(paste(sum(aln.hit.seq.aln.new=="-"),
                "residues of the new sequence were not alignable to the",
                SAPCAname,
                "reference MSA so were ignored"))
  }
  
  # Alignment addition as matrix
  aln.add.mat <- rbind(unlist(strsplit(as.character(Biostrings::pattern(aln.add)),"")),
                       unlist(strsplit(as.character(Biostrings::subject(aln.add)),"")))
  
  # Unmathcable resiues removed from aligned sequence
  aln.add2 <- aln.add.mat[2,][aln.add.mat[1,]!="-"]
  aln.add3 <- paste(as.AAstring(aln.add2))
  
  # Gaps in the hit sequence (original and newly aligned)
  gaps.orig       <- unlist(strsplit(as.character(aln.hit.orig),"[A-Z]"))
  gaps.count.orig <- nchar(gaps.orig)
  gaps.lead.orig  <- gaps.count.orig[1]
  gaps.trail.orig <- gaps.count.orig[length(gaps.count.orig)]
  
  gaps.new        <- unlist(strsplit(as.character(gsub("-","",aln.hit.new)),"[A-Z]"))
  gaps.count.new  <- nchar(gaps.new)
  gaps.lead.new   <- gaps.count.new[1]
  gaps.trail.new  <- gaps.count.new[length(gaps.count.new)]
  
  if(length(gaps.count.orig)>length(gaps.count.new)){
    gaps.count.new <- append(gaps.count.new,0)
  }
  
  gaps.discrep     <- suppressWarnings(rbind(gaps.count.orig,gaps.count.new))
  gaps.discrep.num <- gaps.discrep[1,]-gaps.discrep[2,]
  gaps.lead.add    <- gaps.discrep.num[1]
  gaps.trail.add   <- gaps.discrep.num[length(gaps.discrep.num)]
  
  # New alignment
  aln.add4  <- c(rep("-",gaps.lead.add),
                 aln.add2,
                 rep("-",gaps.trail.add))
  
  seq.alignable <- Biostrings::pairwiseAlignment(seq.d,
                                                 as.AAstring(aln.add3, degap=TRUE),
                                                 substitutionMatrix = BLOSUM40,
                                                 gapOpening   = 0,
                                                 gapExtension = 8)
  length(aln.add4)==ncol(MSA)
  
  query     <- aln.add4
  
  aln.final <- rbind(query,MSA)
  output    <- list(MSA             = aln.final,
                    aln.hit.name    = aln.hit.name,
                    aln.hit.seq     = aln.hit.seq,
                    aln.hit.score   = aln.hit.score,
                    aln.all.score   = aln.all/aln.limit,
                    seq.unalignable = seq.alignable)
  output
}



seq.rotate <- function(SAPCA,newseq){
  
  res.props  <- colnames(SAPCA$numerical.alignment$res.prop)
  prop.means <- SAPCA$numerical.alignment$prop.means
  prop.vars  <- SAPCA$numerical.alignment$prop.vars
  # Align new sequence with MSA
  
  # Numericise new sequence
  newseq.num <- numericise_MSA(MSA      = newseq$MSA[c("query",newseq$aln.hit.name),],
                               res.prop = SAPCA$numerical.alignment$res.prop)
  
  # Scale new sequnce using same scaling as SAPCA (gaps as "NA")
  newseq.scale.stack <- NULL
  for (x in 1:length(res.props)) {
    newseq.scale.stack[[res.props[x]]] <- (newseq.num$MSA.num.stack[[res.props[x]]]- prop.means[x]) /
      sqrt(prop.vars[x])
  }
  
  # Reflow into single wide matrix
  newseq.scale.wide <- newseq.scale.stack[[1]]
  for (x in 2:length(res.props)) {
    newseq.scale.wide <- cbind(newseq.scale.wide, newseq.scale.stack[[res.props[x]]])
  }
  
  # Replace gaps (currently "NA") with column average of the scaled SAPCA
  # Create na.colmean function
  gapvalues  <- colMeans(SAPCA$numerical.alignment$MSA.scale.wide)
  na.replace <- function(x,y){
    x[is.na(x)] <- y
    x
  }
  newseq.scale.wide.g <- NULL
  for(i in 1:ncol(newseq.scale.wide)){
    newseq.scale.wide.g <- cbind(newseq.scale.wide.g,
                                 na.replace(newseq.scale.wide[,i],gapvalues[i]))
  }
  
  # Rotate scaled sequence into same space as SAPCA
  newseq.rot <- scale(newseq.scale.wide.g,
                      SAPCA$seq.space.PCA$centre,
                      SAPCA$seq.space.PCA$scale) %*% SAPCA$seq.space.PCA$loadings
  
  # Output
  output <- list(seq       = newseq,
                 seq.num   = newseq.num$MSA.num.wide,
                 seq.scale = newseq.scale.wide.g,
                 seq.rot   = newseq.rot)
  output
}



seq.clust.add <- function(SAPCA,newseq.r){
  
  SAPCA.c  <- mclustrev(SAPCA)
  newseq.c <- mclust::predict.Mclust(SAPCA.c,newseq.r$seq.rot[,SAPCA$call$clusterPCs])
  newseq.c
}



mclustrev <- function(SAPCA){
  output <- SAPCA$seq.space.clusters$other
  
  output$classification <- SAPCA$seq.space.clusters$classification
  output$G              <- SAPCA$seq.space.clusters$optimal
  output$z              <- SAPCA$seq.space.clusters$likelihoods
  
  class(output) <- "Mclust"
  output
}



seq.SAPCA.add <- function (SAPCA,newseq.r,newseq.c){
  output <- SAPCA
  output$numerical.alignment$seq.names      <- rownames(newseq.r$seq$MSA)
  output$numerical.alignment$MSA            <- newseq.r$seq$MSA
  output$numerical.alignment$MSA.num.wide   <- rbind(newseq.r$seq.num[1,],
                                                     SAPCA$numerical.alignment$MSA.num.wide)
  output$numerical.alignment$MSA.scale.wide <- rbind(newseq.r$seq.scale[1,],
                                                     SAPCA$numerical.alignment$MSA.scale.wide)
  output$numerical.alignment$MSA.num.stack  <- NULL
  output$numerical.alignment$MSA.scale.stack<- NULL
  
  output$seq.space.PCA$coordinates          <- rbind(newseq.r$seq.rot[1,],
                                                     SAPCA$seq.space.PCA$coordinates)
  output$seq.space.clusters$likelihoods     <- rbind(newseq.c$z[1,],
                                                     SAPCA$seq.space.clusters$likelihoods)
  output$seq.space.clusters$classification  <- c(newseq.c$classification[1],
                                                 SAPCA$seq.space.clusters$classification)
  
  rownames(output$numerical.alignment$MSA.num.wide)[1]   <- "query"
  rownames(output$numerical.alignment$MSA.scale.wide)[1] <- "query"
  rownames(output$seq.space.PCA$coordinates)[1]          <- "query"
  rownames(output$seq.space.clusters$likelihoods)[1]     <- "query"
  
  output
}



seq.add.full <- function (SAPCA,sequence,SAPCAname=NULL){
  
  newseq   <- seq.MSA.add(SAPCA,sequence,SAPCAname)
  newseq.r <- seq.rotate(SAPCA,newseq)
  newseq.c <- seq.clust.add(SAPCA,newseq.r)
  SAPCA2   <- seq.SAPCA.add(SAPCA,newseq.r,newseq.c)
  
  SAPCA2
}



#BLOSUM40
blosum <- function(file="C:\\Users\\T\\OneDrive\\0-Sequences\\2-PCA\\0-Raw data and scalers\\0 - BLOSUM40.csv"){
  BLOSUM40 <- read.csv(file)
  BLOSUM40.names <- BLOSUM40[,1]
  BLOSUM40 <- BLOSUM40[,-1]
  BLOSUM40 <- as.matrix(BLOSUM40)
  rownames(BLOSUM40)<-BLOSUM40.names
  colnames(BLOSUM40)<-BLOSUM40.names
  BLOSUM40
}



percent <- function(x, digits = 1, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}



plot_3Dclusters <- function(SAPCA,
                            plotPCs    = 1:3,
                            col        = "cluster",
                            radius     = 1,
                            labels     = NULL,
                            write      = FALSE,
                            axeslabels = "default"){
  if (!is.null(SAPCA$seq.space.PCA$coordinates)){
    data <- SAPCA$seq.space.PCA$coordinates
  }else{
    data <- SAPCA    
  }
  
  if (all(col=="cluster")){
    colour <- SAPCA$seq.space.clusters$classification
  }else{
    colour <- col
  }
  # Calculate radius size
  rad <- (range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[2]-range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[1])/100
  rad <- rad*radius
  
  if (all(axeslabels=="default")){
    axes <- paste("PC",plotPCs,sep="")
  }else{
    axes <- axeslabels
  }
  if (is.null(axeslabels)){
    axes <- c("","","")
  }
  
  # Plot model-based clusters in 3D
  rgl::plot3d(data[,plotPCs],
              col      = colour,      # colour by clusters
              specular = "black",     # matte lighting
              type     = "s",         # "p" is points, "s" is spheres
              radius   = rad,         # sphere radius if using spheres
              size     = 4,           # point size
              axes     = FALSE,       # draw axes separately
              xlab     = axes[1],
              ylab     = axes[2],
              zlab     = axes[3])       
  # Draw axes
  if (write!=FALSE){
    rgl::axes3d(color = "black", labels = FALSE)                       
  }else{
    rgl::axes3d(color = "black", alpha=0.5, labels = FALSE) 
  }
  
  for (NAME in labels){
    SUB = row.names(SAPCA$seq.space.PCA$coordinates)==NAME      # Label based on its row.name
    rgl::text3d(subset(SAPCA$seq.space.PCA$coordinates[,plotPCs],subset=SUB), 
                text      = paste('---',NAME),   # data label text
                font      = 2,                   # bold
                color     = "black",             # colour
                adj       = -rad/2)              # offset
  }
  
  # Write html for interactive data
  if (write!=FALSE){
    rglwidget::.writeWebGL(write)                          
  }
}

