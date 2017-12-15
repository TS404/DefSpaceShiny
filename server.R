
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

# Import necessary files ----------------------------------------------------------------------
library("shiny")
library("rgl")
library("ggplot2")
library("Biostrings")
library("DECIPHER")
library("mclust")
library("tidyr")
library("rJava")

folder    <- "C:\\Users\\T\\OneDrive\\1-Scripts\\GitHub\\DefSpaceShiny"
setwd(folder)

SAPCA.cis <- readRDS("data\\CisDef.reference.PCA.RDS")
SAPCA.tra <- readRDS("data\\TransDef.reference.PCA.RDS")

view.cis <- readRDS("data\\CisDef.viewangle.RDS")
view.tra <- readRDS("data\\TransDef.viewangle.RDS")

BLOSUM40 <- readRDS("data\\BLOSUM.RDS")

motifs <- readRDS("data\\cysteine_motifs.RDS")

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

