## Assignment: Genomics course - Comparing DNA sequences using Phylogenetic Analysis  

## This script details a phylogenetic analysis to determine the evolutionary relationship between various isolates of H3N2.

## Load libraries

library(stats)
library(ape)
library(ade4)

## Load DNA sequences in FASTA format

dna <- read.dna(file = "https://code.omicslogic.com/resources/FD5Z", format = "fasta")
str(dna)

## Generate the distance matrix from dist, the object class that contains the distances between every pair of sequences

Dist <- dist.dna(dna, model = "TN93")

## Plotting the matrix of pairwise genetic distances in the form of a heat map

temp <- as.data.frame(as.matrix(Dist))
table.paint(temp, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)

## In the heatmap representing the distance matrix:  
## Darker shades of grey represent greater distances and lighter shades represent smaller distance and greater similarity.

## Improving the visualization by adding colors to the Distance matrix

## Convert into matrix

temp1 <- as.matrix(Dist)

## Define position of plot

par(mar = c(0.05, 4, 3.2, 0.05))

## Plot Distance matrix in the form of a colored heatmap

image(x = 1:40, y = 1:40, temp1, col = rev(heat.colors(100)), xaxt = "n", yaxt = "n", xlab = "Samples", ylab = "Samples")
axis(side = 2, at = 1:40, lab = rownames(dna), las = 2, cex.axis = 0.6)
axis(side = 3, at = 1:40, lab = rownames(dna), las = 3, cex.axis = 0.6)

## Darker shades represent greater distances and lighter shades represent smaller distances and greater similarity

## Building phylogenetic tree based on the classical Neighbor-Joining (NJ) method

tree_nj <- nj(Dist)

## Plotting tree

plot (tree_nj, cex = 0.6)

## In the phylogenetic tree, y-axis represents the time-axis and x-axis represents the order of branching with taxon.

## Building phylogenetic tree (dendrogram) using hclust method

tree_hclust <- hclust(Dist)

## Plotting tree

plot(tree_hclust, labels = NULL, hang = 0.1, check = TRUE, cex = 0.6, axes = TRUE, frame.plot = FALSE, ann = TRUE, main = "", sub = NULL, xlab = NULL, ylab = "Height")

