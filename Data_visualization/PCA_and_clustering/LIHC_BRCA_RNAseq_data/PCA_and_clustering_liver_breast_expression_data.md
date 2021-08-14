PCA\_and\_clustering\_liver\_breast\_expression\_data
================
Christina Kosanovic
8/13/2021

This example provides the steps to prepare a Principal Component
Analysis and clustering (kmeans, hierarchical) of a TCGA breast cancer
and liver cancer expression data set

Load necessary libraries for PCA and clustering analysis:

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.0.5

``` r
library(ggfortify)
```

    ## Warning: package 'ggfortify' was built under R version 4.0.5

``` r
library(stats)
library(tools)
```

Load data set into an object and view its format:

``` r
data <- read.table("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/LIHC_BRCA_data1_marked_no0.txt", header = TRUE, row.names=1, stringsAsFactors = F, check.names=F)
View(data)
```

Save the original format of data in an object:

``` r
Origin_data <- data
```

With a quick check we can see names are assigned to first column and
correspond to sample name:

``` r
names(data) 
```

    ##  [1] "TCGA-2V-A95S-01A-11R-A37K-07_LIHC_TP"
    ##  [2] "TCGA-2Y-A9GS-01A-12R-A38B-07_LIHC_TP"
    ##  [3] "TCGA-2Y-A9GU-01A-11R-A38B-07_LIHC_TP"
    ##  [4] "TCGA-BC-A10Q-11A-11R-A131-07_LIHC_NT"
    ##  [5] "TCGA-BC-A10T-11A-11R-A131-07_LIHC_NT"
    ##  [6] "TCGA-BC-A10W-11A-11R-A131-07_LIHC_NT"
    ##  [7] "TCGA-3C-AAAU-01A-11R-A41B-07_BRCA_TP"
    ##  [8] "TCGA-3C-AALJ-01A-31R-A41B-07_BRCA_TP"
    ##  [9] "TCGA-3C-AALK-01A-11R-A41B-07_BRCA_TP"
    ## [10] "TCGA-A7-A0CE-11A-21R-A089-07_BRCA_NT"
    ## [11] "TCGA-A7-A0CH-11A-32R-A089-07_BRCA_NT"
    ## [12] "TCGA-A7-A0D9-11A-53R-A089-07_BRCA_NT"

We are interested in variance among the samples of this data set.  
Assigning the Sample names and Groups (sample + class) into objects:

``` r
Group <- t(data[1,]) # Groups and Samples from the transposed first index of rows stored in Groups object
samples <- colnames(data)
print(Group)
```

    ##                                      class    
    ## TCGA-2V-A95S-01A-11R-A37K-07_LIHC_TP "LIHC_TP"
    ## TCGA-2Y-A9GS-01A-12R-A38B-07_LIHC_TP "LIHC_TP"
    ## TCGA-2Y-A9GU-01A-11R-A38B-07_LIHC_TP "LIHC_TP"
    ## TCGA-BC-A10Q-11A-11R-A131-07_LIHC_NT "LIHC_NT"
    ## TCGA-BC-A10T-11A-11R-A131-07_LIHC_NT "LIHC_NT"
    ## TCGA-BC-A10W-11A-11R-A131-07_LIHC_NT "LIHC_NT"
    ## TCGA-3C-AAAU-01A-11R-A41B-07_BRCA_TP "BRCA_TP"
    ## TCGA-3C-AALJ-01A-31R-A41B-07_BRCA_TP "BRCA_TP"
    ## TCGA-3C-AALK-01A-11R-A41B-07_BRCA_TP "BRCA_TP"
    ## TCGA-A7-A0CE-11A-21R-A089-07_BRCA_NT "BRCA_NT"
    ## TCGA-A7-A0CH-11A-32R-A089-07_BRCA_NT "BRCA_NT"
    ## TCGA-A7-A0D9-11A-53R-A089-07_BRCA_NT "BRCA_NT"

``` r
print(samples)
```

    ##  [1] "TCGA-2V-A95S-01A-11R-A37K-07_LIHC_TP"
    ##  [2] "TCGA-2Y-A9GS-01A-12R-A38B-07_LIHC_TP"
    ##  [3] "TCGA-2Y-A9GU-01A-11R-A38B-07_LIHC_TP"
    ##  [4] "TCGA-BC-A10Q-11A-11R-A131-07_LIHC_NT"
    ##  [5] "TCGA-BC-A10T-11A-11R-A131-07_LIHC_NT"
    ##  [6] "TCGA-BC-A10W-11A-11R-A131-07_LIHC_NT"
    ##  [7] "TCGA-3C-AAAU-01A-11R-A41B-07_BRCA_TP"
    ##  [8] "TCGA-3C-AALJ-01A-31R-A41B-07_BRCA_TP"
    ##  [9] "TCGA-3C-AALK-01A-11R-A41B-07_BRCA_TP"
    ## [10] "TCGA-A7-A0CE-11A-21R-A089-07_BRCA_NT"
    ## [11] "TCGA-A7-A0CH-11A-32R-A089-07_BRCA_NT"
    ## [12] "TCGA-A7-A0D9-11A-53R-A089-07_BRCA_NT"

Reloading data, skipping first row to remove sample IDs from the
expression table:

``` r
data <- read.table("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/LIHC_BRCA_data1_marked_no0.txt", header = TRUE, row.names=1, skip =1, stringsAsFactors = F, check.names=F)
View(data)
samples <- colnames(data)
print(samples)
```

    ##  [1] "LIHC_TP" "LIHC_TP" "LIHC_TP" "LIHC_NT" "LIHC_NT" "LIHC_NT" "BRCA_TP"
    ##  [8] "BRCA_TP" "BRCA_TP" "BRCA_NT" "BRCA_NT" "BRCA_NT"

Transform table to numeric to convert character objects –&gt; numeric
form:

``` r
data <- transform(data, as.numeric())
str(data)  ## check that data is num form 
```

    ## 'data.frame':    16837 obs. of  12 variables:
    ##  $ LIHC_TP: num  0 2.31 5.69 138.3 1561 ...
    ##  $ LIHC_TP: num  0 53.59 5.41 144.07 1297 ...
    ##  $ LIHC_TP: num  0 6.86 6.14 73.93 1423 ...
    ##  $ LIHC_NT: num  0 2 0 104 1454 ...
    ##  $ LIHC_NT: num  0 1.41 2.59 96.89 1125 ...
    ##  $ LIHC_NT: num  0 4.94 1.06 97.03 2128 ...
    ##  $ BRCA_TP: num  0 16.4 12.9 52.2 408.1 ...
    ##  $ BRCA_TP: num  0.907 11.623 9.229 154.297 1360.834 ...
    ##  $ BRCA_TP: num  0 12.1 11.1 143.9 865.5 ...
    ##  $ BRCA_NT: num  0 4.33 3.92 78.92 978.41 ...
    ##  $ BRCA_NT: num  0 4.21 2.19 53.64 970.76 ...
    ##  $ BRCA_NT: num  0 3.06 0 87.58 770.37 ...

Transposing the expressions matrix to move genes from rows –&gt; columns
and samples from columns –&gt; rows. This will ensure samples are
treated as the objects and genes as the features to analyze sample
variation instead of gene variation.

``` r
data_T <- t(data)
View(data_T)
```

## Principal Component Analysis

Performing the Principal Component Analysis and placing it in an
object after scaling and centering the data:

``` r
pca <- prcomp(data_T, scale. = TRUE, center = TRUE)
```

View the results:

``` r
summary(pca)
```

    ## Importance of components:
    ##                            PC1     PC2     PC3     PC4      PC5      PC6
    ## Standard deviation     64.9555 57.6242 46.2684 43.5985 35.54857 33.43112
    ## Proportion of Variance  0.2506  0.1972  0.1272  0.1129  0.07505  0.06638
    ## Cumulative Proportion   0.2506  0.4478  0.5750  0.6878  0.76291  0.82929
    ##                            PC7     PC8      PC9     PC10     PC11     PC12
    ## Standard deviation     31.3300 26.1770 21.29679 20.61981 18.13196 2.47e-13
    ## Proportion of Variance  0.0583  0.0407  0.02694  0.02525  0.01953 0.00e+00
    ## Cumulative Proportion   0.8876  0.9283  0.95522  0.98047  1.00000 1.00e+00

Creates object pca$x which holds PCs:

``` r
pca_result <- data.frame(pca$x, Group)
```

Save PCA table:

``` r
write.table(pca_result, file="PCA_table1.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, append=TRUE)
```

    ## Warning in write.table(pca_result, file = "PCA_table1.txt", row.names = FALSE, :
    ## appending column names to file

Create object based visualization with ggplot using the PCA result:

``` r
ggplot(pca_result, aes(x=PC1, y=PC2, color=Group)) +geom_point() +stat_ellipse(level=0.4) 
```

    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse

    ## Warning: Removed 4 row(s) containing missing values (geom_path).

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Create object based visualization with autoplot using the PCA result:

``` r
plot_pca <- autoplot(pca, label = TRUE, label.size = 3, colour = "red")
plot(plot_pca)
```

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Plot variance for each principal component (PC1 will explain most
variance in our data):

``` r
plot(pca)
```

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## K-means Clustering

Each object will be assigned to a group depending upon similarity or
distance between them.

K-means is the most popular method to partition observations into
clusters. K-means will take first random point in data and from that
point, create a centroid and find the closest element to that point.
Then, those 2 points become one and it searches for the next closest
point to the new centroid and so on.

Using the transposed data set created previously:

``` r
kmeans <- kmeans(data_T,4)
cluster_result <- data.frame(kmeans$cluster, Group)
View(cluster_result)
```

Plotting clusters:

``` r
autoplot(kmeans, data =data_T, label = TRUE, label.size = 3, frame = TRUE, frame.type ='t')
```

    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

K-means doesn’t give us the position but assigns numbers to each one of
our samples.  
We may achieve a better result with an alternate form of clustering.

# Hierarchical Clustering

A good alternative to k-means is hierarchical clustering.  
In this example, ward.d2 is the algorithm for hierarchical clustering.

Reload and define the data:

``` r
data <- read.table ("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/LIHC_BRCA_data1_marked_no0.txt", header = TRUE, row.names = 1)
```

Transpose data to ensure sample rather than gene comparison:

``` r
data1 <- t(data)
DataT <- cbind(data1, Group)
View(DataT)
```

Define clustering and save in an object.  
Specify distance matrix (euclidean / manhattan / maximum / canberra /
binary / minkowski) and linkage(single / complete / average / mean /
centroid / ward.D / ward.D2) type.

``` r
hclustering <- hclust(dist(DataT, method = 'euclidean'), method='ward.D2')
```

    ## Warning in dist(DataT, method = "euclidean"): NAs introduced by coercion

``` r
plot(hclustering)
```

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

Define parameters and save in an object:

``` r
descr1 <- paste ("Distance: ", "euclidean", sep="")
descr2 <- paste ("Linkage: ", "ward.D2", sep="")
```

Plot clusters to create dendrogram:

``` r
plot(hclustering, xlab=descr1, sub=descr2, cex = 0.8)
```

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Print cluster order:

``` r
hclustering$order
```

    ##  [1]  2  9 12 10 11  7  8  5  4  6  1  3

Observe some inconsistency in the clustering results. While the tissue
type (breast and liver) have been separated better than the normal
tissue vs. tumor, one of the liver (tumor) samples appears to be
misplaced.

Often, results can appear skewed with raw data that has not been
filtered or normalized to account for noise/experimental artifacts.

Typically, RNA-Seq data does not have a normal distribution and requires
quantile normalization.  
Let’s observe how the output of PCA and clustering varies with quantile
normalization and removal of lowly expressed genes.

## PCA and clustering with normalization

``` r
## Repeating steps above for loading data into object 
data <- read.table("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/LIHC_BRCA_data1_marked_no0.txt", header = TRUE, row.names=1, stringsAsFactors = F, check.names=F)
Group <- t(data[1,]) # Groups and Samples from the transposed first index of rows stored in Groups object
samples <- colnames(data)
## Reloading data
data <- read.table("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/LIHC_BRCA_data1_marked_no0.txt", header = TRUE, row.names=1, skip =1, stringsAsFactors = F, check.names=F)
samples <- colnames(data)
dim(data)
```

    ## [1] 16837    12

Filter out lowly expressed genes and normalize data:

``` r
# Requiring the mean count to be greater than 1000 :  
datanew = data[rowMeans(data) >1000,] 
# check dimensions of filtered data 
dim(datanew)
```

    ## [1] 4428   12

``` r
# Performing log normalization :
log_data = log(datanew + 1) 
# Transforming data to numeric and transposing expression table :
data_t <- transform(log_data, as.numeric())
data_T <- t(data_t)
```

Performing PCA analysis on normalized data:

``` r
pca <- prcomp(data_T, scale. = TRUE, center = TRUE)
## View the results
summary(pca)
```

    ## Importance of components:
    ##                            PC1     PC2     PC3      PC4      PC5      PC6
    ## Standard deviation     38.3941 33.1799 24.0208 20.46770 16.29017 13.48524
    ## Proportion of Variance  0.3329  0.2486  0.1303  0.09461  0.05993  0.04107
    ## Cumulative Proportion   0.3329  0.5815  0.7118  0.80644  0.86637  0.90744
    ##                             PC7    PC8     PC9    PC10    PC11      PC12
    ## Standard deviation     12.70701 9.4106 8.63974 7.02813 5.98112 5.856e-14
    ## Proportion of Variance  0.03647 0.0200 0.01686 0.01116 0.00808 0.000e+00
    ## Cumulative Proportion   0.94391 0.9639 0.98077 0.99192 1.00000 1.000e+00

``` r
## Creates object pca$x which hold PCs
pca_result <- data.frame(pca$x, Group)

#Save PCA table
write.table(pca_result, file="PCA_table2.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, append=TRUE)
```

    ## Warning in write.table(pca_result, file = "PCA_table2.txt", row.names = FALSE, :
    ## appending column names to file

``` r
## Create object based visualization using the PCA result
ggplot(pca_result, aes(x=PC1, y=PC2, color=Group)) +geom_point() +stat_ellipse(level=0.4) 
```

    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse

    ## Warning: Removed 4 row(s) containing missing values (geom_path).

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
plot_pca <- autoplot(pca, label = TRUE, label.size = 3, colour = "red")
plot(plot_pca)
```

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

Performing k-means clustering on normalized data:

``` r
## Using the transposed dataset created previously 
kmeans <- kmeans(data_T,4)
cluster_result <- data.frame(kmeans$cluster, Group)
View(cluster_result)

#Plot clusters
autoplot(kmeans, data =data_T, label = TRUE, label.size = 3, frame = TRUE, frame.type ='t')
```

    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse
    ## Too few points to calculate an ellipse

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

Performing hierarchical clustering on normalized data:

``` r
## Loading and defining data
data <- read.table("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/LIHC_BRCA_data1_marked_no0.txt", header = TRUE, row.names=1, stringsAsFactors = F, check.names=F)
Group <- t(data[1,]) # Groups and Samples from the transposed first index of rows stored in Groups object
samples <- colnames(data)
## Reloading data  
data <- read.table("https://raw.githubusercontent.com/PineBiotech/omicslogic/master/LIHC_BRCA_data1_marked_no0.txt", header = TRUE, row.names=1, skip =1, stringsAsFactors = F, check.names=F)
samples <- colnames(data)
dim(data)
```

    ## [1] 16837    12

``` r
#Filter out lowly expressed genes:  
datanew = data[rowMeans(data) >1000,] 
dim(datanew)
```

    ## [1] 4428   12

``` r
## Performing log normalization 
log_data = log(datanew + 1) 
## Transforming data to numeric and transposing expression table 
data_t <- transform(log_data, as.numeric())
data_T <- t(data_t)
DataT <- cbind(data_T, Group)
View(DataT)

# Define clustering and save in a object
hclustering <- hclust(dist(DataT, method = 'euclidean'), method='ward.D2')
```

    ## Warning in dist(DataT, method = "euclidean"): NAs introduced by coercion

``` r
plot(hclustering)
```

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
#Define parameters and save in object
descr1 <- paste ("Distance: ", "euclidean", sep="")
descr2 <- paste ("Linkage: ", "ward.D2", sep="")
#Plot clusters to create dendrogram 
plot(hclustering, xlab=descr1, sub=descr2, cex = 0.8)
```

![](PCA_analysis_and_clustering_liver_breast_expression_data_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

``` r
#Print cluster order
hclustering$order
```

    ##  [1]  7  8  9 12 10 11  5  4  6  3  1  2

In this case, filtering and normalization appears to have fixed the
observed issue in previous hierarchical clustering. Now the samples are
correctly separated by tissue type and normal-like vs. tumor, where the
greatest difference is observed between the two tissue types (breast and
liver).
