library(MASS)
library(labdsv)
library(cluster)
library(pvclust)
library(ade4)
library(factoextra)
library(FField)


### This script runs PCA with two different packages. 
### Then, it performs hierarchical clustering.
### Then it calculates GLM


#read file
data<-read.csv(file="abiotic_factors.csv", head=F, row.names=1)

labels = c(rownames(data))

#PCA
data.pca<-pca(data, cor=TRUE)

data.pca.scores<-scores.pca(data.pca, dim=7)

#kmeans of PCA scores
data.pca.scores.cl<-kmeans(data.pca.scores, 3)

res.pca <- dudi.pca(data,
                    scannf = FALSE,   # Hide scree plot
                    nf = 5            # Number of components kept in the results
)

#access PCA results
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs


groups <- as.factor(data.pca.scores.cl$cluster)

fviz_pca_ind(res.pca,
             habillage = c(groups), # Color by the quality of representation
             repel = TRUE,     # Avoid text overlapping
             label="all",
             legend.title="Zone"
)

#plots subjects and components
scatter(res.pca,
        posieig = "none", # Hide the scree plot
        clab.row = 1      # Hide row labels
        
)

s.class(res.pca$li,
        fac = groups,  # color by groups
        cpoint = 5,
        col = c("blue", "darkorange", "salmon", "peru"),
        addaxes = TRUE,
        label=c("Zone 1", "Zone 2", "Zone 3")
)






### HIERARCHICAL CLUSTERING ###

tdata=t(data)

#add tags column for color purposes
#zone_1<-rep("1", 3)
#zone_2<-rep("2", 3)
#zone_3<-rep("3", 2)
#tags<-append(append(zone_1, zone_2),zone_3)
#data.tag<-cbind(data,tags)

#perform hierarchical clustering on raw data
data.hc<-hclust(dist(data, method="euclidean"), method ="ave")
#plot (data.hc, hang = -1, labels=rownames(data))


#data.hc<-hclust(dist(data, method="manhattan"), method ="ave")
data.hc.t<-as.dendrogram(data.hc)

#calculate pvalues of the clusters
data.hc.t<-pvclust(tdata, method.hclust="ave",
               method.dist="euclidean")

plot(data.hc.t, type="triangle", main="All abiotic factors")


### GLM ###

FTO_data<-read.csv(file="all_data.csv", head=1, row.names=1)

FTOglm <- glm(meth~temp + sal + pH + vis, family=gaussian, FTO_data)

summary(FTOglm)

