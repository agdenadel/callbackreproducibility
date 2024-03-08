library(mvtnorm)
library(ggplot2)

# todo
# clustering metrics vignette

if (FALSE) {

n = 30

group1 <- rmvnorm(n=n, mean=c(0,0), sigma=diag(2))
group2 <- rmvnorm(n=n, mean=c(6,6), sigma=diag(2))
group3 <- rmvnorm(n=n, mean=c(0,6), sigma=diag(2))


data <- as.data.frame(rbind(group1, group2, group3))


underclustered <- kmeans(data, centers = 2)$cluster
overclustered <- kmeans(data, centers = 4)$cluster

group <- c(rep(1, n), rep(2, n), rep(3, n))

data$group <- factor(group)
data$overclustered <- factor(overclustered)
data$underclustered <- factor(underclustered)



ggplot(data, aes(x=V1, y=V2, col=group)) + 
  geom_point(size=5) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  NoLegend() +
  scale_color_brewer(palette = "Set1")
 
ggplot(data, aes(x=V1, y=V2, col=underclustered)) + 
  geom_point(size=5) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  NoLegend() +
  scale_color_brewer(palette = "Set1")

ggplot(data, aes(x=V1, y=V2, col=overclustered)) + 
  geom_point(size=5) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  NoLegend() +
  scale_color_brewer(palette = "Set1")


# overclustered
pdfCluster::adj.rand.index(data$group, data$overclustered)
clusteval::jaccard(data$group, data$overclustered)
clevr::fowlkes_mallows(data$group, data$overclustered)

clevr::v_measure(data$group, data$overclustered)
clevr::completeness(data$group, data$overclustered)
clevr::homogeneity(data$group, data$overclustered)
aricode::NMI(data$group, data$overclustered)

# underclustered
pdfCluster::adj.rand.index(data$group, data$underclustered)
clusteval::jaccard(data$group, data$underclustered)
clevr::fowlkes_mallows(data$group, data$underclustered)

clevr::v_measure(data$group, data$underclustered)
clevr::completeness(data$group, data$underclustered)
clevr::homogeneity(data$group, data$underclustered)
aricode::NMI(data$group, data$underclustered)


}