library(Matrix)
library(huge)
library(ggplot2)
library(igraph)
library(matrixcalc)
library(latex2exp)
library(heterocop)

library(seriation)

data <- read.csv2("data.csv",sep=";")[,-1] # data

Sigma <- as.matrix(read.csv2("MLE_cor_RNA_disc.csv",sep=";")[,-1]) #PPMLE
Sig_b <- matrix(as.numeric(as.matrix(read.csv2("R.csv",sep=",")[,-1])),278,278) #bridge

#hist(Sigma)
#hist(Sig_b)

# Inverting the results for a grid of lambda
#copula
lambda <- seq(0.05,0.75,0.01)
res_PD <- matrix(nearPD(Sigma,corr=T,maxit=1000)$mat@x,278,278)
res_Pd_b <- matrix(nearPD(Sig_b,corr=T,maxit=1000)$mat@x,278,278)

HBIC <- c()

#selecting the lambda
for(l in 1:length(lambda)){
  print(l)
  OM_hat <- huge(x=res_PD,lambda=lambda[l],method="glasso")$icov[[1]]
  crit <- matrix.trace(res_PD %*%OM_hat)-log(det(OM_hat))+log(log(250))*(log(278)/250)*sum(OM_hat!=0)
  HBIC <- c(HBIC,crit)
}
# 0.53 for the copula

HBIC_b <- c()

for(l in 1:length(lambda)){
  print(l)
  OM_hat_b <- huge(x=res_Pd_b,lambda=lambda[l],method="glasso")$icov[[1]]
  crit_b <- matrix.trace(res_Pd_b %*%OM_hat_b)-log(det(OM_hat_b))+log(log(250))*(log(278)/250)*sum(OM_hat_b!=0)
  HBIC_b <- c(HBIC_b,crit_b)
}
# for the bridge function lambda is 0.05

# plotting the HBIC curves
data_HBIC <- data.frame(cbind(lambda,HBIC[5:75]))
colnames(data_HBIC) <- c("lambda","HBIC")
HBIC_min <- min(HBIC[-1])
lambda_min <- lambda[which(HBIC==HBIC_min)] #0.53

g <- ggplot(data_HBIC,aes(x=lambda,y=HBIC)) + geom_line()+theme_bw()+xlab(TeX("$\\lambda$"))+ylab(TeX("HBIC($\\lambda$)"))+geom_vline(xintercept=0.53,color="red",linetype="dashed")+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))

ggsave("HBIC.eps",
       plot = g,
       device = "eps",
       width = 9, height = 5)

data_HBICb <- data.frame(cbind(lambda,HBIC_b_new))
colnames(data_HBICb) <- c("lambda","HBIC")
HBIC_min_b <- min(HBIC_b_new)
lambda_min_b <- lambda[which(HBIC_b_new==HBIC_min_b)] #0.05

g <- ggplot(data_HBICb,aes(x=lambda,y=HBIC)) + geom_line()+theme_bw()+xlab(TeX("$\\lambda$"))+ylab(TeX("HBIC($\\lambda$)"))+geom_hline(yintercept=210,color="red",linetype="dashed")+geom_vline(xintercept=0.01,color="red",linetype="dashed")

ggsave("HBIC_bridge.eps",
       plot = g,
       device = "eps",
       width = 9, height = 5)


#same for gLasso has given 0.51

HBIC_gl <- read.csv("HBIC_glasso.csv",sep=",")[,2]

data_HBICb <- data.frame(cbind(lambda,HBIC_gl))
colnames(data_HBICb) <- c("lambda","HBIC")
g <- ggplot(data_HBICb,aes(x=lambda,y=HBIC)) + geom_line()+theme_bw()+xlab(TeX("$\\lambda$"))+ylab(TeX("HBIC($\\lambda$)"))+geom_hline(yintercept=312,color="red",linetype="dashed")+geom_vline(xintercept=0.51,color="red",linetype="dashed")

ggsave("HBIC_glasso.eps",
       plot = g,
       device = "eps",
       width = 9, height = 5)



## Computing the matrices
OM <- huge(x=res_PD,lambda=0.53,method="glasso")$icov[[1]]
colnames(OM) <- colnames(data)
rownames(OM) <- colnames(data)


OMb <- huge(x=res_Pd_b,lambda=0.05,method="glasso")$icov[[1]]
colnames(OMb) <- colnames(data)
rownames(OMb) <- colnames(data)


OM_gl<- huge(x=as.matrix(data),lambda=0.51,method="glasso")$icov[[1]]
colnames(OM_gl) <- colnames(data)
rownames(OM_gl) <- colnames(data)



# values for Table 4

liens <- function(M){
RR <- ((sum(M[1:108,1:108]!=0)-108)/2)
PP <- ((sum(M[109:216,109:216]!=0)-108)/2)
MM <- ((sum(M[217:278,217:278]!=0)-62)/2)
RP <- sum(M[1:108,109:216]!=0)
RM <- sum(M[1:108,217:278]!=0)
PM <- sum(M[109:216,217:278]!=0)
return(c(RR,PP,MM,RP,RM,PM))

}

liens(OM)

# Figure 5 of the Supplmentary : number of links detected depending on lambda.
RR <- c()
PP <- c()
MM <- c()
RP <- c()
RM <- c()
PM <- c()

for(l in 1:length(lambda)){
  print(l)
  OM_hat <- huge(x=res_PD,lambda=lambda[l],method="glasso")$icov[[1]]
  
  RR <- c(RR, (sum(OM_hat[1:108,1:108]!=0)-108)/2)
  PP <- c(PP, (sum(OM_hat[109:216,109:216]!=0)-108)/2)
  MM <- c(MM, (sum(OM_hat[217:278,217:278]!=0)-62)/2)
  RP <- c(RP, sum(OM_hat[1:108,109:216]!=0))
  RM <- c(RM, sum(OM_hat[1:108,217:278]!=0))
  PM <- c(PM, sum(OM_hat[109:216,217:278]!=0))

}

data_links <- data.frame(rep(lambda,6),c(RR/5778,PP/5778,MM/1891,RP/11664,RM/6696,PM/6696),rep(c("RNA-RNA", "protein-protein","mutation-mutation","RNA-protein","RNA-mutation","protein-mutation"),each=71))
colnames(data_links) <- c("Lambda","Proportion","Type")

g <- ggplot(data_links,aes(x=Lambda,y=Proportion,color=Type))+xlab(TeX("$\\lambda$"))+ scale_color_brewer(palette="Dark2") + geom_line(size=1.2)+ylim(0,0.5)+theme_bw()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15),legend.text = element_text(size=15),legend.title = element_text(size=15))

ggsave("nb_links.eps",
       plot = g,
       device = "eps",
       width = 10, height = 5)

# Graph for Figure 6

OM2 <- OM
OM2[OM2!=0] <- 1

#creating a diamond shape for RNA-seq
my_diamond <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex","color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1.2 * params("vertex","size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  for (i in seq_len(nrow(coords))) {
    x <- coords[i, 1]
    y <- coords[i, 2]
    r <- vertex.size[i]
    
    xs <- x + c(0, -r, 0, r)
    ys <- y + c(r, 0, -r, 0)
    
    polygon(xs, ys, col=vertex.color[i], border="black")
  }
}

add_shape(shape="diamond", clip=shape_noclip, plot=my_diamond)

shapes <- c(rep("diamond",108),rep("circle",108),rep("square",62))


par(mar=c(0,0,0,0)+.1)
g <- graph_from_adjacency_matrix(round(OM2,1),mode="undirected",diag=F) %>% set_vertex_attr("color",value=c(rep(1,108),rep(2,108),rep(3,62)))%>% set_vertex_attr("shape",value=shapes)
plot(g,vertex.size=6,vertex.shape=V(g)$shape, vertex.label.cex=0.8, vertex.color=V(g)$color,vertex.label.color="black")


# Figure 7 of the main and of the Supplementary
dg<- sort(degree(g),decreasing=T)

d <- data.frame(dg)
g <- ggplot(d,aes(x=dg))+geom_histogram(fill="white",color="black")+theme_bw()+xlab("Degree")+theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))


ggsave("degree.eps",
       plot = g,
       device = "eps",
       width = 10, height = 5)

# by omics type
d_p <- degree(g)[1:108]
d_p <- data.frame(d_p)
g_d <- ggplot(d_p,aes(x=d_p))+geom_histogram(fill="white",color="black")+theme_bw()+xlab("Degree")+theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))

ggsave("degree_protein.eps",
       plot = g_d,
       device = "eps",
       width = 10, height = 5)

d_r <- degree(g)[109:216]
d_r <- data.frame(d_r)
g_d <- ggplot(d_r,aes(x=d_r))+geom_histogram(fill="white",color="black")+theme_bw()+xlab("Degree")+theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))

ggsave("degree_arn.eps",
       plot = g_d,
       device = "eps",
       width = 10, height = 5)


d_m <- degree(g)[217:278]
d_m <- data.frame(d_m)
g_d <- ggplot(d_m,aes(x=d_m))+geom_histogram(fill="white",color="black")+theme_bw()+xlab("Degree")+theme(axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))

ggsave("degree_mutations.eps",
       plot = g_d,
       device = "eps",
       width = 10, height = 5)


# Figure 12 of the Supplementary
OM2b <- OMb
OM2b[OM2b!=0] <- 1
rownames(OM2b) <- colnames(data)
gb <- graph_from_adjacency_matrix(OM2b,mode="undirected",diag=F) %>% set_vertex_attr("color",value=c(rep(1,108),rep(2,108),rep(3,62)))


# Figure 8

high_degree_nodes <- V(g)[degree(g) > 10]
neighbors <- unlist(neighborhood(g, order = 1, nodes = high_degree_nodes))
neighbors <- unique(neighbors)

subgraph <- induced_subgraph(g, neighbors)
V(subgraph)$color[(V(subgraph)$name%in%names(high_degree_nodes))] <- "red"

plot(subgraph, layout = layout_with_fr(subgraph),vertex.size=8,vertex.label.color="black")



# Figure 9
G <- ego(g, order=1, nodes = c("PTEN","ACACA","GAB2","RAF1","SMAD1","RPS6KB1","IGFBP2","EIF4EBP1","YAP1"), mode = "all", mindist = 0)
subG <- induced_subgraph(g, c("PTEN","ACACA","GAB2","RAF1","SMAD1","RPS6KB1","IGFBP2","EIF4EBP1","YAP1","PTEN.1","ACACA.1","GAB2.1","RAF1.1","SMAD1.1","RPS6KB1.1","IGFBP2.1","EIF4EBP1.1","YAP1.1"))
plot(subG,layout=layout_with_fr(subG),vertex.label.color="black")

save.image("mutations_genes.RData")

## Obtaining Table 5

# we assimilate the mutations to the affected genes
edb <- EnsDb.Hsapiens.v86

seqlevelsStyle(edb) <- "UCSC"

#selecting the genes for which the mutations have been kept
gene_mut <- read.csv("gene_mut.csv")[,-1]
gene_mut <- gene_mut %>% dplyr::filter(mut.icgc_mutation_id%in%colnames(data))
#%>% dplyr::filter(mut.consequence_type=="exon_variant")
gene_mut <- unique(gene_mut[,c(1,2)])

# il y a seulement 30 mutations qui correspondent à des exon_variant

#test avec autre chose que les exons
genes_mut <- c()
for(i in 1:dim(gene_mut)[1]){
  if(length(ensembldb::select(edb, keys = ~ gene_id == gene_mut[i,2], columns = "SYMBOL")$SYMBOL)>0){
    genes_mut <- c(genes_mut,ensembldb::select(edb, keys = ~ gene_id == gene_mut[i,2], columns = "SYMBOL")$SYMBOL)
  }else{
    genes_mut <- c(genes_mut,"None")
  }
}

mut_gene <- cbind(gene_mut[,1],genes_mut)


mut_gene[which(mut_gene[,1]%in%names(V(subG))),]



library(ggraph)
library(reshape2)
library(ggplot2)
library(graphlayouts)

#Figure 8
comp <- components(g)
giant <- induced_subgraph(g, which(components(g)$membership ==
                                     which.max(components(g)$csize)))
E(g)$weight <- c(rep(1,3),6,rep(1,313))
E(g)$color <- c(rep("black",3),"red",rep("black",313))

plot(giant,layout=layout_with_stress(giant),vertex.size=6, vertex.label.cex=0.8, vertex.label.color="black",edge.width=E(giant)$weight,edge.color=E(giant)$color)

ggraph(giant, layout = "stress") +
  geom_edge_link(
    aes(width = weight),
    alpha = 0.6
  ) +
  scale_edge_width(range = c(0.4, 2.5)) +
  
  geom_node_point(
    aes(shape = I(V(giant)$shape)),
    colour = I(V(giant)$color),
    size = 4,
    stroke = 1
  ) +
  
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "none"
  )

# Figure 9
A <- as.matrix(as_adjacency_matrix(giant, sparse = FALSE))


ord <- seriate(A, method="PCA")

order_index <- get_order(ord)

M_reordered <- A[order_index, order_index, drop=FALSE]

df <- melt(M_reordered)

g_d <- ggplot(df, aes(Var1, Var2)) +
  geom_tile(aes(fill = value), color = "grey80") +
  scale_fill_gradient(low = "white", high = "black") +
  coord_equal() +
  labs(x = "", y = "", fill = "Poids") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
        axis.text.y = element_text(size=6),
        legend.position="none")

df_red <- subset(df, (Var1 == "CCNE1" & Var2 == "ANLN") | (Var2 == "CCNE1" & Var1 == "ANLN"))  

g_d <- g_d + geom_tile(data = df_red, aes(x = Var1, y = Var2),
                       fill = "red", color = "grey80")

ggsave("adj_matrix.eps",
       plot = g_d,
       device = "eps",
       width = 10, height = 10)
