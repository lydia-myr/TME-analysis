#jacard index hclust
literature_sig<-read.xlsx("TME_related_pathway_done.xlsx",sheet=1)
literature_geneSet<-list()
for(i in seq(1,nrow(literature_sig))){
  tmp<-as.vector(unlist(literature_sig[i,]))
  title<-tmp[1]
  genes<-unique(tmp[-1])
  if(length(grep(TRUE,is.na(genes)))>0){
    genes<-genes[-(which(is.na(genes)==TRUE))]
  }
  literature_geneSet[[title]]<-genes
}

load("leading_edge_filtered_pathway.rda")#final_path_list
#ls()
load("literature_keyword_pathway.rda")
pathway_nes_ensembl<-c()
path_ensembl_name<-c()
for(key in names(keyword_path_list)){
    path<-intersect(keyword_path_list[[key]],row.names(pahtway_nes))
    if(length(path)==1){
      path_ensembl_name<-c(path_ensembl_name,key)
      pathway_nes_ensembl<-rbind(pathway_nes_ensembl,pahtway_nes[path,])
    }
    if(length(path)>1){
      path_ensembl_name<-c(path_ensembl_name,key)
      tmp<-apply(pahtway_nes[path,],2,function(x) (median(x)))
      pathway_nes_ensembl<-rbind(pathway_nes_ensembl,tmp)
    }
}
row.names(pathway_nes_ensembl)<-path_ensembl_name

pathway_item_liter<-read.table("literature_item_list.txt",header=F,sep="\t")

literature_geneSet_ense<-list()
for(i in seq(1,nrow(pathway_item_liter))){
  tmp<-as.vector(unlist(pathway_item_liter[i,]))
  title<-tmp[1]
  genes<-unique(tmp[-1])
  if(length(genes)<4){
    genes<-genes[-(length(genes))]
  }
  if(length(genes)==4){
    if(nchar(genes[4])==0)
      genes<-genes[-(length(genes))]
  }
  literature_geneSet_ense[[title]]<-genes
}

#pathway binary
keyword_path_list_all<-c(keyword_path_list,literature_geneSet_ense)
#write.table(names(keyword_path_list),file="E:/project/remodeling tumor microenvironment/NMF_test/GSVA_path_item.txt",sep="\t",quote=F)
keyword_path_list_new<-keyword_path_list_all
keyword_path_list_new[["chemokine"]]<-unique(c(keyword_path_list_new[["chemokine"]],keyword_path_list_new[["chemokines"]]))
keyword_path_list_new[["cytokine"]]<-unique(c(keyword_path_list_new[["cytokine"]],keyword_path_list_new[["cytokines"]]))
keyword_path_list_new[["Epithelial"]]<-unique(c(keyword_path_list_new[["Epithelial"]],keyword_path_list_new[["epithelial"]]))
keyword_path_list_new[["Mesenchymal"]]<-unique(c(keyword_path_list_new[["Mesenchymal"]],keyword_path_list_new[["mesenchymal"]]))
keyword_path_list_new[["mTOR"]]<-unique(c(keyword_path_list_new[["mTOR"]],keyword_path_list_new[["MTOR"]]))
keyword_path_list_new[["vessel"]]<-unique(c(keyword_path_list_new[["vessel"]],keyword_path_list_new[["vessels"]]))
keyword_path_list_new1<-list()
for(key in names(keyword_path_list_new)){
    if(!(key%in%c("chemokines","cytokines","epithelial","mesenchymal","MTOR","vessels")))
    keyword_path_list_new1[[key]]<-keyword_path_list_new[[key]]
}
#length(keyword_path_list_new1)
keyword_path_list_all<-keyword_path_list_new1

#keyword_path_list_all<-literature_geneSet_ense
pathway_nes_binary<-c()
path_ensembl_name<-c()
pahtway_nes_all<-rbind(pahtway_nes,pahtway_nes_liter)
for(key in names(keyword_path_list_all)){
    path<-intersect(keyword_path_list_all[[key]],row.names(pahtway_nes_all))
    if(length(path)>2){
     # path_ensembl_name<-c(path_ensembl_name,key)
      tmp<-apply(pahtway_nes_all[path,],2,function(x) {x[which.max(x)]=1;x[which(x!=1)]=0;return(x)})
      row.names(tmp)<-paste(key,row.names(tmp),sep="|")
      pathway_nes_binary<-rbind(pathway_nes_binary,tmp)
    }
}
apply(pathway_nes_binary,1,function(x) sum(x))->pathway_result_static

clust_function(pathway_nes_binary,"TME_cluster","complete")->TME_clust_complete

TME_clust_complete$n_clust
table(TME_clust_complete$clust)

#which(TME_clust$clust==1)
#Sample clust assignment
#Sample score normalization 
# pahtway_nes_all<-pathway_fra_normalization

path<-as.vector(unlist(sapply(row.names(pathway_nes_binary),function(x) as.vector(unlist(strsplit(x,"\\|")))[2])))
names(path)<-row.names(pathway_nes_binary)
#nes_score<-pathway_nes_norm[path,]
nes_score_norm<-t(pahtway_nes_all[path,])
colnames(nes_score_norm)<-names(path)
score_var<-apply(nes_score_norm,2,function(x) var(x))
#nes_score_norm<-t(apply(nes_score_norm,1,function(x) x/sum(x)))
clust_abun<-c()
clust_item<-TME_clust_complete$clust
for(i in seq(1,max(clust_item))){
  ft<-names(which(clust_item==i))
  tmp<-apply(nes_score_norm[,ft],1,function(x) {mean(x)})
  clust_abun<-rbind(clust_abun,tmp)
}
row.names(clust_abun)<-paste("TME_clust_",seq(1,max(clust_item)),sep="")
#clust_abun<-apply(clust_abun,2,function(x) x/sum(x))
TME_clust_assign<-apply(clust_abun,2,function(x) names(which.max(x)))
TME_clust_assign_final<-TME_clust_assign
wilcox_re<-c()
for(i in seq(1,length(TME_clust_assign))){
  clust<-as.vector(unlist(strsplit(TME_clust_assign[i],"_")))[3]
  ft<-names(which(clust_item==clust))
  sam<-colnames(TME_clust_assign)[i]
  score<-nes_score_norm[i,ft]
  score_other<-nes_score_norm[i,setdiff(colnames(nes_score_norm),ft)]
  wilcox.test(as.numeric(score,score_other))->wilcox_va
  if((wilcox_va$p.value>0.15)){ #|| (!(names(which.max(nes_score_norm[i,]))%in%ft))){
    TME_clust_assign_final[i]<-"NA"
  }
}
final_clust<-names(table(TME_clust_assign_final))[-c(1,5)]

final_clust_path<-list()

included_sample<-c()
included_path<-c()
clust_item_final<-c()
for(i in final_clust){
  clust<-as.vector(unlist(strsplit(i,"_")))[3]
  ft<-names(which(clust_item==clust))
  clust_item_final<-c(clust_item_final,clust_item[which(clust_item==clust)])
  final_clust_path[[i]]=ft
  included_sample<-c(included_sample,names(which(TME_clust_assign_final==i)))
  included_path<-c(included_path,ft)
}
clust_abun_heatmap<-t(nes_score_norm[included_sample,included_path])
clust_item<-clust_item_final


ranks <- 2:15
NMF_re <- nmf(clust_abun_heatmap,ranks, nrun=50)
coph <- NMF_re$measures$cophenetic
#plot(2:6,coph, type="b", col="purple")
coph_diff <- NULL
for (i in 2:length(coph)) 
{
  coph_diff <- c(coph_diff, coph[i-1]-coph[i])
}
k.best <- which.max(coph_diff)+1

df_filter<-clust_abun[final_clust,colnames(clust_abun_heatmap)]
df_filter<-apply(df_filter,2,function(x) x/sum(x))
clust_boxplot_df<-data.frame(as.vector(unlist(df_filter)))
clust_boxplot_df$Clust<-rep(row.names(df_filter),times=ncol(df_filter))
clust_boxplot_df$Subtype<-rep(as.vector(unlist(TME_clust_assign_final[colnames(df_filter)])),each=nrow(df_filter))


TME_distribution<-data.frame(table(TME_clust_assign_final[colnames(df_filter)]))
TME_distribution$Group<-as.vector(unlist(sapply(as.vector(unlist(TME_distribution$Var1)),function(x) paste("TME",strsplit(x,"_")[[1]][3],sep=""))))
TME_distribution<-TME_distribution[,-1]
colnames(TME_distribution)<-c("Freq","Subtype")
TME_distribution$Subtype[which(TME_distribution$Subtype=="TME10")]="TME1"
TME_distribution$Subtype[which(TME_distribution$Subtype=="TME2")]="TME3"
TME_distribution$Subtype[which(TME_distribution$Subtype=="TME12")]="TME2"
TME_distribution$Subtype[which(TME_distribution$Subtype=="TME8")]="TME7"
Subtype<-TME_distribution$Subtype


clust_function<-function(binary_feature,fig_name,method_choose){
  jacard_mt_fun(binary_feature)->Jaccard_mt_update
  hc = hclust(as.dist(1-Jaccard_mt_update), method = method_choose)
  n_clust = choose_clusters(Jaccard_mt_update, "jaccard", range = 2:(nrow(Jaccard_mt_update) - 1))
  #print(n_clust)
  #n_clust = 5
  clust =  hclusCut(Jaccard_mt_update, n_clust)$cluster
  new_ft<-c()
  #n_index<-intersect(names(which(table(clust)>2)),names(which(table(clust)<11)))
  n_index<-names(which(table(clust)>2))
  for(i in n_index){
    new_ft<-c(new_ft,(names(which(clust==i))))
  }
  #Jaccard_mt_update<-Jaccard_mt_update[row.names(pathway_nes_binary_filter),row.names(pathway_nes_binary_filter)]
  Jaccard_mt_update<-Jaccard_mt_update[new_ft,new_ft]
  hc = hclust(as.dist(1-Jaccard_mt_update), method = method_choose)
  n_clust = choose_clusters(Jaccard_mt_update, fig_name, range = 2:(nrow(Jaccard_mt_update) - 1))
  clust =  hclusCut(Jaccard_mt_update, n_clust)$cluster
  return(list(clust_detail=clust,n_clust=n_clust))
}

phyper_fun<-function(x,y){
  stat<-x+y
  all<-length(x)
  a<-sum(x)
  b<-sum(y)
  c<-all-max(a,b)
  d<-length(which(stat==2))
  pva<-1-phyper(d-1,max(a,b),c,min(a,b))
  return(pva)
}
jacard_mt_fun<-function(binary_feature){
  Jaccard_mt<-c()
  for(i in seq(1,nrow(binary_feature))){
    apply(binary_feature,1,function(x) jaccard::jaccard(as.numeric(x),as.numeric(binary_feature[i,])))->re
    Jaccard_mt<-rbind(Jaccard_mt,re)
    #Jaccard_mt_update[i,which(re>0.01)]=0
  }
  row.names(Jaccard_mt)<-row.names(binary_feature)
  Jaccard_mt_update<-Jaccard_mt
  for(i in seq(1,nrow(binary_feature))){
    apply(binary_feature,1,function(x) phyper_fun(x,binary_feature[i,]))->re
    Jaccard_mt_update[i,which(re>0.01)]=0
    #Jaccard_mt_update[which(re>0.01,i)]=0
  }
  return(Jaccard_mt_update)
}
hclusCut <- function(x, k, ...) list(cluster = cutree(hclust(as.dist(1-x), method = "average", ...), k=k))
choose_clusters <- function(data, name, range = 2:10)
{
        #gap = clusGap(jaccard, hclusCut, 50, B = 100)
        silh <- data.frame(K = range, Silhouette = sapply(range, function(k){
                sil <<- silhouette(hclusCut(data, k)$cluster, as.dist(1-data))
                tmp <<- summary(sil)
                tmp$avg.width
        }))

        g2 <- ggplot(silh, aes(x = K, y = Silhouette)) +
                geom_point() +
                #ylab("Sparseness (coef)") +
                geom_line() +
                geom_vline(xintercept = silh[which.max(silh$Silhouette), 1], lty = 2, colour = "red") +
                theme_bw() +
                theme(panel.grid = element_blank()) +
                theme(aspect.ratio = 1)
                #facet_wrap(CellType, ncol = 4) +

        pdf(file.path("result/", paste0("nclusters_", name, ".pdf")), width = 4, height = 4)
        plot(g2)
        tmp = dev.off()
        png(file.path("result/", paste0("nclusters_", name, ".png")), width = 4, height = 4, res = 300, units = "in")
        plot(g2)
        tmp = dev.off()
        silh[which.max(silh$Silhouette), 1]
}





