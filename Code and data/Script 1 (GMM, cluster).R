##################################################################################################
### Script 1 - Read landmark data, GMM analysis, cluster analysis and ontogenetic shape curves ###
##################################################################################################

library(dispRity)
library(drc)
library(geomorph)
library(Morpho)
library(vioplot)

#Read landmark dataset from Stayton's papers

landmark_data = read.table('turtle_landm-1.csv',sep = ',',header = T) #Stayton et al. 2018
updated_lm = read.table('AllTurtleDataFixed_updated.csv',sep = ',',header = T) #Stayton 2024

original_specimens <- paste0(landmark_data$Museum,'_',landmark_data$Museum.Number)
updated_specimens <- paste0(updated_lm$Museum,'_',updated_lm$Number)


#Get additional specimens from Stayton 2024 to include in the original dataset

new_specimens <- !updated_specimens %in% original_specimens
new_specimens <- updated_lm[new_specimens , -4]
  colnames(new_specimens) <- colnames(landmark_data)

#Bind two datasets
landmark_data <- rbind(new_specimens,landmark_data)

#Read additional specimens

new.dir <- paste0(getwd(),'/Stayton landmarks/')

newLMs <-  list.files(new.dir, pattern='-files')

newLM_list <- list()

for ( i in 1:length(newLMs)){
  
  dir.tmp <- paste0(new.dir,newLMs[i],'/')
  not <- c('.am','.icol','.ply') #add here (and in the next line) file extensions to be ignored
  LM.tmp <- list.files(dir.tmp)[!grepl(not[1],x=list.files(dir.tmp)) & !grepl(not[2],x=list.files(dir.tmp)) & !grepl(not[3],x=list.files(dir.tmp))]
  
  from <- '@1'
  list.tmp <-  lapply(LM.tmp, function(x) readLines(paste0(dir.tmp,x)))
  
  for ( j in 1:length(list.tmp)){
    if (length(na.omit(list.tmp[[j]][which(list.tmp[[j]]==from)+1:length(list.tmp[[j]])])) == 0 ){
      list.tmp[[j]] <- matrix(NA,
                              ncol=3,byrow = F)
    } else {
      list.tmp[[j]] <- na.omit(list.tmp[[j]][which(list.tmp[[j]]==from)+1:length(list.tmp[[j]])])
      list.tmp[[j]] <- list.tmp[[j]] [-length(list.tmp[[j]])]  
      list.tmp[[j]] <- matrix(as.numeric(do.call(rbind,strsplit( list.tmp[[j]],' '))),
                              ncol=3,byrow = F)}
  }
  
  newLM_list[[i]] <- list.tmp
  
  
  if ( !length(strsplit(newLMs[i],'-')[[1]]) > 2 ){
    names(newLM_list)[[i]] <- strsplit(newLMs[i],'-')[[1]][1]
  } else names(newLM_list)[[i]] <- paste( strsplit(newLMs[i],'-')[[1]][-length(strsplit(newLMs[i],'-')[[1]])],collapse = '_')
  
  names(newLM_list[[i]]) <-  LM.tmp
  
  print(names(newLM_list)[[i]])
  
}

new_turtles <- newLM_list

right_order <- paste0('Landmark_',1:length(new_turtles[[1]]))

new_turtles <-  lapply(new_turtles,function(x) x[right_order] )


new_turtles <- lapply(new_turtles,function(x) do.call(rbind,x) )

new_turtles_LMs <- array(NA,dim= c(nrow(new_turtles[[1]]),3,length(new_turtles) ),
      dimnames = list(c(1:nrow(new_turtles[[1]])),c('x','y','z'),names(new_turtles)))

for ( i in 1:length(new_turtles)){
  new_turtles_LMs[,,i] <- new_turtles[[i]]
    
}

new_turtles_LMs <- two.d.array(new_turtles_LMs)
  rownames(new_turtles_LMs) <- NULL

new_info <- read.csv('new_info.csv',header = T, row.names = NULL)  

new_turtles_final <- cbind(new_info,new_turtles_LMs)
colnames(new_turtles_final) <- colnames(landmark_data)
rownames(new_turtles_final) <- as.character(nrow(landmark_data)+(1:length(newLM_list)))

#Rescale some specimens because they are not in the same unit of others

#Chelonia_mydas_Sketchfab_72A
id <- grep('72A',new_turtles_final$Museum.Number)
new_turtles_final[id ,-c(1:7) ] <- new_turtles_final[id ,-c(1:7) ]*1000

#Chelonia_mydas_NMB_152
id <- grep('152',new_turtles_final$Museum.Number)
new_turtles_final[id ,-c(1:7) ] <- new_turtles_final[id ,-c(1:7) ]*100

#Trachemys_scripta_Sketchfab_unnumbered
id <- which(new_turtles_final$Name=='Trachemys_scripta' & new_turtles_final$Museum=='Sketchfab' & new_turtles_final$Museum.Number=='unnumbered')
new_turtles_final[id ,-c(1:7) ] <- new_turtles_final[id ,-c(1:7) ]*25.4 #originally in inches

#Trachemys_scripta_Sketchfab_unnumbered2
id <- which(new_turtles_final$Name=='Trachemys_scripta' & new_turtles_final$Museum=='Sketchfab' & new_turtles_final$Museum.Number=='unnumbered2')
new_turtles_final[id ,-c(1:7) ] <- new_turtles_final[id ,-c(1:7) ]*25.4 #originally in inches

#Chelus_fimbriata_UF-herp_117469
id <- grep('117469',new_turtles_final$Museum.Number)
new_turtles_final[id ,-c(1:7) ] <- new_turtles_final[id ,-c(1:7) ] /1000


#######

#Assing length data to specimen with length ==0 (Trachemys gaigeae KU 51203)

mtmp <- matrix(landmark_data[landmark_data$Length==0,colnames(landmark_data) %in% c('X1','Y1','Z1','X12','Y12','Z12')],
       ncol=3,byrow = T)

landmark_data[landmark_data$Length==0,'Length'] <- round(c(dist(mtmp)))

#Delete specimen with inconsistent landmark placement (Macrochelys temminckii CM 96008)
landmark_data <- landmark_data[-grep('96008',landmark_data$Museum.Number),]

#Delete largest Cyclemys oldhamii specimen (USNM 328006) because it is outside the range of the species (TTWG 2021)
landmark_data <- landmark_data[-grep('328006',landmark_data$Museum.Number),]

#Delete largest Podocnemis vogli specimen (FMNH 73416) because it is outside the range of the species (TTWG 2021)
landmark_data <- landmark_data[-grep('73416',landmark_data$Museum.Number),]

#Delete largest Pelomedusa subrufa specimen (USNM 196497) because it is outside the range of the species (TTWG 2021)
landmark_data <- landmark_data[-grep('196497',landmark_data$Museum.Number),]


#Bind original dataset with sea turtle new specimens
landmark_data <- rbind(new_turtles_final,landmark_data)

coords = landmark_data[,c(1,8:ncol(landmark_data))]

#make 3D array  
coords <- geomorph::arrayspecs(coords[,-1],p = ncol(coords[,-1])/3,k=3)
dimnames(coords)[[3]] <- paste0(landmark_data$Name,'_',landmark_data$Museum,'_',landmark_data$Museum.Number)

#Generalized Procrustes Alignment

GPA = gpagen(coords,ProcD = F)

#which landmark connects to which
lines_list = list(
  c(1,13,37,38,39,40,24,12),
  c(1:12),
  c(13,14,25,15:17,26,18,19,27,20,21,28,22,29,23,24),
  c(25,30:36),
  c(2,14),c(3,15),c(4,16),c(5,17),c(6,18),c(7,19),c(8,20),c(9,21),c(10,22),c(11,23),
  c(30,37),c(32,38),c(34,39),c(36,40),
  c(52,43:41),c(53,45:48),c(43:45),
  c(41,49,50,51,48),
  c(49,44),c(50,45),c(51,46),
  c(31,26),c(33,27),c(35,28),c(36,29)
)

#GPA shape consensus
plot3d(GPA$consensus , size=1.1,col='lightgrey',lit=F ,add=F,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$consensus [lines_list[[i]],] , lwd=1, col='lightgrey') }

#text3d(GPA$consensus[1:40,], texts = 1:40, pos=3,cex=1)
#text3d(GPA$consensus[-c(1:40),], texts = 41:nrow(GPA$consensus), pos=3,cex=1)


#Principal Component Analysis

PCA <- gm.prcomp(GPA$coords)

plot(PCA$x,pch=c(21,22)[as.numeric(as.factor(landmark_data$Habitat))],
     col='grey90', 
     bg=sapply( c('blue','red'),adjustcolor, alpha.f=0.4)[as.numeric(as.factor(landmark_data$Habitat))]   )

dev.off()

#Summarize: get ranges of species with at least 8 specimens for downstream exploratory analysis

summary_range <- cbind(sapply(lapply(sapply(names(which(table(landmark_data$Name)>=8)), function(x) landmark_data[ grep(x,landmark_data$Name),'Length'] ),range), function(x) x[2]/x[1] ),
                       table(landmark_data$Name)[which(table(landmark_data$Name)>=8)])
colnames(summary_range) <- c('range_size','N')
summary_range <- summary_range[ order(summary_range[,1],decreasing = T),]
head(summary_range)

summary_range[,1] <- round(summary_range[,1],1)

#Take species which have a range size of at least 3-fold
#(i.e., largest specimen is at least 3x larger than the smallest)

to_take <- rownames(summary_range)[which(summary_range[,1]>=3)]
to_take <- to_take[!to_take %in% c('Graptemys_barbouri','Graptemys_pulchra',
                                   'Hardella_thurjii','Elseya_irwini')]

summary_range <- summary_range[to_take,]

CAC.list <- list()
PredLine.list <- list()
stats.list <- list()
clust.list <- list()
clust.pred <- list()
groups.list <- list()
groups.pred <- list()

for (i in 1:nrow(summary_range)){
  
#taxon
    take <- rownames(summary_range)[i]

#Clustering based on size
#size info from table
size.tmp <- log10(GPA$Csize [grep(take,landmark_data$Name)])
names(size.tmp) <- dimnames(coords)[[3]]  [grep(take,landmark_data$Name)]
tmp <- !duplicated(dimnames(coords)[[3]]  [grep(take,landmark_data$Name)])
real.size <- log10(landmark_data$Length [grep(take,landmark_data$Name)])
  
  CAC.list[[i]] <- procD.lm(GPA$coords[,,(grep(take,landmark_data$Name)[tmp]) ] ~ size.tmp[tmp])
  stats.list[[i]] <- round(unlist(c(summary(CAC.list[[i]])$table[c('Rsq','Z','Pr(>F)')][1,])),4)
  PredLine.list[[i]] <- plotAllometry(CAC.list[[i]],size = size.tmp[tmp],logsz = F,method = 'PredLine')$PredLine
  CAC.list[[i]] <- plotAllometry(CAC.list[[i]],size = size.tmp[tmp],logsz = F,method = 'CAC')$CAC
  
  
  smallest <- which.min(size.tmp[tmp])
  largest <- which.max(size.tmp[tmp])
  
  if(PredLine.list[[i]][smallest] > PredLine.list[[i]][largest]){
    PredLine.list[[i]] <- cbind(PredLine.list[[i]]*(-1) , real.size[tmp])}
  else 
    PredLine.list[[i]] <- cbind(PredLine.list[[i]] , real.size[tmp])
  
 # if(CAC.list[[i]][smallest] > CAC.list[[i]][largest]){
  #  CAC.list[[i]] <- cbind(CAC.list[[i]]*(-1) , real.size[tmp])}
  #else 
    CAC.list[[i]] <- cbind(CAC.list[[i]] , real.size[tmp])
  
    dev.off()
  
  clust.list[[i]] <- hclust(dist(cbind(CAC.list[[i]][,1:2])),method = 'ward.D')
  groups.list[[i]] <- as.factor(cutree(clust.list[[i]],k=3))
  
  clust.pred[[i]] <- hclust(dist(cbind(PredLine.list[[i]][,1:2])),method = 'ward.D')
  groups.pred[[i]] <- as.factor(cutree(clust.pred[[i]],k=3))
  
  print(take)
}

names(CAC.list) <- names(PredLine.list) <- names(stats.list) <- names(clust.list) <- names(groups.list) <- names(clust.pred) <- names(groups.pred) <- rownames(summary_range)

#change group names to 'small', 'intermediate' and 'large'
for ( i in 1:length(groups.list)){
  maxx <- names(which.max(CAC.list[[i]][,2]))
  minn <- names(which.min(CAC.list[[i]][,2]))
  
  levels(groups.list[[i]])[levels(groups.list[[i]]) %in% groups.list[[i]][maxx]] <- 'large'
  levels(groups.list[[i]])[levels(groups.list[[i]]) %in% groups.list[[i]][minn]] <- 'small'
  
  interm <- names(groups.list[[i]][groups.list[[i]] !='large' & groups.list[[i]] !='small' ])
  levels(groups.list[[i]])[levels(groups.list[[i]]) %in% groups.list[[i]][interm]] <- 'interm'
  
  groups.list[[i]] <- factor(groups.list[[i]],levels = c('small','interm','large'))
  groups.list[[i]] <- droplevels(groups.list[[i]])
  
  
  levels(groups.pred[[i]])[levels(groups.pred[[i]]) %in% groups.pred[[i]][maxx]] <- 'large'
  levels(groups.pred[[i]])[levels(groups.pred[[i]]) %in% groups.pred[[i]][minn]] <- 'small'
  
  interm <- names(groups.pred[[i]][groups.pred[[i]] !='large' & groups.pred[[i]] !='small' ])
  levels(groups.pred[[i]])[levels(groups.pred[[i]]) %in% groups.pred[[i]][interm]] <- 'interm'
  
  groups.pred[[i]] <- factor(groups.pred[[i]],levels = c('small','interm','large'))
  groups.pred[[i]] <- droplevels(groups.pred[[i]])
  
  
}


#For some reason, some species CAC values are not positively correlated with SCL. Invert
revert <- c('Chelodina_oblonga','Elseya_dentata','Emydura_macquarii')

for ( i in revert){
  CAC.list[[i]][,1] <- CAC.list[[i]][,1]*(-1)
  PredLine.list[[i]][,1] <- PredLine.list[[i]][,1]*(-1)
  }


supp_table_1 <- cbind('sample size'=summary_range[,2], do.call(rbind,stats.list))
supp_table_1 <- supp_table_1[sort(rownames(supp_table_1)),]
write.csv(supp_table_1,file='Supplementary Table 1.csv')

#PDF 1
#Plot significant ones with large individuals sampled
taxa_order <- c('Aldabrachelys_gigantea','Astrochelys_radiata','Centrochelys_sulcata',
                'Chelonia_mydas','Chelus_fimbriata','Cyclemys_oldhamii',
                'Dermatemys_mawii','Elseya_novaeguineae','Heosemys_grandis',
                'Heosemys_annandalii','Kinosternon_creaseri','Leucocephalon_yuwonoi',
                'Manouria_emys','Mauremys_rivulata','Mauremys_sinensis',
                'Orlitia_borneensis','Pelusios_rhodesianus','Podocnemis_expansa',
                'Podocnemis_sextuberculata','Podocnemis_unifilis',
                'Podocnemis_vogli','Stigmochelys_pardalis','Trachemys_ornata',
                'Trachemys_scripta','Trachemys_yaquia',
                'Chelodina_oblonga','Chelodina_parkeri','Cuora_flavomarginata',
                'Elseya_dentata',
                'Emydura_macquarii','Pelomedusa_subrufa','Phrynops_hilarii',
                'Trachemys_callirostris')
taxa_order <- sort(taxa_order)

pdf('1_CAC_plots_significant_large_specimens_sampled.pdf',width = 7, height = 7, useDingbats = F)

par(mfrow=c(3,3))
for ( i in taxa_order){
  
  mydata <- as.data.frame(CAC.list[[i]])
  mydata$V2 <- 10^mydata$V2
  lm.tmp <- drm(V1~V2,data=mydata,fct=G.4())
  
  
  specimens <- sapply(strsplit(rownames(mydata),'_'), function(x) x[4] )
  sex <- subset(landmark_data,Name==i)
  sex <- sex$Sex[ sex$Museum.Number %in% specimens ]
  females <- grep('F',sex)
    females <- setNames(rep('F',length(females)),females)
  males <- grep('M',sex)
    males <- setNames(rep('M',length(males)),males)
  unknown <- which(sex=='?')
    unknown <- setNames(rep('uncertain',length(unknown)),unknown)
  juvenile <- which(sex=='J')
    juvenile <- setNames(rep('uncertain',length(juvenile)),juvenile)
  
  myorder <- c(females,males,unknown,juvenile)[as.character(1:length(sex))]
  
  mypch <- setNames(c(16,17,18),c('F','M','uncertain'))
  mypch <- mypch[names(mypch) %in% unique(myorder)]
  
  plot(x=10^CAC.list[[i]][,2],y=CAC.list[[i]][,1], xlab='Centroid size',ylab='Predicted shape',
       pch=mypch [as.numeric(as.factor(myorder))] , 
       col=c('#ccebc5' , '#4eb3d3' , '#084081') [as.numeric(groups.list[[i]])] ,
       main=i, bty='l')
  
  
  new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l')
  
  
  legend('topleft',cex=0.7,legend=levels(groups.list[[i]]),
         fill=c('#ccebc5' , '#4eb3d3' , '#084081'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 

}

dev.off()

#PredLine
par(mfrow=c(3,3))
for ( i in taxa_order){
  
  mydata <- as.data.frame(PredLine.list[[i]])
  mydata$V2 <- 10^mydata$V2
  lm.tmp <- drm(V1~V2,data=mydata,fct=G.4())
  
  
  specimens <- sapply(strsplit(rownames(mydata),'_'), function(x) x[4] )
  sex <- subset(landmark_data,Name==i)
  sex <- sex$Sex[ sex$Museum.Number %in% specimens ]
  females <- grep('F',sex)
  females <- setNames(rep('F',length(females)),females)
  males <- grep('M',sex)
  males <- setNames(rep('M',length(males)),males)
  unknown <- which(sex=='?')
  unknown <- setNames(rep('uncertain',length(unknown)),unknown)
  juvenile <- which(sex=='J')
  juvenile <- setNames(rep('uncertain',length(juvenile)),juvenile)
  
  myorder <- c(females,males,unknown,juvenile)[as.character(1:length(sex))]
  
  mypch <- setNames(c(16,17,18),c('F','M','uncertain'))
  mypch <- mypch[names(mypch) %in% unique(myorder)]
  
  plot(x=10^PredLine.list[[i]][,2],y=PredLine.list[[i]][,1], xlab='Centroid size',
       ylab='Multivariate shape',
       pch=mypch [as.numeric(as.factor(myorder))] , 
       col=c('#ccebc5' , '#4eb3d3' , '#084081') [as.numeric(groups.pred[[i]])] ,
       main=i, bty='l')
  
  
  new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l')
  
  
  legend('topleft',cex=0.7,legend=levels(groups.pred[[i]]),
         fill=c('#ccebc5' , '#4eb3d3' , '#084081'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 
  
}

dev.off()







#PDF 2
#Plot significant ones, but very large individuals missing from sample
taxa_order <- c('Chelonoidis_denticulata','Chelonoidis_nigra','Chelydra_acutirostris',
                'Chelydra_serpentina','Geochelone_elegans',
                'Geoclemys_hamiltonii','Mauremys_reevesii',
                'Pelusios_gabonensis','Pelusios_niger','Pelusios_sinuatus',
                'Podocnemis_lewyana',
                'Pseudemys_gorzugi','Testudo_marginata','Trachemys_grayi',
                'Elseya_dentata','Chelodina_oblonga','Emydura_macquarii','Pelomedusa_subrufa',
                'Chelodina_rugosa','Trachemys_callirostris','Elseya_irwini',
                'Cuora_flavomarginata','Chelodina_parkeri','Phrynops_hilarii',
                'Pelusios_castaneus','Podocnemis_erythrocephala')
taxa_order <- sort(taxa_order)


pdf('2_CAC_plots_significant_full_adult_range_unsampled.pdf',width = 7, height = 7, useDingbats = F)

par(mfrow=c(3,3))
for ( i in taxa_order){
  
  mydata <- as.data.frame(CAC.list[[i]])
  mydata$V2 <- 10^mydata$V2
  lm.tmp <- drm(V1~V2,data=mydata,fct=G.4())
  
  
  specimens <- sapply(strsplit(rownames(mydata),'_'), function(x) x[4] )
  sex <- subset(landmark_data,Name==i)
  sex <- sex$Sex[ sex$Museum.Number %in% specimens ]
  females <- grep('F',sex)
  females <- setNames(rep('F',length(females)),females)
  males <- grep('M',sex)
  males <- setNames(rep('M',length(males)),males)
  unknown <- which(sex=='?')
  unknown <- setNames(rep('uncertain',length(unknown)),unknown)
  juvenile <- which(sex=='J')
  juvenile <- setNames(rep('uncertain',length(juvenile)),juvenile)
  
  myorder <- c(females,males,unknown,juvenile)[as.character(1:length(sex))]
  
  mypch <- setNames(c(16,17,18),c('F','M','uncertain'))
  mypch <- mypch[names(mypch) %in% unique(myorder)]
  
  plot(x=10^CAC.list[[i]][,2],y=CAC.list[[i]][,1], xlab='SCL',ylab='Predicted shape',
       pch=mypch [as.numeric(as.factor(myorder))] , 
       col=c('#ccebc5' , '#4eb3d3' , '#084081') [as.numeric(groups.list[[i]])] ,
       main=i, bty='l')
  
  
  new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l')
  
  
  legend('topleft',cex=0.7,legend=levels(groups.list[[i]]),
         fill=c('#ccebc5' , '#4eb3d3' , '#084081'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 
  
}

dev.off()


#PDF 3
#Plot non-significant results

pdf('3_CAC_plots_nonsignificant.pdf', width = 7, height = 7, useDingbats = F)

taxa_order <- names(CAC.list)[ which(sapply(stats.list, function(x) x[3]>0.05 ))]
taxa_order <- sort(taxa_order)

par(mfrow=c(3,3))
for ( i in taxa_order){
  
  mydata <- as.data.frame(CAC.list[[i]])
  mydata$V2 <- 10^mydata$V2
  lm.tmp <- drm(V1~V2,data=mydata,fct=G.4())
  
  
  specimens <- sapply(strsplit(rownames(mydata),'_'), function(x) x[4] )
  sex <- subset(landmark_data,Name==i)
  sex <- sex$Sex[ sex$Museum.Number %in% specimens ]
  females <- grep('F',sex)
  females <- setNames(rep('F',length(females)),females)
  males <- grep('M',sex)
  males <- setNames(rep('M',length(males)),males)
  unknown <- which(sex=='?')
  unknown <- setNames(rep('uncertain',length(unknown)),unknown)
  juvenile <- which(sex=='J')
  juvenile <- setNames(rep('uncertain',length(juvenile)),juvenile)
  
  myorder <- c(females,males,unknown,juvenile)[as.character(1:length(sex))]
  
  mypch <- setNames(c(16,17,18),c('F','M','uncertain'))
  mypch <- mypch[names(mypch) %in% unique(myorder)]
  
  plot(x=10^CAC.list[[i]][,2],y=CAC.list[[i]][,1], xlab='SCL',ylab='CAC',
       pch=mypch [as.numeric(as.factor(myorder))] , 
       col=c('#ccebc5' , '#4eb3d3' , '#084081') [as.numeric(groups.list[[i]])] ,
       main=i, bty='l')
  

  new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l')
  
  
  legend('topleft',cex=0.7,legend=levels(groups.list[[i]]),
         fill=c('#ccebc5' , '#4eb3d3' , '#084081'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 
  
}


dev.off()




###############
#Shape deformations (individual species)


take='Cuora_flavomarginata'

quarts <- quantile(CAC.list[[take]][,1],probs=c(0.025,0.975))

PCwarps = shape.predictor(GPA$coords[,, rownames(CAC.list[[take]]) ],
                          x=CAC.list[[take]] [,1],Intercept = T,
                          min=min(quarts),max=max(quarts))

#Simple plot

open3d()
par3d(windowRect = c(0,0,500,500))
Sys.sleep(1)
mfrow3d(nr = 1, nc = 3, byrow = TRUE, sharedMouse = TRUE)


choose='min'
plot3d(PCwarps[[choose]] , size=1.2,col='black' ,add=F,type='s',
       box=F,axes=F,lit=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(PCwarps[[choose]] [lines_list[[i]],] , lwd=1, col='black') }
#title3d(choose)

choose='max'
plot3d(PCwarps[[choose]] , size=1.2,col='black' ,add=F,type='s',
       box=F,axes=F,lit=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(PCwarps[[choose]] [lines_list[[i]],] , lwd=1, col='black') }
#title3d(choose)





#The code below plot the normalized Euclidean distance between
#the minimum and maximum CAC-associated shape changes, colour-coded as a gradient
#The redder the colour, greater morphological difference between two points


#### Custom function for calculation of Euclidean distance ####

Edist <- function ( x , Y ) { ( sum( ( x - Y ) ^ 2 ) ) ^ 0.5 }

#plot
#open3d()


# define the colours you want
point.distance.scale <- colorRamp(c("lightgrey" ,'blue')) 

point.distances <- c()
for ( i in 1:nrow(PCwarps$max) ){
  
  # here you calculate the Euclidean distance between the points in each of your shape matrices ('shape1' and 'shape2')
  point.distances[i] <- Edist(PCwarps$min[i,] , PCwarps$max[i,]) 
}	

#point.distances

# normalise the distances so they range from 0 to 1
point.distances.norm <- (point.distances - min(point.distances)) / max ( point.distances - min(point.distances))

# and then you're able to create a colorRamp that goes from 0 (grey) to 1 (red)
point.colours <- point.distance.scale(point.distances.norm)


choose = "min"

#min
plot3d(PCwarps[[choose]] , size=1.1,col='lightgrey',lit=F ,add=F,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(PCwarps[[choose]] [lines_list[[i]],] , lwd=1, col='lightgrey') }


plot3d(PCwarps$max, size = 1.4, lit=F ,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso',add=T,
       col=rgb(point.colours,maxColorValue = 255) )
for ( i in 1:length(lines_list)){
  lines3d(PCwarps$max [lines_list[[i]],] , lwd=1, col=rgb(point.colours,maxColorValue = 255) ) }



#######



#######

