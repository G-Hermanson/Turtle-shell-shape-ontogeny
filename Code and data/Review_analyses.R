#Review 1

#####################


#Multivariate approach of inferring size thresholds

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

#Maximum sizes from TTWG (2021), but these reflect the largest sex represented in these species-specific datasets.
#In some cases (e.g., Kinosternon, Chelus, Manouria), max. sizes are identical between females/males
max.sizes <- c(1380,395,860,310,274,1410,437,190,
               260,479,330,229,368,
               500,346,125,278,600,244,271,800,197,
               255,408,
               1090,250,500,
               360, 655,315,
               359,328,309)
names(max.sizes) <- taxa_order


thresholds <- c()
thresholds.CAC <- c()

#Derive a maximum GPA Csize from real sizes
weird <- which(landmark_data$X1<0.1 & landmark_data$X1>0) #Remove weird Csizes from analysis

linear.size <- lm(log10(GPA$Csize)[-weird] ~ log10(landmark_data$Length)[-weird] )
#coef(linear.size)

preds.list <- list() #store shape vectors/size from PredLine estimates
preds.CAC <- list() #store shape vectors/size from CAC estimates

for ( taxon in 1:length(taxa_order)){

take = taxa_order[taxon]

coords.tmp <- GPA$coords[,,names(groups.list[[take]])]

size.tmp <- GPA$Csize[names(groups.list[[take]])]
size.tmp <- size.tmp[size.tmp>1]

coords.tmp <- coords.tmp[,,names(size.tmp)]

proctmp <- procD.lm( coords.tmp ~ log10(size.tmp)    )
#summary(proctmp)

#PredLine method

preds <- plotAllometry(proctmp,logsz = F,method = 'PredLine',size=proctmp$X[,2])$PredLine
dev.off()

#smallest <- which.min(size.tmp)
#largest <- which.max(size.tmp)

#if(preds[smallest] > preds[largest]){
 # preds <- preds*(-1)}
#else 
 # preds <- preds

max.GPA <- c(1,log10(max.sizes[take])) %*% matrix(coef(linear.size)) #Maximum GPA Csize transformed from maximum real size
max.GPA <- c(max.GPA)
sizes <- seq(min(proctmp$X[,2]) , max.GPA, length.out=1000)

shapes <- lapply ( sizes , function(x) shape.predictor(A=coords.tmp,
                                             x=preds, Intercept = T, shape=x)  ) 

shapes <- lapply(shapes,function(x) x[[1]] )

ldmk.tmp <-  array(NA,dim=c(53,3,length(shapes)),
                   dimnames = list(c(1:53),c('x','y','z'),paste0('size_',1:1000)  ))
for ( i in 1:length(shapes)){
  
  ldmk.tmp[,,i] <- shapes[[i]]
}




#Introduce some 'noise'
noise <- morphol.disparity(proctmp,print.progress = F)^0.5
for ( i in 1:length(shapes)){
  
  #noise <- seq(-noise,noise,length.out = 100)
  #noise <- sample(noise,1)
  noise <- noise/((1:1000)^0.1)
  
  ldmk.tmp[,,i] <- ldmk.tmp[,,i]+noise[i]
}
  
dists <- 1-apply(ldmk.tmp,3,ProcDist,shape1=ldmk.tmp[,,1000])
dists <- (dists-min(dists)) / (max(dists-min(dists)))

thr <- 0.85

minimum <- which(dists>=thr)[1]

thresholds[taxon] <-  10^sizes[minimum] / 10^max(sizes)


#proctmp2 <- procD.lm(ldmk.tmp ~ sizes)
#summary(proctmp2)

#Also use real maximum as the maximum in the sequence below
min.real <- min(subset(landmark_data,Name==take)$Length)
real.sizes <- seq(log10(min.real) , log10(max.sizes[take]) ,
                  length.out=1000)

preds <- plotAllometry(proctmp2,size=10^sizes, method = 'PredLine',logsz = F)$PredLine
dev.off()

preds.list[[taxon]] <- data.frame('Predicted_shape'=preds,'CS'=10^sizes, 'SCL'=10^real.sizes)
######

#CAC method
preds <- plotAllometry(proctmp,logsz = F,method = 'CAC',size=proctmp$X[,2])$CAC
dev.off()

#smallest <- which.min(size.tmp)
#largest <- which.max(size.tmp)

#if(preds[smallest] > preds[largest]){
# preds <- preds*(-1)}
#else 
# preds <- preds

max.GPA <- c(1,log10(max.sizes[take])) %*% matrix(coef(linear.size)) #Maximum GPA Csize transformed from maximum real size
max.GPA <- c(max.GPA)
sizes <- seq(min(proctmp$X[,2]) , max.GPA, length.out=1000)


shapes <- lapply ( sizes , function(x) shape.predictor(A=coords.tmp,
                                                       x=preds, Intercept = T, shape=x)  ) 

shapes <- lapply(shapes,function(x) x[[1]] )

ldmk.tmp <-  array(NA,dim=c(53,3,length(shapes)),
                   dimnames = list(c(1:53),c('x','y','z'),paste0('size_',1:1000)  ))
for ( i in 1:length(shapes)){
  
  ldmk.tmp[,,i] <- shapes[[i]]
}

#Introduce some 'noise'
noise <- morphol.disparity(proctmp,print.progress = F)^0.5
for ( i in 1:length(shapes)){
  
  #noise <- seq(-noise,noise,length.out = 100)
  #noise <- sample(noise,1)
  noise <- noise/((1:1000)^0.1)
  
  ldmk.tmp[,,i] <- ldmk.tmp[,,i]+noise[i]
}


dists <- 1-apply(ldmk.tmp,3,ProcDist,shape1=ldmk.tmp[,,1000])
dists <- (dists-min(dists)) / (max(dists-min(dists)))

thr <- 0.85

minimum <- which(dists>=thr)[1]

thresholds.CAC[taxon] <-  10^sizes[minimum] / 10^max(sizes)

#proctmp2 <- procD.lm(ldmk.tmp ~ sizes)
#summary(proctmp2)

#Also use real maximum as the maximum in the sequence below
min.real <- min(subset(landmark_data,Name==take)$Length)
real.sizes <- seq(log10(min.real) , log10(max.sizes[take]) ,
                  length.out=1000)

preds <- plotAllometry(proctmp2,size=10^sizes, method = 'CAC',logsz = F)$CAC
dev.off()

preds.CAC[[taxon]] <- data.frame('Predicted_shape'=preds,'CS'=10^sizes, 'SCL'=10^real.sizes)
######

print(take)

}

names(preds.list) <- names(preds.CAC) <- names(thresholds) <- names(thresholds.CAC) <-  taxa_order


thresholds
thresholds.CAC



###new plots?

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


#pdf('Review_CAC_new plot_large specimens sampled.pdf',width = 7,height = 7,useDingbats = F)

pdf('Review_CAC&MultiShape_large specimens sampled.pdf',width = 8,height = 8,useDingbats = F)

thresholds.Gomp.CAC <- c()
thresholds.Gomp.Pred <- c()

par(mfrow=c(4,4))
for ( i in taxa_order){
  
  tmp <- !duplicated(dimnames(coords)[[3]]  [grep(i,landmark_data$Name)])
  real.size <- log10(landmark_data$Length [grep(i,landmark_data$Name)])
  
  mydata <- as.data.frame(CAC.list[[i]])
  #mydata$V2 <- 10^mydata$V2
  mydata$V2 <- 10^(real.size[tmp])
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
  
  plot(x=mydata$V2,y=mydata$V1, xlab='Carapace length',ylab='CAC',
       pch=mypch [as.numeric(as.factor(myorder))] , 
       xlim=c(min(mydata$V2), max(preds.list[[i]][,'SCL'])  ),
       col=c('#ccebc5' , '#4eb3d3' , '#084081') [as.numeric(groups.list[[i]])] ,
       main=gsub('_',' ',i),
       bty='l'
       #,xaxt='n'
       )
  
  #ticks <- seq(min(mydata$V2), max(preds.list[[i]][,'CS']),length.out=5)
  #ticks2 <- seq(min(preds.list[[i]][,'SCL']), max(preds.list[[i]][,'SCL']),length.out=5)
  #axis(1,at=ticks,lab=ticks2)
  
  #new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
  new <- seq(min(mydata$V2), max(preds.list[[i]][,'SCL']), length.out=1000)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l')
  
  legend('topleft',cex=0.7,legend=c('S','I','L'),
         fill=c('#ccebc5' , '#4eb3d3' , '#084081'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 
  
  #bounds <- coef(lm.tmp)[3]
  #alfa=0.1
  #bounds <- bounds * c(1-alfa,1+alfa)
  bounds <- confint(lm.tmp,level = 0.9)[3,]

  segments(x0=max(preds.list[[i]][,3])*0.4,x1=max(preds.list[[i]][,3])*1.2,
           y0=bounds[1],y1=bounds[1],lty=2)
  #segments(x0=max(preds.list[[i]][,3])*0.4,x1=max(preds.list[[i]][,3])*1.2,
   #        y0=bounds[2],y1=bounds[2],lty=2)
  
  
  #real.thresh <- (thresholds[i] /max(preds.list[[i]][,'CS'])) * (max.sizes[i])
  
  #real.thresh.CAC <- thresholds.CAC[i] /max(preds.CAC[[i]][,'CS']) * (max.sizes[i])
  real.thresh.CAC <- thresholds.CAC[i]*max.sizes[i]
  
  rect(xleft = real.thresh.CAC, xright = max(preds.CAC[[i]][,'CS'])*1.2,
       ybottom =min(preds.CAC[[i]][,1])*1.2,ytop =max(preds.CAC[[i]][,1])*1.2,
       border=NA, col=adjustcolor('grey',alpha.f = 0.2))
  
  #abline(v=real.thresh, col='red') #PredLine threshold
  abline(v=real.thresh.CAC, col='red',lty=1) #CAC threshold
  
  
  #Extract Gompertz lower thresholds
  
  #Coefficients
  b <- coef(lm.tmp)[1]
  c <- coef(lm.tmp)[2]
  d <- coef(lm.tmp)[3]
  e <- coef(lm.tmp)[4]
  
  # Bound
  y <- bounds[1]

  # Get SCL value for lower bound
  x <- (1 / b) * log(-log((y - c) / (d - c))) + e
  
  thresholds.Gomp.CAC[i] <- x
  
  
#PredLine plots
  tmp <- !duplicated(dimnames(coords)[[3]]  [grep(i,landmark_data$Name)])
  real.size <- log10(landmark_data$Length [grep(i,landmark_data$Name)])
  
  mydata <- as.data.frame(PredLine.list[[i]])
  #mydata$V2 <- 10^mydata$V2
  mydata$V2 <- 10^(real.size[tmp])
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
  
  plot(x=mydata$V2,y=mydata$V1, xlab='Carapace length',ylab='Multivariate shape',
       pch=mypch [as.numeric(as.factor(myorder))] , 
       xlim=c(min(mydata$V2), max(preds.list[[i]][,'SCL'])  ),
       col=c('#ccebc5' , '#4eb3d3' , '#084081') [as.numeric(groups.pred[[i]])] ,
       main=gsub('_',' ',i),
       bty='l'
       #,xaxt='n'
  )
  
  #new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
  new <- seq(min(mydata$V2), max(preds.list[[i]][,'SCL']), length.out=1000)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l')
  
  legend('topleft',cex=0.7,legend=c('S','I','L'),
         fill=c('#ccebc5' , '#4eb3d3' , '#084081'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 
  
  #bounds <- coef(lm.tmp)[3]
  #alfa=0.1
  #bounds <- bounds * c(1-alfa,1+alfa)
  bounds <- confint(lm.tmp,level = 0.9)[3,]
  
  segments(x0=max(preds.list[[i]][,3])*0.4,x1=max(preds.list[[i]][,3])*1.2,
           y0=bounds[1],y1=bounds[1],lty=2)
  #segments(x0=max(preds.list[[i]][,3])*0.4,x1=max(preds.list[[i]][,3])*1.2,
  #        y0=bounds[2],y1=bounds[2],lty=2)
  
  
  #real.thresh <- (thresholds[i] /max(preds.list[[i]][,'CS'])) * (max.sizes[i])
  real.thresh <- thresholds[i]*max.sizes[i]
  
  rect(xleft = real.thresh, xright = max(preds.list[[i]][,'CS'])*1.2,
       ybottom =min(preds.list[[i]][,1])*1.2,ytop =max(preds.list[[i]][,1])*1.2,
       border=NA, col=adjustcolor('grey',alpha.f = 0.2))
  
  abline(v=real.thresh, col='red') #PredLine threshold
  #abline(v=real.thresh.CAC, col='red',lty=1) #CAC threshold
  
  
  #Extract Gompertz lower thresholds
  
  #Coefficients
  b <- coef(lm.tmp)[1]
  c <- coef(lm.tmp)[2]
  d <- coef(lm.tmp)[3]
  e <- coef(lm.tmp)[4]
  
  # Bound
  y <- bounds[1]
  
  # Get SCL value for lower bound
  x <- (1 / b) * log(-log((y - c) / (d - c))) + e
  
  thresholds.Gomp.Pred[i] <- x
  
  
  
}



dev.off()


#######




#Summary of thresholds
#Some Gompertz values failed to converge, so transform them to NAs
thresholds.Gomp.CAC[thresholds.Gomp.CAC<0 | is.nan(thresholds.Gomp.CAC)] <- NA
thresholds.Gomp.Pred[thresholds.Gomp.Pred<0 | is.nan(thresholds.Gomp.Pred)] <- NA



par(mfrow=c(2,2))

#CAC- simulated
hist((thresholds.CAC*100),
     xlim=c(0,150),breaks=4,
     border=NA,col=adjustcolor('orange',alpha.f = 0.4),ylim=c(0.05,12),
     main='Thresholds\n(CAC)',xlab='% of maximum recorded size')
#PredLine- simulated
hist((thresholds *100),
     xlim=c(0,150),breaks=4,
     border=NA,col=adjustcolor('purple',alpha.f = 0.4),ylim=c(0.05,12),
     main='Thresholds\n(Multivariate shape)',xlab='% of maximum recorded size')
#Gompertz-derived (CAC)
hist((thresholds.Gomp.CAC / max.sizes)*100,
     xlim=c(0,150),breaks=4,
     border=NA,col=adjustcolor('cyan4',alpha.f = 0.4),ylim=c(0.05,12),
     main='Thresholds\n(Gompertz-derived CAC)',xlab='% of maximum recorded size')

#Gompertz-derived (Multivariate shape)
hist((thresholds.Gomp.Pred / max.sizes)*100,
     xlim=c(0,150),breaks=4,
     border=NA,col=adjustcolor('cyan4',alpha.f = 0.4),ylim=c(0.05,12),
     main='Thresholds\n(Gompertz-derived Multivariate shape)',xlab='% of maximum recorded size')



threshold <- data.frame('PredLine'=thresholds,
                        'CAC'=thresholds.CAC,
                        'Gompertz_derived_CAC'=thresholds.Gomp.CAC/max.sizes,
                        'Gompertz_derived_MultiSh'=thresholds.Gomp.Pred/max.sizes,
                          'Max_sizes'=max.sizes)
apply(threshold[,1:4],2,summary)

t(apply(threshold,1,function(x) x[5]*(x[1:4]) ))

dev.off()


supp_table_2 <- cbind(threshold , 
                      t(apply(threshold,1,function(x) x[5]*(x[1:4]) )))
write.csv(supp_table_2,file='Supplementary Table 2 (thresholds).csv')
##########

