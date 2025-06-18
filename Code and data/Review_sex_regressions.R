####ProcD regressions of shape ~ size + sex

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

sex.info <- list()

for ( i in taxa_order){
  
  take = i
  
  size.tmp <- log10(GPA$Csize [grep(take,landmark_data$Name)])
  tmp <- !duplicated(dimnames(coords)[[3]]  [grep(take,landmark_data$Name)])
  
  #sex info
  specimens <- sapply(strsplit(names(size.tmp),'_'), function(x) x[4] )
  sex <- subset(landmark_data,Name==take)
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
  
  sex.info[[i]] <- data.frame(specimen=names(size.tmp), sex=myorder)
  
  
}

sex.info

taxa.test <- c('Astrochelys_radiata','Centrochelys_sulcata',
               'Chelodina_oblonga','Chelus_fimbriata',
               'Cuora_flavomarginata',
               'Dermatemys_mawii',
               'Emydura_macquarii',
               'Heosemys_grandis',
               'Heosemys_annandalii','Kinosternon_creaseri','Leucocephalon_yuwonoi','Manouria_emys',
               'Mauremys_sinensis','Mauremys_rivulata','Orlitia_borneensis',
               'Pelomedusa_subrufa','Phrynops_hilarii',
               'Podocnemis_vogli','Stigmochelys_pardalis',
               'Trachemys_ornata',
               'Trachemys_scripta')




test.list <- list()

for ( i in taxa.test){
  
  dat <- sex.info[[i]]
  use <- which(dat$sex !='uncertain')
    use <- dat$specimen[use]
  
  size.tmp <- log10(GPA$Csize [grep(i,landmark_data$Name)])
  size.tmp <- size.tmp[use]
  sex.tmp <- dat$sex[dat$sex!='uncertain']
  tmp <- !duplicated(use)
  
  proc.compl <-  procD.lm(GPA$coords[,,use[tmp] ] ~ size.tmp[use[tmp]] + sex.tmp[tmp] ,
                          SS.type = 'II')
  
  #sex.tmp[tmp]
  print(i)
  
  test.list[[i]] <- summary(proc.compl)$table
  
  
}

test.list

sex.sign <- sapply(test.list, function(x) x$`Pr(>F)`[2]<0.05 )
size.sign <- sapply(test.list, function(x) x$`Pr(>F)`[1]<0.05 )


pdf('size_sex_effects_comparisons.pdf',width = 7, height = 7, useDingbats = F)

par(mfrow=c(2,2))

#Differece between Z-scores of 'size' and 'sex'

barX <- barplot(rev(sapply ( test.list , function(x) x[,'Z'][1] - abs(x[,'Z'][2]) )),
        cex.names=0.25 , xlab = 'Z(size) - Z(sex)', 
        col = rev(c('grey90','red2')[sex.sign+1]) ,
        xlim=c(-2.5,3.5), horiz = T, 
        border = rev(c('grey90','black')[size.sign+1]), 
        names=NA, ylim=c(0.5,30))

sapply(1:length(taxa.test), 
       function(x) text(y=barX[x,]+0,x=-2.65, 
                        lab=rev(taxa.test)[x],
                        cex=0.5, pos=4 )  )

legend('topright',fill='grey90',title = 'Size statistically significant', 
       legend = c('no','yes'), border=c('grey90','black'), bty='n',ncol=2,cex=0.6)

legend(x=1.45,y=28,fill=c('grey90','red2'),title = 'Sex statistically significant', 
       legend = c('no','yes'), border=NA, bty='n',ncol=2,cex=0.6)


#Variation of Z-scores of 'size' vs. 'sex'
boxplot(do.call(rbind,lapply(test.list, function(x) x$Z[1:2] )),
        names=c('size','sex'), ylab='Z (effect-size)')

dev.off()


supp_table_3 <- lapply(lapply(test.list, function(x) x[1:2,4:7] ), function(x) cbind(x[1,],x[2,]) )
supp_table_3 <- do.call(rbind,supp_table_3)
write.csv(supp_table_3,file='Supplementary Table 3 (size-sex results).csv')

##Fancier plot

pdf('size_sex_effects_comparisons_Final.pdf',width = 7, height = 7, useDingbats = F)


par(mfcol=c(2,2))

#Differece between Z-scores of 'size' and 'sex'

barX <- barplot(rev(sapply ( test.list , function(x) x[,'Z'][1] - abs(x[,'Z'][2]) )),
                cex.names=0.25 , xlab = 'Z(size) - Z(sex)', 
                col = rev(c('grey90','turquoise')[sex.sign+1]) ,
                xlim=c(-1.5,3.5), horiz = T, ylab='species',
                border = rev(c('grey90','black')[size.sign+1]), 
                names=NA, ylim=c(0.5,32),cex.axis=0.8)

sapply(1:length(taxa.test), 
       function(x) text(y=barX[x,]+0,x=-1.55, 
                        lab=rev(1:length(taxa.test))[ x],
                        cex=0.8, pos=4 )  )

legend(x=1.45,y=33,fill='grey90',title = 'Size statistically significant', 
       legend = c('no','yes'), border=c('grey90','black'), bty='n',ncol=2,cex=0.6)

legend(x=1.45,y=29.5,fill=c('grey90','turquoise'),title = 'Sex statistically significant', 
       legend = c('no','yes'), border=NA, bty='n',ncol=2,cex=0.6)


#Variation of Z-scores of 'size' vs. 'sex'

plot('n',xlim=c(0.5,2.5),ylim=c(-1,5),ylab='Z (effect-size)',
     xaxt='n',xlab=NA,cex.axis=0.8)
axis(1,at=1:2,labels=c('size','sex (F/M)'),cex.axis=0.8)

vioplot(do.call(rbind,lapply(test.list, function(x) x$Z[1:2] )),
        xlab=NA,ylab=NA,add=T,col=adjustcolor('grey',alpha.f = 0.2),
        rectCol=adjustcolor('grey20',alpha.f = 0.3),
        colMed=adjustcolor('grey20',alpha.f = 0.3),lineCol=adjustcolor('grey20',alpha.f = 0.2),
        border=F)

for ( i in 1:2 ){
  
  cols <- c('grey60','turquoise')
  
  stripchart(do.call(rbind,lapply(test.list, function(x) x$Z[1:2] ))[,i], 
             col = adjustcolor(cols[i],alpha.f=0.5),#col='grey70',
             method = 'jitter',add=T,at=i,lwd=0.5,
             pch=16,vertical = T,cex=1.2)
  
}

vioplot(do.call(rbind,lapply(test.list, function(x) x$Z[1:2] )),
        xlab=NA,ylab=NA,add=T,col=adjustcolor('grey',alpha.f = 0.1),
        rectCol=adjustcolor('grey20',alpha.f = 0.4),
        colMed=adjustcolor('grey20',alpha.f = 0.4),lineCol=adjustcolor('grey20',alpha.f = 0.6),
        border=F)


dev.off()



