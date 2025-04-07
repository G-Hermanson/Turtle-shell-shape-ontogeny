############################################################################
### Script 2 - Disparity analysis between pre-defined ontogenetic stages ###
############################################################################


#Shape analysis

specimens_to_use <- rownames(do.call(rbind,CAC.list))

PCA_subset <- procSym(coords[,,specimens_to_use],sizeshape = F)

groups.temp <- unname(groups.list)
groups.temp <- unlist(groups.temp)
groups.temp <- groups.temp[rownames(PCA_subset$PCscores)]

landmark_data_new <- paste(landmark_data$Name,landmark_data$Museum,landmark_data$Museum.Number,sep='_')

landmark_data_new <- landmark_data[ landmark_data_new %in% specimens_to_use,]

clades <- landmark_data_new$Family
names(clades) <- paste(landmark_data_new$Name,landmark_data_new$Museum,
                       landmark_data_new$Museum.Number,sep='_')

clades <- clades[specimens_to_use]

clades[which(clades=='Podocnemidae')] <- 'Podocnemididae'

clades[which(clades=='Chelydridae' | clades=='Dermatemydidae' | clades=='Kinosternidae')] <- 'Chelydroidea'

clades[which(clades=='Pelomedusidae' | clades=='Podocnemididae')] <- 'Pelomedusoides'



############

pdf('disparity_ontogenetic_groups.pdf',width = 7,height = 7)

layout(matrix(c(1,2,3,
                1,2,3,
                4,4,4,
                4,4,4,
                0,0,0),5,byrow = T), widths = rep(1/3,3), heights = c(0.3,0.3,0.4,0.6,0.2))



#plot(PCA_subset$PCscores[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4))


#A 'small' specimens
PCA.small <- names(groups.temp)[groups.temp=='small']

PCA.small <- PCA_subset$PCscores[PCA.small,]

plot(PCA_subset$PCscores[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4))
points(PCA.small, pch=21,cex=1.35,
       col = adjustcolor('black',alpha.f=0.55),
       bg = adjustcolor('#fcc5c0',alpha.f=0.55))


#B 'intermediate' specimens
PCA.interm <- names(groups.temp)[groups.temp=='interm']

PCA.interm <- PCA_subset$PCscores[PCA.interm,]

plot(PCA_subset$PCscores[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4))
points(PCA.interm, pch=21,cex=1.35,
       col = adjustcolor('black',alpha.f=0.55),
       bg = adjustcolor('#dd3497',alpha.f=0.55))

#C 'large' specimens
PCA.large <- names(groups.temp)[groups.temp=='large']

PCA.large <- PCA_subset$PCscores[PCA.large,]

plot(PCA_subset$PCscores[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4))
points(PCA.large, pch=21,cex=1.35,
       col = adjustcolor('black',alpha.f=0.55),
       bg = adjustcolor('#7a0177',alpha.f=0.55))



#Disparity
disp_groups <- list('small'=names(groups.temp)[groups.temp=='small'],
                    'interm'=names(groups.temp)[groups.temp=='interm'],
                    'large'=names(groups.temp)[groups.temp=='large'])

disp_data <- custom.subsets(PCA_subset$PCscores,
                            group = disp_groups)
set.seed(123)
disp_data <- boot.matrix(disp_data,bootstraps = 500, rarefaction = 'min')
disp_data_range <- dispRity(disp_data,metric = c(sum,ranges))


boot_disp <- lapply(disp_data_range$disparity, function(x) x[[2]])
boot_disp <- t(do.call(rbind,boot_disp))



#plot
plot('n',xlim=c(0.5,3.5),ylim=c(4.5,6),ylab='Sum of ranges',
     xaxt='n',xlab=NA)
axis(1,at=1:3,labels=c('small','intermediate','large'))

vioplot(boot_disp,
        xlab=NA,ylab=NA,add=T,col=adjustcolor('grey',alpha.f = 0.2),
        rectCol=adjustcolor('grey20',alpha.f = 0.3),
        colMed=adjustcolor('grey20',alpha.f = 0.3),lineCol=adjustcolor('grey20',alpha.f = 0.2),
        border=F)

for ( i in 1:ncol(boot_disp) ){
  
  cols <- c('#fcc5c0','#dd3497','#7a0177')
  
  stripchart(boot_disp[,i], bg = adjustcolor(cols[i],alpha.f=0.5),col='grey70',
             method = 'jitter',add=T,at=i,lwd=0.5,
             pch=21,vertical = T,cex=1.2)
  
}


vioplot(boot_disp,
        xlab=NA,ylab=NA,add=T,col=adjustcolor('grey',alpha.f = 0.1),
        rectCol=adjustcolor('grey20',alpha.f = 0.4),
        colMed=adjustcolor('grey20',alpha.f = 0.4),lineCol=adjustcolor('grey20',alpha.f = 0.3),
        border=F)


dev.off()


#Test whether disparity 95% CIs are different

combs <- t(combn(1:ncol(boot_disp),2))

Z.comparisons <- list()
for ( i in 1:nrow(combs)){
  
#Disparity 95% CIs
  
  ci_x <- apply(boot_disp,2,quantile,probs=c(0.05,0.95))[,combs[i,]][,1]
  ci_y <- apply(boot_disp,2,quantile,probs=c(0.05,0.95))[,combs[i,]][,2]
  
#standard errors (SE)
  n_x <- nrow(boot_disp)
  n_y <- nrow(boot_disp)
  
  SE_x <- (ci_x[2] - ci_x[1]) / (2 * qt(0.975, df = n_x - 1))
  SE_y <- (ci_y[2] - ci_y[1]) / (2 * qt(0.975, df = n_y - 1))
  
# Z-score from Zou method
  Z <- abs((mean(ci_x) - mean(ci_y)) / sqrt(SE_x^2 + SE_y^2))
  Z <- round(Z,3)
  
#p-value
  p_value <- 2 * (1 - pnorm(Z))
  
#Store results in list
  Z.comparisons[[i]] <- setNames(c(Z,p_value), c('Z-score','p-value'))
  
  
}

names(Z.comparisons) <- apply(t(combn(names(disp_groups),2)),1,function(x) paste0(x[1],'_',x[2]) )

Z.comparisons

#Adjusted p-values (Bonferroni)
p.adjust(sapply(Z.comparisons, function(x) x[2]),method = 'bonferroni')


#PCs 1-2 shape variation


open3d()
par3d(windowRect = c(0,0,500,500))
Sys.sleep(1)
mfrow3d(nr = 2, nc = 2, byrow = TRUE, sharedMouse = TRUE)


PCwarps = shape.predictor(GPA$coords[,, specimens_to_use ],
                          x=PCA_subset$PCscores[,1],Intercept = T,
                          min=min(PCA_subset$PCscores[,1]),
                          max=max(PCA_subset$PCscores[,1]))

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




PCwarps = shape.predictor(GPA$coords[,, specimens_to_use ],
                          x=PCA_subset$PCscores[,2],Intercept = T,
                          min=min(PCA_subset$PCscores[,2]),
                          max=max(PCA_subset$PCscores[,2]))

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




#########

######EXTRA STUFF#########
##########################


#Plot by taxonomic groups


par(mfrow=c(2,2))

#'small' specimens

clades.small <- clades
clades.small <- clades.small[groups.temp=='small']


plot(PCA_subset$PCscores[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4),main='small')

points(PCA.small[,c(1,2)], 
       pch=c(3,4,8,15,16,17,18)[as.numeric(as.factor(clades.small))],
       cex=1.35,
       col = adjustcolor('black',alpha.f=0.55),
       bg = adjustcolor('#fcc5c0',alpha.f=0.8))
legend('bottomright', legend = levels(as.factor(clades.small)),
       pch=c(3,4,8,15,16,17,18), bty='n',cex=0.5)

#'intermediate' specimens

clades.interm <- clades
clades.interm <- clades.interm[groups.temp=='interm']


plot(PCA_subset$PCscores[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4),main='intermediate')

points(PCA.interm, 
       pch=c(3,4,8,15,16,17,18)[as.numeric(as.factor(clades.interm))],
       cex=1.35,
       col = adjustcolor('black',alpha.f=0.55),
       bg = adjustcolor('#fcc5c0',alpha.f=0.8))
legend('bottomright', legend = levels(as.factor(clades.interm)),
       pch=c(3,4,8,15,16,17,18), bty='n',cex=0.5)


#'intermediate' specimens

clades.large <- clades
clades.large <- clades.large[groups.temp=='large']


plot(PCA_subset$PCscores[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4),main='large')

points(PCA.large, 
       pch=c(3,4,8,15,16,17,18)[as.numeric(as.factor(clades.large))],
       cex=1.35,
       col = adjustcolor('black',alpha.f=0.55),
       bg = adjustcolor('#fcc5c0',alpha.f=0.8))
legend('bottomright', legend = levels(as.factor(clades.large)),
       pch=c(3,4,8,15,16,17,18), bty='n',cex=0.5)




