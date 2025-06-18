############################################################################
### Script 2 - Disparity analysis between pre-defined ontogenetic stages ###
############################################################################


#Shape analysis

specimens_to_use <- rownames(do.call(rbind,CAC.list))

PCA_subset <- PCA$x[specimens_to_use,]

groups.temp <- unname(groups.list)
groups.temp <- unlist(groups.temp)
groups.temp <- groups.temp[rownames(PCA_subset)]

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
                5,5,5),4,byrow = T), widths = rep(1/3,3), heights = c(0.3,0.3,0.6,0.6))

#layout.show(5)

#A 'small' specimens
PCA.small <- names(groups.temp)[groups.temp=='small']

PCA.small <- PCA_subset[PCA.small,]

plot(PCA_subset[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4), xaxt='n',yaxt='n',
     xlab='PC1',ylab='PC2')
points(PCA.small, pch=21,cex=1.25,
       col = '#ccebc5',
       bg = adjustcolor('#ccebc5',alpha.f=0.85))


#B 'intermediate' specimens
PCA.interm <- names(groups.temp)[groups.temp=='interm']

PCA.interm <- PCA_subset[PCA.interm,]

plot(PCA_subset[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4),xaxt='n',yaxt='n',
     xlab='PC1',ylab='')
points(PCA.interm, pch=21,cex=1.25,
       col = '#4eb3d3',
       bg = adjustcolor('#4eb3d3',alpha.f=0.85))


#C 'large' specimens
PCA.large <- names(groups.temp)[groups.temp=='large']

PCA.large <- PCA_subset[PCA.large,]

plot(PCA_subset[,c(1,2)], pch=16, col=adjustcolor('grey',alpha.f = 0.4),xaxt='n',yaxt='n',
     xlab='PC1',ylab='')
points(PCA.large, pch=21,cex=1.25,
       col = '#084081',
       bg = adjustcolor('#084081',alpha.f=0.85))



#Disparity
disp_groups <- list('small'=names(groups.temp)[groups.temp=='small'],
                    'interm'=names(groups.temp)[groups.temp=='interm'],
                    'large'=names(groups.temp)[groups.temp=='large'])

#Sum of ranges
disp_data <- custom.subsets(PCA_subset,
                            group = disp_groups)
set.seed(123)
disp_data <- boot.matrix(disp_data,bootstraps = 1000, rarefaction = 'min')
disp_data_range <- dispRity(disp_data,metric = c(sum,ranges))

boot_disp <- t(do.call(rbind,lapply(disp_data_range$disparity, function(x) x[[2]])))

#Mean Procrustes distance
reps=1000

boot_disp_pw <- matrix(,nrow=reps,ncol = 3)
  colnames(boot_disp_pw) <- levels(groups.temp)

for ( i in 1:reps){
  
  minsize <- min(sapply(disp_groups,length))-1
  
  samples <- lapply( disp_groups, function(x) sample(x,minsize)  )
  
  boot_disp_pw[i,'small'] <- c(geomorph::morphol.disparity(GPA$coords[,,samples$small] ~ 1,print.progress = F) )
  boot_disp_pw[i,'interm'] <- c(geomorph::morphol.disparity(GPA$coords[,,samples$interm] ~ 1,print.progress = F) )
  boot_disp_pw[i,'large'] <- c(geomorph::morphol.disparity(GPA$coords[,,samples$large] ~ 1,print.progress = F) )
  
  
}

#Observed
c(geomorph::morphol.disparity(GPA$coords[,,disp_groups$small] ~ 1 ) )
c(geomorph::morphol.disparity(GPA$coords[,,disp_groups$interm] ~ 1 ) )
c(geomorph::morphol.disparity(GPA$coords[,,disp_groups$large] ~ 1 ) )

apply(boot_disp_pw,2,quantile,probs=c(0.025,0.975))


#plot (sum of ranges disparity)
plot('n',xlim=c(0.5,3.5),ylim=c(5.2,6.6),ylab='Sum of ranges',
     xaxt='n',xlab=NA)
axis(1,at=1:3,labels=c('small','intermediate','large'))

vioplot(boot_disp,
        xlab=NA,ylab=NA,add=T,col=adjustcolor('grey',alpha.f = 0.2),
        rectCol=adjustcolor('grey20',alpha.f = 0.3),
        colMed=adjustcolor('grey20',alpha.f = 0.3),lineCol=adjustcolor('grey20',alpha.f = 0.2),
        border=F)

for ( i in 1:ncol(boot_disp) ){
  
  #cols <- c('#fcc5c0','#dd3497','#7a0177')
  #cols <- c('#d0d1e6' , '#3690c0' , '#023858')
  cols <- c('#ccebc5' , '#4eb3d3' , '#084081')
  
  stripchart(boot_disp[,i], col = adjustcolor(cols[i],alpha.f=0.25),#col='grey70',
             method = 'jitter',add=T,at=i,lwd=0.5,
             pch=16,vertical = T,cex=1.2)
  
}


vioplot(boot_disp,
        xlab=NA,ylab=NA,add=T,col=adjustcolor('grey',alpha.f = 0.1),
        rectCol=adjustcolor('grey20',alpha.f = 0.8),
        colMed=adjustcolor('grey20',alpha.f = 0.8),lineCol=adjustcolor('grey20',alpha.f = 0.6),
        border=F)


#plot (Procrustes variance disparity)
plot('n',xlim=c(0.5,3.5),ylim=c(0.0155,0.0205),ylab='Procrustes variance',
     xaxt='n',xlab=NA)
axis(1,at=1:3,labels=c('small','intermediate','large'))

vioplot(boot_disp_pw,
        xlab=NA,ylab=NA,add=T,col=adjustcolor('grey',alpha.f = 0.2),
        rectCol=adjustcolor('grey20',alpha.f = 0.3),
        colMed=adjustcolor('grey20',alpha.f = 0.3),lineCol=adjustcolor('grey20',alpha.f = 0.2),
        border=F)

for ( i in 1:ncol(boot_disp_pw) ){
  
  #cols <- c('#fcc5c0','#dd3497','#7a0177')
  #cols <- c('#d0d1e6' , '#3690c0' , '#023858')
  cols <- c('#ccebc5' , '#4eb3d3' , '#084081')
  
  stripchart(boot_disp_pw[,i], col = adjustcolor(cols[i],alpha.f=0.25),#col='grey70',
             method = 'jitter',add=T,at=i,lwd=0.5,
             pch=16,vertical = T,cex=1.2)
  
}


vioplot(boot_disp_pw,
        xlab=NA,ylab=NA,add=T,col=adjustcolor('grey',alpha.f = 0.1),
        rectCol=adjustcolor('grey20',alpha.f = 0.8),
        colMed=adjustcolor('grey20',alpha.f = 0.8),lineCol=adjustcolor('grey20',alpha.f = 0.6),
        border=F)


dev.off()


#Test whether disparity 95% CIs are different

#Sum of ranges

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
  
# Z-score from Zou's method (Zou 2007)
  Z <- abs((mean(ci_x) - mean(ci_y)) / sqrt(SE_x^2 + SE_y^2))
  Z <- round(Z,3)
  
#p-value
  p_value <- 2 * (1 - pnorm(Z))
  
#Store results in a list
  Z.comparisons[[i]] <- setNames(c(Z,p_value), c('Z-score','p-value'))
  
  
}

names(Z.comparisons) <- apply(t(combn(names(disp_groups),2)),1,function(x) paste0(x[1],'_',x[2]) )

Z.comparisons

#Adjusted p-values (Bonferroni)
p.adjust(sapply(Z.comparisons, function(x) x[2]),method = 'bonferroni')


#Same for Procrustes variances
combs <- t(combn(1:ncol(boot_disp_pw),2))

Z.comparisons <- list()
for ( i in 1:nrow(combs)){
  
  #Disparity 95% CIs
  
  ci_x <- apply(boot_disp_pw,2,quantile,probs=c(0.05,0.95))[,combs[i,]][,1]
  ci_y <- apply(boot_disp_pw,2,quantile,probs=c(0.05,0.95))[,combs[i,]][,2]
  
  #standard errors (SE)
  n_x <- nrow(boot_disp_pw)
  n_y <- nrow(boot_disp_pw)
  
  SE_x <- (ci_x[2] - ci_x[1]) / (2 * qt(0.975, df = n_x - 1))
  SE_y <- (ci_y[2] - ci_y[1]) / (2 * qt(0.975, df = n_y - 1))
  
  # Z-score from Zou's method (Zou 2007)
  Z <- abs((mean(ci_x) - mean(ci_y)) / sqrt(SE_x^2 + SE_y^2))
  Z <- round(Z,3)
  
  #p-value
  p_value <- 2 * (1 - pnorm(Z))
  
  #Store results in a list
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
                          x=PCA_subset[,1],Intercept = F,
                          min=min(PCA_subset[,1]),
                          max=max(PCA_subset[,1]))

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
                          x=PCA_subset[,2],Intercept = F,
                          min=min(PCA_subset[,2]),
                          max=max(PCA_subset[,2]))

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
                
