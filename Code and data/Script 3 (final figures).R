#############################################################################
### Script 3 - Create figures in the main text and supplementary material ###
#############################################################################

###############
## MAIN TEXT ##
###############


# Figure 1 #
#Shell shape changes in selected taxa


taxa_plot <- c('Podocnemis_expansa','Eretmochelys_imbricata',
               'Chelonoidis_nigra','Kinosternon_creaseri')

###Shape deformations (loop through taxa to plot)

open3d()
par3d(windowRect = c(0,0,500,500))
Sys.sleep(1)
mfrow3d(nr = 4, nc = 3, byrow = TRUE, sharedMouse = TRUE)

for ( i in taxa_plot){
  
  take=i
  
  quarts <- quantile(CAC.list[[take]][,1],probs=c(0.05,0.5,0.95))
  
  PCwarps = shape.predictor(GPA$coords[,, rownames(CAC.list[[take]]) ],
                            x=CAC.list[[take]] [,1],Intercept = T,
                            min=min(quarts),max=max(quarts),mean=quarts[2])
  
  
  
  
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
  
  
}



# Figure 2 #
#Ontogenetic shell shape curves of selected taxa
#CAC

par(mfrow=c(3,2))

#A-D

#taxa_plot <- c('Kinosternon_creaseri','Trachemys_ornata','Manouria_emys','Heosemys_grandis')


taxa_plot <- c('Kinosternon_creaseri','Trachemys_ornata','Manouria_emys','Podocnemis_expansa')

par(mfrow=c(2,2))

for ( i in taxa_plot){
  
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
       xlim=c(min(mydata$V2), max.sizes[i]),
       col=c('#ccebc5' , '#4eb3d3' , '#084081') [as.numeric(groups.list[[i]])] ,
       main=i, bty='l')
  
  
  new <- seq(min(mydata$V2), max.sizes[i], length.out=100)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l')
  
  
  legend('topleft',cex=0.7,legend=levels(groups.list[[i]]),
         fill=c('#ccebc5' , '#4eb3d3' , '#084081'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 
  
#Add thresholds
  
real.thresh.CAC <- thresholds.CAC[i] * (max.sizes[i])
  
  rect(xleft = real.thresh.CAC, xright = max(preds.CAC[[i]][,'CS'])*1.2,
       ybottom =min(preds.CAC[[i]][,1])*1.2,ytop =max(preds.CAC[[i]][,1])*1.2,
       border=NA, col=adjustcolor('grey',alpha.f = 0.4))
  
  abline(v=real.thresh.CAC, col='red',lty=1) #CAC threshold
  
  bounds <- confint(lm.tmp,level = 0.9)[3,]
  
  segments(x0=max(preds.list[[i]][,3])*0.4,x1=max(preds.list[[i]][,3])*1.2,
           y0=bounds[1],y1=bounds[1],lty=2)
  
  
  
}



# Figure 3 #
#Ontogenetic shell shape curves of selected taxa
#Multivariate shape

#A-D

taxa_plot <- c('Chelus_fimbriata','Trachemys_scripta','Manouria_emys','Orlitia_borneensis')

par(mfrow=c(2,2))

for ( i in taxa_plot){
  
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
  
  plot(x=10^PredLine.list[[i]][,2],y=PredLine.list[[i]][,1], xlab='SCL',ylab='Multivariate shape',
       pch=mypch [as.numeric(as.factor(myorder))] , 
       xlim=c(min(mydata$V2), max.sizes[i]),
       col=c('#ccebc5' , '#4eb3d3' , '#084081') [as.numeric(groups.list[[i]])] ,
       main=i, bty='l')
  
  
  new <- seq(min(mydata$V2), max.sizes[i], length.out=100)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l')
  
  
  legend('topleft',cex=0.7,legend=levels(groups.pred[[i]]),
         fill=c('#ccebc5' , '#4eb3d3' , '#084081'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 
  
  #Add thresholds
  
  real.thresh <- thresholds[i] * (max.sizes[i])
  
  rect(xleft = real.thresh, xright = max(preds.list[[i]][,'CS'])*1.2,
       ybottom =min(preds.list[[i]][,1])*1.2,ytop =max(preds.list[[i]][,1])*1.2,
       border=NA, col=adjustcolor('grey',alpha.f = 0.4))
  
  abline(v=real.thresh, col='red',lty=1) #CAC threshold
  
  bounds <- confint(lm.tmp,level = 0.9)[3,]
  
  segments(x0=max(preds.list[[i]][,3])*0.4,x1=max(preds.list[[i]][,3])*1.2,
           y0=bounds[1],y1=bounds[1],lty=2)
  
  
  
}


dev.off()


# Figure 4 #
#Effects of size and sex on shell shape variation

#This figure was later modified in Adobe Illustrator


par(mfcol=c(2,2))

#Differece between Z-scores of 'size' and 'sex'

#Panel A
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

#Panel B
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


#Panels C-E
#Leucocephalon comparison

#Female vs male

fem <- 'Leucocephalon_yuwonoi_FMNH_261568'
mal <- 'Leucocephalon_yuwonoi_FLMNH_109835'

#Simple plot

open3d()
par3d(windowRect = c(0,0,500,500))
Sys.sleep(1)
mfrow3d(nr = 3, nc = 3, byrow = TRUE, sharedMouse = TRUE)


choose=fem
plot3d(GPA$coords[,,choose] , size=1.2,col='black' ,add=F,type='s',
       box=F,axes=F,lit=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,choose] [lines_list[[i]],] , lwd=1, col='black') }
#title3d(main='female')

choose=mal
plot3d(GPA$coords[,,choose] , size=1.2,col='black' ,add=F,type='s',
       box=F,axes=F,lit=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,choose] [lines_list[[i]],] , lwd=1, col='black') }
#title3d(main='male')



#### Custom function for calculation of Euclidean distance ####

Edist <- function ( x , Y ) { ( sum( ( x - Y ) ^ 2 ) ) ^ 0.5 }

#plot
#open3d()


# define the colours you want
point.distance.scale <- colorRamp(c("lightgrey" ,'blue')) 

point.distances <- c()
for ( i in 1:nrow(GPA$coords[,,fem]) ){
  
  # here you calculate the Euclidean distance between the points in each of your shape matrices ('shape1' and 'shape2')
  point.distances[i] <- Edist(GPA$coords[i,,fem] , GPA$coords[i,,mal]) 
}	

#point.distances

# normalise the distances so they range from 0 to 1
point.distances.norm <- (point.distances - min(point.distances)) / max ( point.distances - min(point.distances))

# and then you're able to create a colorRamp that goes from 0 (grey) to 1 (red)
point.colours <- point.distance.scale(point.distances.norm)


choose = fem

#min
plot3d(GPA$coords[,,choose] , size=1.1,col='lightgrey',lit=F ,add=F,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,choose] [lines_list[[i]],] , lwd=1, col='lightgrey') }

#title3d(main='female -> male')

plot3d(GPA$coords[,,mal], size = 1.4, lit=F ,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso',add=T,
       col=rgb(point.colours,maxColorValue = 255) )
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,mal] [lines_list[[i]],] , 
          lwd=1, col=rgb(point.colours,maxColorValue = 255) ) }

################################



#Smallest male vs largest male

#Female vs male

small <- 'Leucocephalon_yuwonoi_AMNH_145108'
large <- 'Leucocephalon_yuwonoi_FLMNH_111310'

#Simple plot

#open3d()
#par3d(windowRect = c(0,0,500,500))
#Sys.sleep(1)
#mfrow3d(nr = 1, nc = 3, byrow = TRUE, sharedMouse = TRUE)


choose=small
plot3d(GPA$coords[,,choose] , size=1.2,col='black' ,add=F,type='s',
       box=F,axes=F,lit=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,choose] [lines_list[[i]],] , lwd=1, col='black') }
#title3d(main='small')

choose=large
plot3d(GPA$coords[,,choose] , size=1.2,col='black' ,add=F,type='s',
       box=F,axes=F,lit=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,choose] [lines_list[[i]],] , lwd=1, col='black') }
#title3d(main='large')



#### Custom function for calculation of Euclidean distance ####

Edist <- function ( x , Y ) { ( sum( ( x - Y ) ^ 2 ) ) ^ 0.5 }

#plot
#open3d()


# define the colours you want
point.distance.scale <- colorRamp(c("lightgrey" ,'blue')) 

point.distances <- c()
for ( i in 1:nrow(GPA$coords[,,small]) ){
  
  # here you calculate the Euclidean distance between the points in each of your shape matrices ('shape1' and 'shape2')
  point.distances[i] <- Edist(GPA$coords[i,,small] , GPA$coords[i,,large]) 
}	

#point.distances

# normalise the distances so they range from 0 to 1
point.distances.norm <- (point.distances - min(point.distances)) / max ( point.distances - min(point.distances))

# and then you're able to create a colorRamp that goes from 0 (grey) to 1 (red)
point.colours <- point.distance.scale(point.distances.norm)


choose = small

#min
plot3d(GPA$coords[,,choose] , size=1.1,col='lightgrey',lit=F ,add=F,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,choose] [lines_list[[i]],] , lwd=1, col='lightgrey') }

#title3d(main='small -> large')

plot3d(GPA$coords[,,large], size = 1.4, lit=F ,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso',add=T,
       col=rgb(point.colours,maxColorValue = 255) )
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,large] [lines_list[[i]],] , 
          lwd=1, col=rgb(point.colours,maxColorValue = 255) ) }


####

#Interspecific comparison

#Largest Leucocephalon male vs largest Hieremys annandalii male (closest relative with large adults of same sex)

#large_Ha <- 'Heosemys_annandalii_USNM_63039'
large_Ha <- 'Notochelys_platynota_FMNH_151017'
large <- 'Leucocephalon_yuwonoi_FLMNH_111310'

#Simple plot

#open3d()
#par3d(windowRect = c(0,0,500,500))
#Sys.sleep(1)
#mfrow3d(nr = 1, nc = 3, byrow = TRUE, sharedMouse = TRUE)


choose=large_Ha
plot3d(GPA$coords[,,choose] , size=1.2,col='black' ,add=F,type='s',
       box=F,axes=F,lit=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,choose] [lines_list[[i]],] , lwd=1, col='black') }
#title3d(main='Hieremys')

choose=large
plot3d(GPA$coords[,,choose] , size=1.2,col='black' ,add=F,type='s',
       box=F,axes=F,lit=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,choose] [lines_list[[i]],] , lwd=1, col='black') }
#title3d(main='Leucocephalon')


#### Custom function for calculation of Euclidean distance ####

Edist <- function ( x , Y ) { ( sum( ( x - Y ) ^ 2 ) ) ^ 0.5 }

#plot
#open3d()


# define the colours you want
point.distance.scale <- colorRamp(c("lightgrey" ,'blue')) 

point.distances <- c()
for ( i in 1:nrow(GPA$coords[,,large_Ha]) ){
  
  # here you calculate the Euclidean distance between the points in each of your shape matrices ('shape1' and 'shape2')
  point.distances[i] <- Edist(GPA$coords[i,,large_Ha] , GPA$coords[i,,large]) 
}	

#point.distances

# normalise the distances so they range from 0 to 1
point.distances.norm <- (point.distances - min(point.distances)) / max ( point.distances - min(point.distances))

# and then you're able to create a colorRamp that goes from 0 (grey) to 1 (red)
point.colours <- point.distance.scale(point.distances.norm)


choose = large_Ha

#min
plot3d(GPA$coords[,,choose] , size=1.1,col='lightgrey',lit=F ,add=F,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,choose] [lines_list[[i]],] , lwd=1, col='lightgrey') }

#title3d(main='Hieremys -> Leucocephalon')

plot3d(GPA$coords[,,large], size = 1.4, lit=F ,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso',add=T,
       col=rgb(point.colours,maxColorValue = 255) )
for ( i in 1:length(lines_list)){
  lines3d(GPA$coords[,,large] [lines_list[[i]],] , 
          lwd=1, col=rgb(point.colours,maxColorValue = 255) ) }



# Figure 5 #
#Disparity#

#This figure was created using "Script 2" and later modified in Adobe Illustrator



########################
## SUPPLEMENTARY TEXT ##
########################


# Supplementary Figure 1 #
#Landmark configuration 

plot3d(GPA$consensus , size=1.1,col='lightgrey',lit=F ,add=F,type='s',
       box=F,axes=F,xlab='',zlab='',ylab='',aspect = 'iso')
for ( i in 1:length(lines_list)){
  lines3d(GPA$consensus [lines_list[[i]],] , lwd=1, col='lightgrey') }

text3d(GPA$consensus[1:40,], texts = 1:40, pos=3,cex=1)
#text3d(GPA$consensus[-c(1:40),], texts = 41:nrow(GPA$consensus), pos=3,cex=1)

rgl.snapshot('filename.png',fmt = 'png')


# Supplementary Figures #
#These refer to species-specific landmark configurations showing CAC changes from smaller to larger individuals

#Change the object 'take' with the name of a specific taxon

take='Pelusios_sinuatus'

quarts <- quantile(CAC.list[[take]][,1],probs=c(0.05,0.95))

PCwarps = shape.predictor(GPA$coords[,, rownames(CAC.list[[take]]) ],
                          x=CAC.list[[take]] [,1],Intercept = T,
                          min=min(quarts),max=max(quarts))

open3d()
par3d(windowRect = c(0,0,500,500))
Sys.sleep(1)
mfrow3d(nr = 2, nc = 2, byrow = TRUE, sharedMouse = TRUE)


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



#### Custom function for calculation of Euclidean distance ####

Edist <- function ( x , Y ) { ( sum( ( x - Y ) ^ 2 ) ) ^ 0.5 }

#plot
#open3d()


# define the colours you want
point.distance.scale <- colorRamp(c("lightgrey" ,'red')) 

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


# Supplementary Figures 10-15 #
#These refer to the bivariate plots between CAC/Multivariate shape and SCL values



# Supplementary Figure 16 #
#Growth rate comparison Chelonia mydas
#Data from Bjorndal & Bolten (1988) and Kordikova (2002)

pdf('CAC vs growth rate Chelonia mydas.pdf',width = 6,height = 6,useDingbats = F)

growth_data <- read.csv('sea_turtle_growth_Kordikova.csv',sep=',',header = T,row.names = NULL)
growth_data$Size_class <- growth_data$Size_class*10

#Chelonia mydas
take='Chelonia_mydas'

mydata <- as.data.frame(CAC.list[[take]])
mydata$V2 <- 10^mydata$V2

cols = c('#ccebc5' , '#4eb3d3' , '#084081')

plot(mydata[,2:1], xlab='SCL',ylab='CAC',bty='l',cex=1.5,
     pch=18, col=cols[as.numeric(groups.list[[take]])],
     xlim=c(min(mydata$V2),max.sizes[take]))
lm.tmp <- drm(V1~V2,data=mydata,fct=G.4())
new <- seq(min(mydata$V2), max.sizes[take], length.out=100)
points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l',lty=1)

abline(v= threshold[take,'CAC'] * max.sizes[take], col='red')

growth_tmp <- growth_data[ grep(take,growth_data$Taxon),c('Size_class','Growth_rate')]
growth_tmp <- growth_tmp[complete.cases(growth_tmp),]

par(new=T)
plot(growth_tmp,
     type='o',pch=16,axes=F,xlab='',ylab='',
     xlim=c(min(mydata$V2),max.sizes[take]),
     ylim=c(0,9),lwd=1.5,
     col=adjustcolor('grey',0.5),cex=sqrt(growth_data$N)/2)
mtext("Growth rate (cm/year)",side=4,line=2) 
axis(4, ylim=c(0,9),las=1,labels = (0:3)*3, at=(0:3)*3)

for ( i in 1:nrow(growth_data)){
  
  segments(x0=growth_tmp[i,1],x1=growth_tmp[i,1],
           y0=growth_tmp[i,2]-growth_data[i,'SD'], y1=growth_tmp[i,2]+growth_data[i,'SD'],
           col=adjustcolor('grey',0.3),lwd=2)
  
}


dev.off()

dev.off()
