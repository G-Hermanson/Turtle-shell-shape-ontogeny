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
  
  
}



# Figure 2 #
#Ontogenetic shell shape curves of selected taxa

par(mfrow=c(3,2))

#A-D

taxa_plot <- c('Kinosternon_creaseri','Trachemys_ornata','Manouria_emys','Heosemys_grandis')

for ( i in taxa_plot){
  
  mydata <- as.data.frame(CAC.list[[i]])
  mydata$V2 <- 10^mydata$V2
  lm.tmp <- drm(V1~V2,data=mydata,fct=G.4())
  
  
  specimens <- sapply(strsplit(rownames(mydata),'_'), function(x) x[4] )
  sex <- subset(landmark_data,Name==i)
  sex <- sex$Sex[ sex$Museum.Number %in% specimens ]
  females <- which(sex=='F')
  females <- setNames(rep('F',length(females)),females)
  males <- which(sex=='M')
  males <- setNames(rep('M',length(males)),males)
  unknown <- grep('?',sex,fixed = T)
  unknown <- setNames(rep('uncertain',length(unknown)),unknown)
  juvenile <- which(sex=='J')
  juvenile <- setNames(rep('uncertain',length(juvenile)),juvenile)
  
  myorder <- c(females,males,unknown,juvenile)[as.character(1:length(sex))]
  
  mypch <- setNames(c(16,17,18),c('F','M','uncertain'))
  mypch <- mypch[names(mypch) %in% unique(myorder)]
  
  plot(x=10^CAC.list[[i]][,2],y=CAC.list[[i]][,1], xlab='SCL',ylab='CAC',
       pch=mypch [as.numeric(as.factor(myorder))] , 
       col=c('#fcc5c0','#dd3497','#7a0177') [as.numeric(groups.list[[i]])] ,
       main=i, bty='l')
  
  
  new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l')
  
  
  legend('topleft',cex=0.7,legend=levels(groups.list[[i]]),
         fill=c('#fcc5c0','#dd3497','#7a0177'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 
  
  
}


#E-F

taxa_plot <- c('Podocnemis_vogli','Trachemys_scripta')

for ( i in taxa_plot){
  
  
  mydata <- as.data.frame(CAC.list[[i]])
  mydata$V2 <- 10^mydata$V2
  lm.tmp <- drm(V1~V2,data=mydata,fct=G.4())
  
  
  specimens <- sapply(strsplit(rownames(mydata),'_'), function(x) x[4] )
  sex <- subset(landmark_data,Name==i)
  sex <- sex$Sex[ sex$Museum.Number %in% specimens ]
  females <- which(sex=='F')
  females <- setNames(rep('F',length(females)),females)
  males <- which(sex=='M')
  males <- setNames(rep('M',length(males)),males)
  unknown <- grep('?',sex,fixed = T)
  unknown <- setNames(rep('uncertain',length(unknown)),unknown)
  juvenile <- which(sex=='J')
  juvenile <- setNames(rep('uncertain',length(juvenile)),juvenile)
  
  myorder <- c(females,males,unknown,juvenile)[as.character(1:length(sex))]
  
  mypch <- setNames(c(16,17,18),c('F','M','uncertain'))
  mypch <- mypch[names(mypch) %in% unique(myorder)]
  
  plot(x=10^CAC.list[[i]][,2],y=CAC.list[[i]][,1], xlab='SCL',ylab='CAC',
       pch=mypch [as.numeric(as.factor(myorder))] , 
       col=c('#fcc5c0','#dd3497','#7a0177') [as.numeric(groups.list[[i]])] ,
       main=i, bty='l')
  
  legend('topleft',cex=0.7,legend=levels(groups.list[[i]]),
         fill=c('#fcc5c0','#dd3497','#7a0177'),bty='n',border = NA ) 
  
  mypch <- levels(as.factor(myorder))
  mypch <- c(16,17,18) [c('F','M','uncertain')%in% mypch]
  
  legend('bottomright',cex=0.7,legend=levels(as.factor(myorder)),
         pch=mypch,
         bty='n' ) 
  
  #female curves
  female <- mydata[which(subset(landmark_data,Name==i)[,'Sex'] !='M'),]
  female <- female[complete.cases(female),]
  
  lm.tmp <- drm(V1~V2,data=female,fct=G.4())
  new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l',lty=3,lwd=1.25)
  
  
  #male curves
  male <- mydata[which(subset(landmark_data,Name==i)[,'Sex'] =='M'),]
  male <- male[complete.cases(male),]
  
  lm.tmp <- drm(V1~V2,data=male,fct=G.4())
  new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
  points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l',lty=2)
  
  
}


dev.off()


# Figure 3 #
#Disparity (sum of ranges) per ontogenetic stage

#This figure can be created using "Script 2", from lines 8-128, and was later modified in Adobe Illustrator



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


# Supplementary Figures 2-9 #
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
#These refer to the bivariate plots between CAC and SCL values

#Plots were made in "Script 1", lines 217-614, and later modified in Adobe Illustrator


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

plot(mydata[,2:1], xlab='SCL',ylab='CAC',bty='l',cex=1.5,
     pch=18, col=cols[as.numeric(groups.list[[take]])])
lm.tmp <- drm(V1~V2,data=mydata,fct=G.4())
new <- seq(min(mydata$V2), max(mydata$V2), length.out=100)
points( predict(lm.tmp, newdata = data.frame(V2=new)) ~ new, type='l',lty=1)

growth_tmp <- growth_data[ grep(take,growth_data$Taxon),c('Size_class','Growth_rate')]
growth_tmp <- growth_tmp[complete.cases(growth_tmp),]

par(new=T)
plot(growth_tmp,
     type='o',pch=16,axes=F,xlab='',ylab='',xlim=range(mydata$V2),ylim=c(0,9),lwd=1.5,
     col=adjustcolor('grey',0.5),cex=sqrt(growth_data$N)/2)
mtext("Growth rate (cm/year)",side=4,line=2) 
axis(4, ylim=c(0,9),las=1,labels = (0:3)*3, at=(0:3)*3)

for ( i in 1:nrow(growth_data)){
  
  segments(x0=growth_tmp[i,1],x1=growth_tmp[i,1],
           y0=growth_tmp[i,2]-growth_data[i,'SD'], y1=growth_tmp[i,2]+growth_data[i,'SD'],
           col=adjustcolor('grey',0.3),lwd=2)
  
}


dev.off()
