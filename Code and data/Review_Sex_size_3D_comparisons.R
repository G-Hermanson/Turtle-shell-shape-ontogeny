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

