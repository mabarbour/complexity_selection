######################################################################################
# Function to compute trait evenness index TED (Fontana et al. 2015)                 
#                                                                                     
# It requires R libraries 'geozoo', 'flexmix' and 'geometry'                                                              #                                                
# input:                                                                             
# - data frame 'traits': every row represents an individual, every column a trait    
# - NA are not allowed                                                                         
# - ADVICE: it is generally meaningful to standardise trait values (mean=0 and sd=1)
######################################################################################



#####################
##    TED index    ##
#####################
# comparison between sample and an assemblage of equidistant points (with the same number of individuals)

library(geozoo)
library(flexmix)
library(geometry)


#############################
# List of possible REFERENCES
#############################

# Open biggest sample
# Change filename to name of file with maximum points
control.survivors <- filter(gall_selection.df, Treatment.focus=="Control", gall_survival==1) %>% 
  select(Plant_Position, sc.Gall_Height_mm, sc.gall_individuals, sc.Density_per_100_shoots)
nrow(control.survivors)

treatment.survivors <- filter(gall_selection.df, Treatment.focus=="Ectoparasitoid exclusion", gall_survival==1) %>% 
  select(Plant_Position, sc.Gall_Height_mm, sc.gall_individuals, sc.Density_per_100_shoots)
nrow(treatment.survivors)

traits.max <- treatment.survivors[,-1] #read.csv("File_with_maximum_points.csv",sep=",",header=T)

# Define maximum number of points (max1) and number of traits under consideration (dim1)
# Alternatively, it is possible to manually define max1 and dim1!!
max1 <- nrow(traits.max)
dim1 <- ncol(traits.max)

ref.matrix<-matrix(ncol=2,nrow=max1)
if (dim1 == 1) {
  i=0.9 } else { i=1.9 }
n <- 0
rows1<-0

while(rows1<max1){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
  
}

k <- i+1
while(i<k){
  i=i+0.1
  n=n+1
  traits.ref <- sphere.solid.grid(p=dim1, n=i)
  rows1<-nrow(traits.ref$points)
  ref.matrix[n,1]<-i
  ref.matrix[n,2]<-rows1
}

ref.matrix<-na.omit(ref.matrix)


##############################
##  TED index calculation   ##
##############################

TED.index <- function(traitdat){
  
  ##########################################################################################
  # Find the best REFERENCE (minimum number of individuals >= individuals in the sample)
  ##########################################################################################
  
  n.sample<-nrow(traitdat)
  
  diff1<-matrix(ncol=2,nrow=length(ref.matrix)/2)
  diff1[,1] <- ref.matrix[,1]
  diff1[,2] <- ref.matrix[,2]-n.sample
  min.diff1<-min(diff1[,2][diff1[,2]>=0])
  select.i<-diff1[diff1[,2]==min.diff1,][1]
  traits.ref <- sphere.solid.grid(p=dim1, n=select.i)
  
  ###################################
  # Transform REFERENCE in data frame
  ###################################
  
  traits.ref <- as.vector(traits.ref$points)
  ind<-length(traits.ref)/dim1
  reference<-matrix(ncol=dim1,nrow=ind)
  for (j in 1:dim1){
    reference[,j] <- traits.ref[((j-1)*ind+1):(j*ind)]
  }
  traits.ref <- as.data.frame(reference)
  
  
  ############################################################################
  # Ev. delete individuals in order to have the same number as in the sample
  ############################################################################
  
  x <- nrow(traits.ref)-nrow(traitdat)
  
  if (x!=0){
    
    # coordinates of the center of gravity of the vertices (Gv)
    baryv<-apply(traits.ref,2,mean)
    
    # euclidian dstances to Gv (dB)
    distbaryv<-rep(0,nrow(traits.ref))
    for (j in 1:nrow(traits.ref))
      distbaryv[j]<-( sum((traits.ref[j,]-baryv)^2) )^0.5
    
    merge1<-data.frame(traits.ref,distbaryv)
    
    #sort by distbaryv (descending)
    sort1 <- merge1[order(-distbaryv),]
    traits.ref<-sort1[-1:-x,-(ncol(sort1))]
    
  }
  
  
  #######################
  # Compare with sample
  #######################
  
  Distance.method <- "euclidean"
  D1 <- dist(traitdat, method=Distance.method)
  density.D <- density(D1)$y
  rm(D1)
  D.ref <- dist(traits.ref, method=Distance.method)
  density.D.ref <- density(D.ref)$y
  rm(D.ref)
  
  results <- KLdiv(cbind(density.D, density.D.ref))
  
  value <- results[1,2]
  
  TED <- 1-log10(value+1)
  TED
  
}

TED.index(traitdat = control.survivors[,-1])
TED.index(traitdat = treatment.survivors[,-1])

uniq.plants <- unique(control.survivors$Plant_Position)
uniq.plants <- names(which(table(control.survivors$Plant_Position)>2)) # ensure sufficient sample size

control.TEDs <- c()
for(i in 1:length(uniq.plants)){
  tmp.traits <- filter(control.survivors, Plant_Position == uniq.plants[i])
  control.TEDs[i] <- TED.index(tmp.traits[,-1])
}

uniq.plants <- unique(treatment.survivors$Plant_Position)
uniq.plants <- names(which(table(treatment.survivors$Plant_Position)>2)) # ensure sufficient sample size

treatment.TEDs <- c()
for(i in 1:length(uniq.plants)){
  tmp.traits <- filter(treatment.survivors, Plant_Position == uniq.plants[i])
  treatment.TEDs[i] <- TED.index(tmp.traits[,-1])
}

mean(control.TEDs)
mean(treatment.TEDs)
