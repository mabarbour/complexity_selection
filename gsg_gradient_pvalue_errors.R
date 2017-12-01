## REPRODUCIBLE EXAMPLE OF POTENTIAL ERROR 
## IN CALCULATION OF BOOT-STRAPPED P-VALUES
## Matt Barbour, 1 Dec. 2017

## FOLLOWING THE EXAMPLE IN GAM.GRADIENTS HELP FILE ----

# simulated data (stabilizing selection)
z<-rnorm(200,0,2)
W<-rpois(200,exp(1-0.3*z^2))
d<-as.data.frame(list(W=W,z=z))

# characterize the fitness function
library(mgcv)
ff<-gam(W~s(z),family='poisson',data=d)

# derive selection gradients
gradients <- gam.gradients(mod=ff,phenotype="z",se.method='boot.para',standardized=FALSE)

## CURRENT METHOD ----
# Bootstrapped estimates are compared to zero, not the observed estimate.

# B-z: P = 0.934
2*min(sum((gradients$boot[ ,1]>0)+0), sum((gradients$boot[ ,1]<0)+0))/length(gradients$boot[ ,2])
gradients$ests$P.value[1]

# G-z: P = 0.000
2*min(sum((gradients$boot[ ,2]>0)+0), sum((gradients$boot[ ,2]<0)+0))/length(gradients$boot[ ,2])
gradients$ests$P.value[2]

## PROPOSED ALTERNATIVE ----
# Bootstrapped estimates need to be compared to the observed estimate, not zero.

# P = 0.746
2*min(sum(gradients$boot[ ,1]>gradients$ests$estimates[1]), sum(gradients$boot[ ,1]<gradients$ests$estimates[1]))/length(gppr.size_indiv$boot[ ,1])

# P = 0.574
2*min(sum(gradients$boot[ ,2]>gradients$ests$estimates[2]), sum(gradients$boot[ ,2]<gradients$ests$estimates[2]))/length(gppr.size_indiv$boot[ ,2])

