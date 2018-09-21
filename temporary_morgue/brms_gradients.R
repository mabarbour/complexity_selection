
predict(control.major.brm, summary = T)[ ,"Estimate"]
control.major.gam$gam$model

## TEST FUNCTION FOR CALCULATING GRADIENTS WITH BRMS

Wbar <- function(x, mod, phenotype,covariates) {
  new.d <- as.data.frame(mod$data[,c(phenotype,covariates)])
  names(new.d)<-c(phenotype,covariates)
  new.d2<-new.d
  for (i in 1:length(x)) {
    new.d[,as.character(phenotype[i])] <-
      new.d[,as.character(phenotype[i])]+x[i]
  }
  
  p <- predict(
    object=mod,
    newdata=new.d,
    #newdata.guaranteed=TRUE,
    type="response"
  )[ ,"Estimate"]
  return(mean(p))
}

gradients <- function(m, phenotype, covariates) {
  nTraits = length(phenotype)
  first.derivatives <- grad(func=Wbar, x=rep(0, nTraits), 
                            mod=m, phenotype=phenotype, covariates=covariates)
  second.derivatives <- hessian(func=Wbar, x=rep(0, nTraits), 
                                mod=m, phenotype=phenotype, covariates=covariates)
  denom <- Wbar(x=rep(0,nTraits), mod=m, phenotype=phenotype,
                covariates=covariates)
  
  beta <- first.derivatives  / denom
  gamma <- second.derivatives / denom
  
  #if(standardized){
  #  sds<-apply(as.matrix(mod$model[,phenotype]),2,sd)
  #  beta<-beta*sds
  #  gamma<-gamma*outer(sds,sds)
  #}
  return( list(
    beta  = beta,
    gamma = gamma
  ))
}
library(numDeriv)
gradients(m = control.major.brm, phenotype = "term1", covariates = NULL)

# well the first round didn't work...