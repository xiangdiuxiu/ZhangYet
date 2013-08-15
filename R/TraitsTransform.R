library(MASS)

trait.trf <- function (trait,covt,type){
  if(type==3) trait <- as.ordered (trait)
  data <- data.frame (trait,covt)
  if(type==1){
    if(!is.matrix(covt)){
      res.lm <- lm (trait~1,data=data)
    }else{
      res.lm <- lm (trait~.,data=data)
    }
    index <- which(!is.na(trait))
    output <- rep(NA,length(trait))
    output[index] <- res.lm$residuals
    return (output)
  }else if (type==2){
    if(!is.matrix(covt)){ 
      res.glm <- glm(trait~1,data=data,family="binomial")
    }else{
      res.glm <- glm(trait~.,data=data,family="binomial")  
    }
    if(is.matrix(covt)){
      predictor <- cbind(1,covt)%*%res.glm$coefficients                    
    }else{
      predictor <- res.glm$coefficients          
    }
    output <- trait-exp(predictor)/(1+exp(predictor))
    return(output)
  }else if (type==3){      
    if(!is.matrix(covt)){
      res.plr <- polr(trait~1.,method="probit")
    }else{
      res.plr <- polr(trait~.,data=data,method="probit")        
    }
    index <- which(!is.na(trait))
    fitted.values <- rep(NA,length(trait))
    fitted.values[index] <- res.plr$fitted.values
    n.level <- ncol(fitted.values)
    level <- as.numeric(res.plr$lev)
    output <- apply(data.frame(fitted.values,as.numeric(trait)-1),1,
                    function(x){
                      obs <- x[n.level+1]
                      if(obs==level[1]){
                        gamma1 <- x[1]
                        gamma2 <- 0
                      }else{
                        gamma1 <- sum(x[1:which(level==obs)])
                        gamma2 <- sum(x[1:(which(level==obs)-1)])                              
                      }
                      return(1-gamma1-gamma2)
                    }
    )
    return(output)
  }
}

file.trf <- function(phenofile=NULL,famfile=NULL,phenotypes=NULL, pheno_type=NULL,covar_name=NULL) {  
    dat <- read.table(phenofile,header=T)
    fam <- read.table(famfile,header=F)
    com.ind <- intersect(dat[,2],fam[,2])
    index <- dat[,2]%in%com.ind
    dat <- dat[index,]
    t_index <- names(dat)%in%phenotypes
    covt_index <- names(dat)%in%covar_name   
    
    trait <- dat[,t_index]
    covt <-  dat[,covt_index]
    trans.y <- NULL
    n.trait <- length(pheno_type)
    
    trait <- as.data.frame(trait)
    if(!is.null(covt)){covt <- as.matrix(covt)} 
    
    for (k in 1:n.trait){
        tmp <- trait.trf(trait[,k],covt,pheno_type[k]) # transform the trait #
        trans.y <- cbind(trans.y,tmp)
    }
        
    dat[,t_index] <- trans.y
    output.name <- paste(phenofile,"_trf",sep="")
    write.table(dat, file=output.name,row.names=F,quote=F)
}

