# functions
loglikelihood.gompertz <- function(x, param1, param2){
   # param1 is gamma; param2 is lambda
   if ((param1 <= 0) | (param2 <= 0)){ value <- -Inf 
   }else{
      n <- length(x)
      sum <- sum(exp(param1 * x) - 1)
      value <- n * log(param2) + param1 * sum(x) - (param2/param1)*sum
   }
   value
}

#n * log(b) + a*sum(x) - (b/a)*sum(exp(a*x)-1)


prior.gompertz <- function(param1, param2){
   # param1 is gamma; param2 is lambda
   if ((param1 <= 0) | (param2 <= 0)){ value <- 0 }else{ value <- 1 }
   value
}

posterior.gompertz <- function(x, param1, param2){
   # param1 is gamma; param2 is lambda
   if ((param1 <= 0) | (param2 <= 0)){ value <- NA 
   }else{ 
      l <- loglikelihood.gompertz(x, param1, param2)
      p <- prior.gompertz(param1, param2)
      value <- l + p
   }
   value
}

#########
loglikelihood.gompertzmakeham <- function(x, param1, param2, param3){
   # param1 is gamma; param2 is lambda; param 3 is c
   if ((param1 <= 0) | (param2 <= 0)){ value <- 0 
   }else{
      n <- length(x)
      sum <- sum(exp(param1 * x) - 1)
      sum2 <- sum(log(param3 + param2 * exp(param1 * x)))
      value <- sum2 - (param2/param1) * sum - param3 * sum(x)
   }
   value
}

prior.gompertzmakeham <- function(param1, param2, param3){
   # param1 is gamma; param2 is lambda; param 3 is c
   if ((param1 <= 0) | (param2 <= 0) | (param3 <= 0)){ value <- 0 }else{ value <- 1 }
   value
}

posterior.gompertzmakeham <- function(x, param1, param2, param3){
   # param1 is gamma; param2 is lambda; param 3 is c
   if ((param1 <= 0) | (param2 <= 0)){ value <- NA 
   }else{ 
      l <- loglikelihood.gompertzmakeham(x, param1, param2, param3)
      p <- prior.gompertzmakeham(param1, param2, param3)
      value <- l + p
   }
   value
}
