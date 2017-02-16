`srvhaz` <-
function(x,a,b,c=0,s=0,i=1){x<-x*i;if(is.na(c)) c <- 0; if(is.na(s)) s <- 0;c+a*exp(b*x)/(1+s*a*(exp(b*x)-1)/b);}

