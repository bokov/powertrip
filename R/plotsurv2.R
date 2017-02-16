plotsurv2 <- function(sfit, table=T, returns=F, xlabs='\nTime', ylabs='Survival Probability\n',
        ystratalabs = NULL,ystrataname=NULL,timeby=100,main='Kaplan-Meier Plot\n',pvaltxt='',
                      censorshape=4,censorsize=3,xlim,lwd=1,gridtheme=element_blank(),breakby=200,...){

    require(ggplot2)
    require(survival)
    require(gridExtra)
    require(plyr)


#     Create Kaplan Meier plot
    if(is.null(ystratalabs)) ystratalabs <- as.character(levels(summary(sfit)$strata))
    (m<-max(nchar(ystratalabs)))
    if(is.null(ystrataname)) ystrataname <- 'Strata'
    times  <- seq(0, max(sfit$time), by=timeby)
    df  <- data.frame(
            time    =   sfit$time,
            n.risk  =   sfit$n.risk,
            n.event =   sfit$n.event,
            n.censor=   sfit$n.censor,
            surv    =   sfit$surv,
            strata  =   summary(sfit, censored=T)$strata,
            upper   =   sfit$upper,
            lower   =   sfit$lower
            )
    if(missing(xlim)) xlim = c(0,max(sfit$time));
    xbreaks = seq(xlim[1],xlim[2],by=breakby);
    #browser();
    levels(df$strata) <- ystratalabs;
    zeros <- data.frame(time=0, surv=1, strata = ystratalabs,
            upper=1, lower=1)
    df  <- rbind.fill(zeros, df)
    d <- length(levels(df$strata))
    p  <- ggplot(df, aes(time, surv, group=strata))+
        geom_step(aes(color=strata), size=lwd)+
        #theme_bw()+opts(axis.title.x=theme_text(vjust=0.5))+
        theme_bw()+theme(axis.title.x=element_text(vjust=0.5))+          
        scale_x_continuous(xlabs, breaks=xbreaks, limits = xlim)+
        scale_y_continuous(ylabs, limits = c(0,1))+
        #opts(panel.grid.minor = theme_blank())+
        theme(panel.grid.minor = element_blank())+
        #opts(panel.grid.major = gridtheme)+
        #opts(legend.position=c(ifelse(m<10,.28,.35),ifelse(d < 4,.25,.35)))+
        #opts(legend.key=theme_rect(colour=NA))+
        theme(panel.grid.major = gridtheme)+
        theme(legend.position=c(ifelse(m<10,.28,.35),ifelse(d < 4,.25,.35)))+
        theme(legend.key=element_rect(colour=NA))+
        labs(linetype=ystrataname,title=main)+
        #opts(plot.margin = unit(c(0, 1, .5, ifelse(m<10,1.5,2.5)),'lines'))+
        theme(plot.margin = unit(c(0, 1, .5, ifelse(m<10,1.5,2.5)),'lines'))
        #opts(title=main)
    sdiff <- survdiff(eval(sfit$call$formula), data=eval(sfit$call$data));
    if(pvaltxt=='auto'){
      pval <- pchisq(sdiff$chisq, length(sdiff$n)-1, lower.tail=F)
      pvaltxt <- ifelse(pval < 0.0001, 'P < 0.0001',
                        paste('P =', signif(pval,2)))
    }
    p <- p+ geom_point(data=subset(df, n.censor>0), aes(time, surv,color=strata),shape=censorshape,size=censorsize);
    
#     Create a blank plot for place-holding
    blank_pic<- ggplot(df, aes(time,surv))+geom_blank()+theme_bw()+
      #opts(axis.text.x=theme_blank(), axis.text.y=theme_blank(),
      theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.ticks=element_blank(), panel.grid.major=element_blank(),
        panel.border=element_blank())

if(table) {
#     Create table graphic to include at-risk numbers
    risk.data  <- data.frame(strata = summary(sfit, times=times,
                    extend=T)$strata,
            time = summary(sfit, times=times, extend=T)$time,
            n.risk = summary(sfit, times=times, extend=T)$n.risk)
    data_table  <-  ggplot(risk.data, aes(x=time, y=strata,
                        label=format(n.risk, nsmall=0)))+#, color=strata))+
                geom_text(size=3.5)+
                theme_bw()+
                scale_y_discrete(breaks=as.character(levels(risk.data$strata)),
                labels=ystratalabs)+
#                 scale_y_discrete(#format1ter = abbreviate,
#                   breaks=1:3,
#                   labels=ystratalabs)+
                scale_x_continuous('Numbers at risk', limits=c(0, max(sfit$time)))+
                ## opts(axis.title.x=theme_text(size=10, vjust=1))+
                ## opts(panel.grid.major   =   theme_blank())+
                ## opts(panel.grid.minor   =   theme_blank())+
                ## opts(panel.border       =   theme_blank())+
                ## opts(axis.text.x        =   theme_blank())+
                ## opts(axis.ticks         =   theme_blank())+
                ## opts(axis.text.y       =   theme_text(face='bold', hjust=1))
                theme(axis.title.x=element_text(size=10, vjust=1))+
                theme(panel.grid.major   =   element_blank())+
                theme(panel.grid.minor   =   element_blank())+
                theme(panel.border       =   element_blank())+
                theme(axis.text.x        =   element_blank())+
                theme(axis.ticks         =   element_blank())+
                theme(axis.text.y       =   element_text(face='bold', hjust=1))

    #data_table  <- data_table + opts(legend.position='none')+
    data_table  <- data_table + theme(legend.position='none')+
                                        xlab(NULL)+ylab(NULL)

    #data_table <- data_table+opts(plot.margin=unit(c(-1.5, 1, 0.1, ifelse(m<10,2.5,3.5)-0.28*m),'lines'))
    data_table <- data_table+theme(plot.margin=unit(c(-1.5, 1, 0.1, ifelse(m<10,2.5,3.5)-0.28*m),'lines'))


#     Plotting the graphs
    p <- ggplotGrob(p)
    p <- addGrob(p, textGrob(x=unit(.8, 'npc'), y=unit(.25,'npc'), label=pvaltxt,
      gp=gpar(fontsize=12)))
    grid.arrange(p, blank_pic, data_table, clip=F, nrow=3, ncol=1,
      heights = unit(c(2,.1,.25), c('null','null','null')))
    if(returns) {
      a <- arrangeGrob(p, blank_pic, data_table, clip=F, nrow=3, ncol=1,
        heights = unit(c(2, .1, .25), c('null','null','null')))
      return(a)
    }
} else {
  #browser();
  #p <- ggplotGrob(p)
  #p <- addGrob(p, textGrob(x=unit(0.5,'npc'), y=unit(0.23,'npc'),
  #  label=pvaltxt, gp=gpar(fontsize=12)))
  grid.arrange(p)
  if(returns) return(p)
}
}
