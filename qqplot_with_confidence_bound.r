qqplot_conf=function(P_T=NULL, mainT="QQ-Plot"){

    if ( is.null(P_T) ){
        stop('Invalid p-value list!')
    }
    tn=length(P_T)
    P_T_l = sort(-log10(P_T),decreasing=TRUE)
    #P_S_l = sort(-log10(P_S),decreasing=TRUE)
    x = seq(1/tn,1-1/tn, length.out=tn-1)
    x_ub = x + 1.96*sqrt(x*(1-x)/(tn-1))
    x_lb = x - 1.96*sqrt(x*(1-x)/(tn-1))
    x = -log10(x)
    x_ub = -log10(x_ub)
    x_lb = -log10(x_lb)
    plot(x,x, xlab="-log10(expected)", ylab="-log10(p-value)", main=mainT)
    lines(x, x_ub,col="grey")
    lines(x, x_lb,col="grey")
    lines(x, P_T_l[2:tn],col="orange" ,lwd=2)
    #lines(x, P_S_l[2:tn],col="red" ,lwd=2)

}
