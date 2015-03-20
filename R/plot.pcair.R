plot.pcair <- function(x, vx=1, vy=2, pch=NULL, col=NULL, xlim=NULL, ylim=NULL, main=NULL, sub=NULL, xlab=NULL, ylab=NULL, ...){
	xlab <- if(is.null(xlab)){ paste("PC-AiR PC",vx) }else{ xlab }
	ylab <- if(is.null(ylab)){ paste("PC-AiR PC",vy) }else{ ylab }
	xy <- xy.coords(x$vectors[,vx], x$vectors[,vy], xlab=xlab, ylab=ylab)
	xlim <- if(is.null(xlim)){ range(xy$x[is.finite(xy$x)]) }else{ xlim }
	ylim <- if(is.null(ylim)){ range(xy$y[is.finite(xy$y)]) }else{ ylim }
	main <- if(is.null(main)){ paste("Plot of PC-AiR PC", vx, "vs. PC", vy) }else{ main }
	dev.hold()
	on.exit(dev.flush())
	plot.new()
	plot.window(xlim=xlim, ylim=ylim, ...)
	if(is.null(pch)){
        IDs <- rownames(x$vectors)
        pch <- rep(20,x$nsamp); pch[which(IDs %in% x$rels)] <- 3
	}
	if(is.null(col)){
        IDs <- rownames(x$vectors)
		col <- rep(1,x$nsamp); col[which(IDs %in% x$rels)] <- 4
	}
	plot.xy(xy, type="p", pch=pch, col=col, ...)
	Axis(x$vectors[,vx], side=1, ...)
	Axis(x$vectors[,vy], side=2, ...)
	box(...)
	title(main=main, sub=sub, xlab=xlab, ylab=ylab, ...)
}

