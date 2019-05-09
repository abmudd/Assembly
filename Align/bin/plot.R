#!/usr/bin/env Rscript
PROGRAM = 'plot.R'
DESCRIP = 'Simple plotting in R'
VERSION = '0.1'

# This script is derived from a script written by Jessen Bredeson <jessenbredeson@berkeley.edu>.


args = strsplit(commandArgs(TRUE),' ');
if (! isTRUE(as.logical(length(args)))) {
    cat('',"\n");
    cat('Program: ',PROGRAM,' (',DESCRIP,")\n", sep='');
    cat('Version: ',VERSION,"\n",  sep='');
    cat('Usage:   ',PROGRAM,' <command> <arguments>',"\n\n",sep='');
    cat('Command: hexbin    Hexbin',"\n");
    cat('         hist      Histogram',"\n");
    cat('         hist2d    2D histograms',"\n");
    cat('         line      Line plot',"\n");
    cat('         scatter   Scatter plot',"\n\n");
    q(status=1);
}


if (isTRUE(as.logical(match(args[[1]],'scatter')))) {
    if (length(args) < 4) {
        cat('Usage: ',PROGRAM,' scatter <in.table> <X.column> <Y.column>',"\n",sep='');
        q(status=1);
    }
    dat = read.table(args[[2]],check.names=T,row.names=NULL);
    x = as.integer(args[[3]]);
    y = as.integer(args[[4]]);
    o = args[[2]]
    if ((ncol(dat) < x) | (ncol(dat) < y)) {
        cat('[',PROGRAM,'] Column out of bounds',"\n",sep='');
        q(status=1);
    }
    pdf(paste(o,'.scatter.pdf',sep=''));
    xlbl = names(dat)[x];
    ylbl = names(dat)[y];
    if (isTRUE(as.logical(match(names(dat)[1],'V1')))) {
        xlbl = paste('Field',x);
        ylbl = paste('Field',y);
    }
    m = lm(dat[,y] ~ dat[,x]);
    plot(dat[,x],dat[,y],type="p", pch='+', 
        main=paste(args[[2]],' (N = ',length(dat[,1]),')',sep=''), 
        xlab=xlbl, ylab=ylbl, 
        sub=paste('Y = ',m$coefficients[1],' + ',m$coefficients[2],'X'));
    abline(a=m$coefficients[1],b=m$coefficients[2],col="red");
    suppressMessages(dev.off());

} else if (isTRUE(as.logical(match(args[[1]],'line')))) {
    if (length(args) < 4) {
        cat('Usage: ',PROGRAM,' line <in.table> <X.column> <Y.column>',"\n",sep='');
        q(status=1);
    }
    dat = read.table(args[[2]],check.names=T,row.names=NULL);
    x = as.integer(args[[3]]);
    y = as.integer(args[[4]]);
    o = args[[2]]
    if ((ncol(dat) < x) | (ncol(dat) < y)) {
        cat('[',PROGRAM,'] Column out of bounds',"\n",sep='');
	q(status=1);
    }
    pdf(paste(o,'.line.pdf',sep=''));
    xlbl = names(dat)[x];
    ylbl = names(dat)[y];
    if (isTRUE(as.logical(match(names(dat)[1],'V1')))) {
        xlbl = paste('Field',x);
        ylbl = paste('Field',y);
    }
    m = lm(dat[,y] ~ dat[,x]);
    plot(dat[,x],dat[,y],type="l",
        main=paste(args[[2]],' (N = ',length(dat[,1]),')',sep=''), 
        xlab=xlbl, ylab=ylbl, 
        sub=paste('Y = ',m$coefficients[1],' + ',m$coefficients[2],'X'));
    suppressMessages(dev.off());

} else if (isTRUE(as.logical(match(args[[1]],'hist')))) {
    if (length(args) < 4) {
        cat('Usage: ',PROGRAM,' hist <in.table> <X.column> <n.bins>',"\n",sep='');
        q(status=1);
    } 
    dat = read.table(args[[2]],check.names=T,row.names=NULL);
    x = as.integer(args[[3]]);
    b = as.integer(args[[4]]);
    o = args[[2]];
    if (ncol(dat) < x) {
       cat('[',PROGRAM,'] Column index out of bounds',"\n",sep='');
       q(status=1);
    }
    if (b == 0) {
        b = 'Sturges';
    }
    xlbl = names(dat)[x];
    if (isTRUE(as.logical(match(names(dat)[1],'V1')))) {
        xlbl = paste('Field',x);
    }
    pdf(paste(o,'.hist.pdf',sep=''));
    hist(dat[,x], breaks=b, col="grey80", 
        xlab=xlbl, main=paste(args[[2]],' (N = ',length(dat[,1]),')',sep=''));
    suppressMessages(dev.off());

} else if (isTRUE(as.logical(match(args[[1]],'hist2d')))) {
    if (length(args) < 5) {
        cat('Usage: ',PROGRAM,' hist2d <in.table> <X.column> <Y.column> <n.bins>',"\n",sep='');
        q(status=1);
    }
    packages = c("RColorBrewer","MASS")
    package.check <- lapply(packages, FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE)
            library(x, character.only = TRUE)
        }
    })
    dat = read.table(args[[2]],check.names=T,row.names=NULL);
    x = as.integer(args[[3]]);
    y = as.integer(args[[4]]);
    b = as.integer(args[[5]]);
    o = args[[2]]
    if ((ncol(dat) < x) | (ncol(dat) < y)) {
        cat('[',PROGRAM,'] Column out of bounds',"\n",sep='');
        q(status=1);
    }
    pdf(paste(o,'.hist2d.pdf',sep=''));
    xlbl = names(dat)[x];
    ylbl = names(dat)[y];
    if (isTRUE(as.logical(match(names(dat)[1],'V1')))) {
        xlbl = paste('Field',x);
        ylbl = paste('Field',y);
    }
    rf <- colorRampPalette(c("royalblue","springgreen","yellow","red"))
    h1 <- hist(dat[,x], breaks=b, plot=F)
    h2 <- hist(dat[,y], breaks=b, plot=F)
    top <- max(h1$counts, h2$counts)
    k <- kde2d(dat[,x], dat[,y], n=b)
    par(mar=c(3,3,1,1))
    layout(matrix(c(2,0,1,3),2,2,byrow=T),c(3,1), c(1,3))
    image(k, col=rf(50))
    par(mar=c(0,2,1,0))
    barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col="grey80")
    par(mar=c(2,0,0,1))
    barplot(h2$counts, axes=F, xlim=c(0, top), space=0, horiz=T, col="grey80")
    mtext(paste(args[[2]],' (N = ',length(dat[,1]),')',sep=''), side=3, outer=TRUE, line=-3)
    mtext(xlbl, side=1, outer=TRUE, line=-1.1, cex=0.8)
    mtext(ylbl, side=2, outer=TRUE, line=-1.1, cex=0.8)
    suppressMessages(dev.off())

} else if (isTRUE(as.logical(match(args[[1]],'hexbin')))) {
    if (length(args) < 5) {
        cat('Usage: ',PROGRAM,' hexbin <in.table> <X.column> <Y.column> <n.bins>',"\n",sep='');
        q(status=1);
    }
    packages = c("hexbin")
    package.check <- lapply(packages, FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
	    install.packages(x, dependencies = TRUE)
	    library(x, character.only = TRUE)
	}
    })
    dat = read.table(args[[2]],check.names=T,row.names=NULL);
    x = as.integer(args[[3]]);
    y = as.integer(args[[4]]);
    b = as.integer(args[[5]]);
    o = args[[2]]
    if ((ncol(dat) < x) | (ncol(dat) < y)) {
        cat('[',PROGRAM,'] Column out of bounds',"\n",sep='');
        q(status=1);
    }
    pdf(paste(o,'.hexbin.pdf',sep=''));
    xlbl = names(dat)[x];
    ylbl = names(dat)[y];
    if (isTRUE(as.logical(match(names(dat)[1],'V1')))) {
        xlbl = paste('Field',x);
        ylbl = paste('Field',y);
    }
    bin<-hexbin(dat[,x], dat[,y], xbins=b)
    plot(bin, 
        main=paste(args[[2]],' (N = ',length(dat[,1]),')',sep=''), 
        xlab=xlbl, ylab=ylbl);
    suppressMessages(dev.off());

} else {
    cat('',"\n");
    cat('Program: ',PROGRAM,' (',DESCRIP,")\n", sep='');
    cat('Version: ',VERSION,"\n",  sep='');
    cat('Usage:   ',PROGRAM,' <command> <arguments>',"\n\n",sep='');
    cat('Command: hexbin    Hexbin',"\n");
    cat('         hist      Histogram',"\n");
    cat('         hist2d    2D histograms',"\n");
    cat('         line      Line plot',"\n");
    cat('         scatter   Scatter plot',"\n\n");
    q(status=1);
}
