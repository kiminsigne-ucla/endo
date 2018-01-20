zcrosscorr <- function(inputfile,filetype,srange=c(50,500),bin=5,acceptflag=F){
  data <- read.aligned.tags(inputfile,filetype)
  zcc <- get.binding.characteristics(data=data,srange=srange,bin=bin,accept.all.tags=acceptflag)
  print(paste("binding peak separation distance =",zcc$peak$x))
  jfile <- inputfile
  print_zcc(zccprofile=zcc,outfile=jfile)
  return(zcc)
}

print_zcc <- function(zccprofile,outfile){
  jpeg(paste(outfile,"_CCPLOT.jpeg",sep=""),width=960,height=960)
  par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8)
  plot(zccprofile$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation")
  abline(v=zccprofile$peak$x,lty=2,col=2)
  dev.off()
}

read.aligned.tags <- function(filename,filetype) {
  tl <- lapply(.Call("read_aligned_tags",filename,filetype),function(d) {
    xo <- order(abs(d$t));
    d$t <- d$t[xo];
    d$n <- d$n[xo];
    return(d);
  });
  # separate tags and quality
  return(list(tags=lapply(tl,function(d) d$t),quality=lapply(tl,function(d) d$n)));
}

get.binding.characteristics <- function(data,srange=c(50,500),bin=5,cluster=NULL,min.tag.count=1e3,acceptance.z.score=3,remove.tag.anomalies=T,anomalies.z=5,accept.all.tags=F){
  if(remove.tag.anomalies) {
    data <- remove.tag.anomalies(data,z=anomalies.z);
  }
  # take highest quality tag bin
  if(!is.null(data$quality) & !accept.all.tags) {
    min.bin <- min(unlist(lapply(data$quality,min)))
    chrl <- names(data$tags); names(chrl) <- chrl;
    otl <- lapply(chrl,function(chr) data$tags[[chr]][data$quality[[chr]]==min.bin]);
  } else {
    otl <- data$tags;
  }
  # remove empty chromosomes
  otl <- otl[unlist(lapply(otl,length))!=0];

  # calculate strand scc
  if(!is.null(cluster)) {
    cc <- clusterApplyLB(cluster,otl,tag.scc,srange=srange,bin=bin);
    names(cc) <- names(otl); 
  } else {
    cc <- lapply(otl,tag.scc,srange=srange,bin=bin);
  }
  ccl<-list(sample=cc);
  ccl.av <- lapply(names(ccl),t.plotavcc,type='l',ccl=ccl,return.ac=T,ttl=list(sample=otl),plot=F)[[1]]
  ccl.av <- data.frame(x=as.numeric(names(ccl.av)),y=as.numeric(ccl.av));
  
  # find peak
  pi <- which.max(ccl.av$y);
  
  # determine width at third-height
  th <- (ccl.av$y[pi]-ccl.av$y[length(ccl.av$y)])/3+ccl.av$y[length(ccl.av$y)]
  whs <- max(ccl.av$x[ccl.av$y>=th]);

  # determine acceptance of different quality bins
  # calculates tag scc for the best tags, and combinations of best tag category with every other category
  # for subsequent selection of acceptable categories
  scc.acceptance.calc <- function() {

    qr <- range(unlist(lapply(data$quality,range)))

    # start with best tags

    # determine half-width for scc calculations
    pi <- which.max(ccl.av$y);

    # determine width at half-height
    th <- (ccl.av$y[pi]-ccl.av$y[length(ccl.av$y)])/2+ccl.av$y[length(ccl.av$y)]
    lwhs <- max(ccl.av$x[ccl.av$y>=th])-ccl.av$x[pi];
    lwhs <- max(c(20,bin*10,lwhs));
    srange <- ccl.av$x[pi]+c(-lwhs,lwhs)

    # calculate chromosome-average scc
    t.scc <- function(tags) {
      if(is.null(cluster)) {
        cc <- lapply(tags,tag.scc,srange=srange,bin=bin);
      } else {
        cc <- clusterApplyLB(cluster,tags,tag.scc,srange=srange,bin=bin); names(cc) <- names(tags);
      }
      return(t.plotavcc(1,type='l',ccl=list(cc),ttl=list(tags),plot=F,return.ac=T))
    }

    # returns info list for a given tag length (lv), mismatch count (nv)
    t.cat <- function(qual) {
      # construct tag set
      if(qual==qr[1]) {
        ts <- otl;
      } else {
        nts <- names(otl); names(nts) <- nts;
        # select tags
        at <- lapply(nts,function(chr) data$tags[[chr]][data$quality[[chr]]==qual]);
        ntags <- sum(unlist(lapply(at,length)));
        if(ntags<min.tag.count) { return(NULL); }
        # append to otl
        ts <- lapply(nts,function(nam) c(otl[[nam]],at[[nam]]));
      }
      return(t.scc(ts));
    }
    # calculate cross-correlation values for each quality bin
    ql <- sort(unique(unlist(lapply(data$quality,unique)))); names(ql) <- ql;
    qccl <- lapply(ql,t.cat);
    # acceptance tests
    ac <- c(T,unlist(lapply(qccl[-1],function(d) if(is.null(d)) { return(F) } else { t.test(d-qccl[[as.character(min.bin)]],alternative="greater")$p.value<pnorm(acceptance.z.score,lower.tail=F) }))); names(ac) <- names(qccl);
    return(list(informative.bins=ac,quality.cc=qccl))
  }

  if(accept.all.tags | is.null(data$quality)) {
    return(list(cross.correlation=ccl.av,peak=list(x=ccl.av$x[pi],y=ccl.av$y[pi]),whs=whs))    
  } else {
    acc <- scc.acceptance.calc();
    return(list(cross.correlation=ccl.av,peak=list(x=ccl.av$x[pi],y=ccl.av$y[pi]),whs=whs,quality.bin.acceptance=acc));
  }
}

tag.scc <- function(tags,srange=c(50,250),bin=1,tt=NULL,llim=10) {
  if(is.null(tt)) {
    tt <- table(sign(tags)*as.integer(floor(abs(tags)/bin+0.5)));
  }
  if(!is.null(llim)) { l <- mean(tt); tt <- tt[tt<llim*l] }
  tc <- as.integer(names(tt));
  tt <- as.numeric(tt);

  pv <- tt; pv[tc<0]<-0;
  nv <- tt; nv[tc>0]<-0;

  pti <- which(tc>0)
  nti <- which(tc<0);

  ptc <- tc[pti];
  ntc <- (-1)*tc[nti];

  ptv <- tt[pti];
  ntv <- tt[nti];

  trng <- range(c(range(ptc),range(ntc)))
  l <- diff(trng)+1;
  rm(tc,tt);

  mp <- sum(ptv)*bin/l;   mn <- sum(ntv)*bin/l;
  ptv <- ptv-mp; ntv <- ntv-mn;
  ss <- sqrt((sum(ptv*ptv)+(l-length(ptv))*mp^2) * (sum(ntv*ntv)+(l-length(ntv))*mn^2));

  t.cor <- function(s) {
    smi <- match(ptc+s,ntc);
    return((sum(ptv[!is.na(smi)]*ntv[na.omit(smi)]) -
           mn*sum(ptv[is.na(smi)]) -
           mp*sum(ntv[-na.omit(smi)]) +
           mp*mn*(l-length(ptv)-length(ntv)+length(which(!is.na(smi)))))/ss);
  }
  shifts <- floor(seq(srange[1],srange[2],by=bin)/bin+0.5);
  scc <- unlist(lapply(shifts,t.cor)); names(scc) <- shifts*bin;
  return(scc);
}

################################################################################################
################################## OTHER ROUTINES ##############################################
################################################################################################

# plot tag cross-correlation
t.plotcc <- function(ac, lab=c(10,5,7), ylab="correlation", xlab="lag", pch=19, grid.i=c(-5:5), grid.s=10, type='b', plot.grid=F, cols=c(1,2,4,"orange",8,"pink"), min.peak.x=NULL, xlim=NULL, plot.147=F, plot.max=T, rmw=1, rescale=F, legendx="right", ltys=rep(1,length(ac)), ...) {
  if(is.list(ac)) {
    cols <- cols[1:length(ac)];  
    if(!is.null(xlim)) {
      vx <- as.numeric(names(ac[[1]])); vx <- which(vx>=xlim[1] & vx<=xlim[2]);
      ac[[1]] <- (ac[[1]])[vx];
    }else{
      xlim <- range(as.numeric(names(ac[[1]])));
    }
  
    plot(as.numeric(names(ac[[1]])),runmean(ac[[1]],rmw),type=type,pch=pch,xlab=xlab,ylab=ylab,lab=lab, col=cols[1], xlim=xlim, lty=ltys[1], ...);
    if(length(ac)>1) {
      for(i in seq(2,length(ac))) {
        irng <- range(ac[[i]]);
        vx <- as.numeric(names(ac[[i]])); vx <- which(vx>=xlim[1] & vx<=xlim[2]);
        if(rescale) {
          lines(as.numeric(names(ac[[i]])[vx]),runmean((ac[[i]][vx]-irng[1])/diff(irng)*diff(range(ac[[1]]))+min(ac[[1]]),rmw),col=cols[i],lty=ltys[i]);
        } else {
          lines(as.numeric(names(ac[[i]]))[vx],runmean(ac[[i]][vx],rmw),col=cols[i],lty=ltys[i]);
        }
      }
    }
    if(is.null(min.peak.x)) {
      m <- as.numeric(names(ac[[1]])[which.max(ac[[1]])]);
    } else {
      sac <- (ac[[1]])[which(as.numeric(names(ac[[1]]))>min.peak.x)]
      m <- as.numeric(names(sac)[which.max(sac)]);
    }
    legend(x="topright",bty="n",legend=c(names(ac)),col=cols,lty=ltys)
  } else {
    if(!is.null(xlim)) {
      vx <- as.numeric(names(ac));
      vx <- which(vx>=xlim[1] & vx<=xlim[2]);
      ac <- ac[vx];
    } else {
      xlim <- range(as.numeric(names(ac)));
    }
    
    plot(names(ac),runmean(ac,rmw),type=type,pch=pch,xlab=xlab,ylab=ylab,lab=lab, xlim=xlim, ...);
    if(is.null(min.peak.x)) {
      m <- as.numeric(names(ac)[which.max(ac)]);
    } else {
      sac <- ac[which(names(ac)>min.peak.x)]
      m <- as.numeric(names(sac)[which.max(sac)]);
    }
  }
  if(plot.147) {
    abline(v=147,lty=2,col=8);
  }
  if(plot.grid) {
    abline(v=m+grid.i*grid.s,lty=3,col="pink");
  }
  if(plot.max) {
    abline(v=m,lty=2,col=2);
    legend(x=legendx,bty="n",legend=c(paste("max at ",m,"bp",sep="")));
    return(m);
  }
}

# plot chromosome-acerage cross-correlation 
t.plotavcc <- function(ci, main=paste(ci,"chromosome average"), ccl=tl.cc, return.ac=F, ttl=tl, plot=T, ... ) {
  cc <- ccl[[ci]];
  if(length(cc)==1)  { return(cc[[1]]) };
  if(length(cc)==0) { return(c()) };
  ac <- do.call(rbind,cc);
  # omit NA chromosomes
  ina <- apply(ac,1,function(d) any(is.na(d)));
  
  tags <- ttl[[ci]]; 
  avw <- unlist(lapply(tags,length));    avw <- avw/sum(avw);    
  ac <- ac[!ina,]; avw <- avw[!ina];
  ac <- apply(ac,2,function(x) sum(x*avw));
  if(plot) {
    m <- t.plotcc(ac, main=main, ...);
    if(!return.ac) { return(m) }
  }
  if(return.ac) { return(ac) }
}

# removes tag positions that have anomalously high counts on both strands
# z - z-score used to determine anomalous bins
# zo - z used to filter out one-strand matches
# trim.fraction - fraction of top bins to discard when calculating overall background density
remove.tag.anomalies <- function(data, bin=1,trim.fraction=1e-3,z=5,zo=3*z) {
  
  t.remove.tag.anomalies <- function(tv,bin=1,trim.fraction=1e-3,z=5,zo=3*z,return.indecies=F) {
    tt <- table(floor(tv/bin));

    # trim value
    stt <- sort(as.numeric(tt));
    stt <- stt[1:(length(stt)*(1-trim.fraction))];
    mtc <- mean(stt); tcd <- sqrt(var(stt));

    thr <- max(1,ceiling(mtc+z*tcd));
    thr.o <- max(1,ceiling(mtc+zo*tcd));
    # filter tt
    tt <- tt[tt>=thr]
    # get + and - tags
    tp <- as.numeric(names(tt));
    pti <- tp>0;
    it <- intersect(tp[pti],(-1)*tp[!pti]);
    # add one-strand matches
    it <- unique(c(it,tp[tt>=thr.o]));
    sit <- c(it,(-1)*it);
    
    if(bin>1) {
      sit <- sit*bin;
      sit <- c(sit,unlist(lapply(1:bin,function(i) sit+i)))
    }
    if(return.indecies) {
      return(!tv %in% sit);
    } else {
      return(tv[!tv %in% sit]);
    }
  }

  vil <- lapply(data$tags,t.remove.tag.anomalies,return.indecies=T,bin=bin,trim.fraction=trim.fraction,z=z,zo=zo);
  chrl <- names(data$tags); names(chrl) <- chrl;
  data$tags <- lapply(chrl,function(chr) data$tags[[chr]][vil[[chr]]]);
  # count tags to remove empty chromosomes
  nt <- unlist(lapply(data$tags,length));
  if(any(nt==0)) {
    data$tags <- data$tags[nt!=0]
  }
  
  if(!is.null(data$quality)) {
    data$quality <- lapply(chrl,function(chr) data$quality[[chr]][vil[[chr]]]);
    data$quality <- data$quality[nt!=0];
  }
  if(!is.null(data$names)) {
    data$names <- lapply(chrl,function(chr) data$names[[chr]][vil[[chr]]]);
    data$names <- data$names[nt!=0];
  }
  
  return(data);
}