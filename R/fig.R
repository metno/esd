#' @name fig
#' @description R-function for producing plots of downscaled ensembles (empirical-statistical downscaling, ESD)
#' @param x a list object which has several elements, each of which contain an 
#' ensemble of downscaled results for a specific CMIP run (e.g. ssp370)
#' @param loc name of location to plot
#' @param obs a station object to replace observations used for calibrating the method. E.g. an updated station series.
#' @param it index time - annual or seasonal - c('annual','DJF', 'MAM', 'JJA', 'SON')
#' @param lev - levels to show in terms of quantiles. 0 is minimum and maximum values 
#' (also use 1 - lev), whereas 0.01 produces the range between the 1-99% percentiles. Several 
#' levels are shown as layers of darker shading.
#' @param ylim - see \code{\link{plot}}
#' @param xlim - see \code{\link{plot}}
#' @param highlight - names of GCMs to highlight
#' @param scenarios names of scenarios to show
#' @param angle - see \code{\link{polygon}}
#' @param density - see \code{\link{polygon}}
#' @param legendpos position of the legend - see \code{\link{legend}}
#' @param col.scen colour of the plumes.
#' @param threshold a threshold used to compare with percentiles from the ensemble spread
#' @param verbose Logical value defaulting to FALSE. If FALSE, do not display
#' comments (silent mode). If TRUE, displays extra information on progress.
#' @return A time series of "zoo" "station" class with additional attributes
#' used for further processing.
#'
#' @author Rasmus Benestad
#' @examples
#' 
#' data(tas)
#' fig(tas)
#' # Highlight NorESM
#' fig(tas,highlight='NorESM')
#' 

# Function to convert hex color to rgb format with optional alpha
#' @export 
hex_to_rgb <- function(hex, alpha = 1,verbose=FALSE) {
  # Remove the '#' character if present
  if (verbose) print(hex)
  hex <- gsub("#", "", hex)
  # Convert hex components to decimal and scale to [0, 1]
  red <- strtoi(substr(hex, 1, 2), base = 16) / 255
  green <- strtoi(substr(hex, 3, 4), base = 16) / 255
  blue <- strtoi(substr(hex, 5, 6), base = 16) / 255
  # Generate the rgb color with alpha
  if (verbose) print(c(red, green, blue, alpha))
  return(rgb(red, green, blue, alpha))
}

## function that generates plumes:
## max-min, q5-q95, q10-q90, q25-q75
#' @export plume
plume <- function(x,t,col="#255c36",alpha=0.2,
                  levs=c(0,0.05,0.10,0.25),density=NULL,angle=NULL,verbose=FALSE) {
  require(RColorBrewer)
  require(ggplot2)
  
  ## 
  if (verbose) print(match.call())
  mu <- apply(x,2,'mean')
  valid.data <- is.finite(mu)
  if (length(mu)!= length(t)) browser()
  if (verbose) print(paste('plume:',sum(valid.data)))
  mus <- smooth.spline(t[valid.data], mu[valid.data])$y
  
  ## Lower and upper quantile
  for (lev in levs) { 
    ql <- apply(x,2,quantile,probs=lev)
    qu <- apply(x,2,quantile,probs=1-lev)
    qls <- smooth.spline(t[valid.data], ql[valid.data])$y
    qus <- smooth.spline(t[valid.data], qu[valid.data])$y
    polygon(c(t[valid.data],rev(t[valid.data])),
            c(qls, rev(qus)),col=hex_to_rgb(col,alpha = alpha),
            angle=angle, density=density, border=NA)
  }
  lines(t,mus,lwd=3,col=col)
}

#' @export fig
fig <- function(x,loc="OSLO BLINDERN",it='annual',main=NULL,obs=NULL,
                levs=c(0,0.05,0.10,0.25),showobstype="b",param=NULL,
                ylim=NULL,xlim=c(1950,2100),highlight=NULL,ylab=NULL,
                scenarios=c('ssp370','rcp45','ssp245','ssp585'),
                angle=NULL, density=NULL,legendpos="bottomright",
                col.scen=c("#255c36","#54ab54","#b9dabb",
                           "#8b316b","#6a5102","#ea5175","#f49857","#e4ca54",
                           "#0c465f","#238582","#007cc1","#8db7cb"),
                bty='n',threshold=NULL,verbose=FALSE) {
  
  require(RColorBrewer)
  require(ggplot2)
  
  allN <- function(x,verbose=FALSE) {
    ## Remove members with short time span
    nt <- length(index(x))
    cls <- class(x)
    if (is.null(dim(x))) {
      x <- combine.stations(x,x)
      attr(x,'location') <- paste0(loc(x),c('','2'))
    }
    if (verbose) {print('Remove members with short time span:'); print(dim(x))}
    nv <- apply(coredata(x),2,FUN='nv')
    if (verbose) print(nv)
    im <- nv >= quantile(nv,probs = 0.75)
    x <- subset(x,im=im,verbose=verbose)
    ## Remove times with incomplete ensemble
    if (verbose) {print('Select interval with complete ensemble:'); print(dim(x))}
    nv <- apply(coredata(x),1,FUN='nv')
    if (verbose) print(nv)
    it <- nv==max(nv,na.rm=TRUE)
    y <- subset(x,it=it,verbose=verbose)
    y <- attrcp(x,y)
    class(y) <- cls
    if (sum(it) != dim(y)[1]) browser()
    ## Adjust the offset of the ensemble to match recent years of the observations
    if (verbose) print('Adjust offset')
    obs <- attr(x,'station')
    yrs <- sort(year(obs),decreasing = TRUE)
    i1 <- is.element(year(obs),year(x))
    i2 <- is.element(year(x),year(obs))
    if (verbose) cat(range(yrs,sum(i1),sum(i2),'\n'))
    mobs <- round(mean(coredata(obs[i2])),2)
    if (verbose) print(mobs)
    for (i in 1:dim(x)[2]) x[,i] <- x[,i] - mean(coredata(x[i1,i])) + mobs
    attr(y,"model_id") <- attr(x,"model_id")[im]
    attr(y,"selected models") <-im
    if (verbose) cat('finished allN - exit \n')
    return(y)
  }
  
  fix.model.id <- function(x) {
    ## Fix missing rip
    y <- x[grep('rip',x)]
    gcms <- rownames(table(sub('.rip','',y)))
    print(gcms)
    for (gcm in gcms) {
      im <- grep(gcm,x)
      for (i in 1:length(im)) x[im[i]] <- sub('rip',paste0('r',i,'i1p1f2'),x[im[i]])
    }
    return(x)
  }
  
  
  
  
  if (verbose) print(match.call())
  #if (is.na(threshold)) threshold <- NULL
  if (length(grep('RCM',names(x)))>0) {
    ## The data is for regional series, so Fig(x) should be used instead
    print('Region curves - use Fig()')
    eval(parse(text=sub('fig','Fig',deparse(match.call()))))
    return()
  }
  iesd <- NULL
  ## Match colour with scenario
  ## Extract the variable name from the list name
  if (is.null(param)) param <- sub('\\..*','',names(x)[1])
  ## Extract the season from the list name
  if (!is.null(it)) { 
    wrongseas <- names(x)[-grep(it,names(x))]
    for (omit in wrongseas) x[[omit]] <- NULL
  }
  if (verbose) print(names(x))
  if (is.null(names)) browser()
  
  ## The locations are stored in separate list elements
  select <- names(x[[1]])[grep(tolower(loc),tolower(names(x[[1]])))]
  
  if (length(select)!=1) {
    print(loc); print(names(x[[1]])); print(select)
    print('fig: Detected zero or several locations with the same name!')
    if (length(select)>1) select <- select[1] else {
      select <- 1
    }
  }
  if (verbose) {print(select); print(names(x[[1]]))}
  x1 <- allN(x[[1]][[select]],verbose=verbose) 
  if (is.null(obs)) y <- attr(x1,'station') else y <- obs
  index(y) <- year(y)
  variable <- switch(param,'tas'='temperature','tam'='temperature','t2m'='temperature',
                     'rr'='precipitation','precip'='precipitation','obs'='observation',
                     'fw'='wet-day frequency','mu'='wet-day mean precipitation',
                     'tp'='precipitation')
  if (verbose) cat(variable,'\n')
  if (is.null(ylab)) {
    ylab <- switch(param,'tas'=expression(degree*C),'tam'=expression(degree*C),
                   't2m'=expression(degree*C),'rr'='mm','precip'='mm','tp'='mm',
                   'fw'=expression(f[w]),'mu'=expression(mu))
  }
  if (is.null(it)) main <- paste(variable,'at',select)
  if (is.null(main)) {
    main <- switch(it,'annual'='Yearly','DJF'='Winter',
                   'MAM'='spring','JJA'='Summer','SON'='Autumn')
    main <- paste(main,variable,'at',select)
  } 

  if (is.null(ylim)) 
    ylim <- range(unlist(lapply(x,function(y) y[[select]])),na.rm=TRUE) * c(0.8,1.2)
  
  plot(index(y),coredata(y),type=showobstype,pch=19,lwd=2,xlim=xlim,ylim=ylim,
       main=main,bty=bty,xlab='',ylab=ylab)
  grid()
  if (verbose) print(range(index(y)))
  
  ## Number of model runs
  NM <- c(); descr <- c(); cols <- c()
  ## Store percentage of ensemble members above a threshold if a given threshold
  qp <- list()
  ## Loop through scenarios
  if (length(names(x[[1]]))==0) {print("Problem in fig()!"); browser()}
  for (scen in names(x)) { 
    sce <- sub('\\..*','',sub(paste0(param,'.'),'',scen,fixed=TRUE))
    if (verbose) cat('sce:',sce,' \n')
    if (sce=='obs') sce <- 'ssp370'
    if (length(grep(sce,scenarios))>0) {  
      if (verbose) print(scen)
      nm <- dim(x[[scen]][[select]])[2]
      t <- xlim[1]:xlim[2]; nt <- length(t)
      Z <- matrix(rep(NA,nt*nm),nm,nt)
      i <- 1
      NM <- c(NM,nm); descr <- c(descr,sce)
      xx <- x[[scen]]
      for (im in 1:nm) {
        z <- xx[[select]]
        z <- z[,im]
        index(z) <- year(z)
        #print(c(range(round(coredata(z))),range(index(z))))
        #lines(index(z),coredata(z),col="grey")
        i1 <- is.element(t,index(z))
        i2 <- is.element(index(z),t)
        if (sum(i1) !=sum(i2)) browser("sum(i1) !=sum(i2)")
        Z[i,i1] <- coredata(z)[i2]; i <- i + 1
      }
      
      ## Code crashes for only one station - so duplicate it
      if (is.null(dim(Z))) Z <- combine.stations(Z,Z)
      nv <- apply(Z,1,'nv') 
      Z <- Z[nv==max(nv),]
      NM[length(NM)] <- sum(nv==max(nv))
      if (sum(!is.finite(Z))>0) {
        nv <- apply(Z,2,'nv')
        Z <- Z[,nv==max(nv)]
        t <- t[nv==max(nv)]
      }
      #print(dim(Z))
      mu <- apply(Z,2,'mean')
      sigma <- apply(Z,2,'sd')
      q10 <- apply(Z,2,quantile,probs=0.1)
      q90 <- apply(Z,2,quantile,probs=0.9)
      
      #print(paste('plume',sce)); print(dim(Z)); print(length(t))
      if (verbose) {print(col.scen[grep(sce,scenarios)]);print(descr)}
      cols <- c(cols,col.scen[grep(sce,scenarios)])
      plume(Z,t,col=col.scen[grep(sce,scenarios)],levs=levs,verbose=verbose)
      #print('...')
      if (!is.null(highlight)) {
        if (verbose) print(paste('Highlight',paste(highlight,collapse=', ')))
        for (im in 1:length(highlight)) { 
          show <- grep(highlight[im],attr(x[[1]],"model_id"))
          for (ii in show) lines(t,Z[ii,])
        }
      }
      if (!is.null(threshold)) {
        if (verbose) print(paste('Threshold provided: n above',threshold,
                                 'sce=',sce))
        ensqp <- zoo(
          apply(Z,2,function(x) 100*sum(x < threshold,na.rm=TRUE)/sum(is.finite(x))),
          order.by=t)
        attr(ensqp,'unit') <- '"%"'
        attr(ensqp,'description') <- paste('Persentage of ensemble below',
                                           threshold,esd::unit(y))
        attr(ensqp,'variable') <- varid(y)
        attr(ensqp,'ensemble_mean') <- mu
        attr(ensqp,'ensemble_sd') <- sigma
        qp[[sce]] <- ensqp
      }
      
      mu.1991.2020 <- round(mean(mu[is.element(t,1991:2020)]),2)
      mu.2041.2070 <- round(mean(mu[is.element(t,2041:2070)]),2)
      if (!is.finite(mu.2041.2070)) browser("!is.finite(mu.2041.2070)")
      mu.2071.2100 <- round(mean(mu[is.element(t,2071:2100)]),2)
      q10.1991.2020 <- round(mean(q10[is.element(t,1991:2020)]),2)
      q90.1991.2020 <- round(mean(q90[is.element(t,1991:2020)]),2)
      ## Future
      q10.2041.2070 <- round(mean(q10[is.element(t,2041:2070)]),2)
      q10.2071.2100 <- round(mean(q10[is.element(t,2071:2100)]),2)
      q90.2041.2070 <- round(mean(q90[is.element(t,2041:2070)]),2)
      q90.2071.2100 <- round(mean(q90[is.element(t,2071:2100)]),2)
      print(paste0("1991-2020: ",mu.1991.2020,
                   " (",q10.1991.2020," - ",q90.1991.2020,")"))
      print(paste0("2041-2070: ",mu.2041.2070,
                   " (",q10.2041.2070," - ",q90.2041.2070,")"))
      print(paste0("2071-2100: ",mu.2071.2100,
                   " (",q10.2071.2100," - ",q90.2071.2100,")"))
    }
  }
  lines(index(y),coredata(y),showobstype,pch=19,lwd=2)
  
  if (verbose) {print(names(x)); print(sce); print(descr); print(rbind(scenarios,NM))}
  scenarios <- descr
  for (ijk in scenarios) {
    #iscen <- toupper(sub(paste0(param,'.'),'',ijk,fixed=TRUE))
    if (length(grep('ssp',scenarios))>0)
      IJK <- toupper(paste(substr(ijk,1,4),substr(ijk,5,nchar(ijk)),sep='-')) else
        IJK <- toupper(ijk)
    descr <- sub(ijk,paste0(NM[grep(ijk,scenarios)],
                            ' GCM runs (',IJK,')'),descr)
  }
  # descr <- sub('esd',paste('ESD for',NM[is.element(scenarios,'esd')],
  #                          'GCM runs (SSP3-70)'),descr)
  # descr <- sub('ssp370',paste(NM[grep(scenarios,'ssp370')],
  #                             'GCM runs (SSP3-70)'),descr)
  # descr <- sub('ssp245',paste(NM[is.element(scenarios,'ssp245')],
  #                             'GCM runs (SSP2-45)'),descr)
  # descr <- sub('ssp126',paste(NM[is.element(scenarios,'ssp126')],
  #                             'GCM runs (SSP1-26)'),descr)
  # descr <- sub('ssp585',paste(NM[is.element(scenarios,'ssp585')],
  #                             'GCM runs (SSP5-85)'),descr)
  # descr <- sub('rcp45',paste(NM[is.element(scenarios,'rcp45')],
  #                            'GCM runs (RCP45)'),descr)
  # descr <- sub('rcp85',paste(NM[is.element(scenarios,'rcp85')],
  #                            'GCM runs (RCP85)'),descr)
  if (!is.null(highlight)) {
    descr <- c(descr,highlight)
    cols <- c(cols,rep('black',length(highlight)))
  }
  if (verbose) {print('legend'); print(cols)}
  legend('topleft',c('Observed',descr),
         col=c('black',cols),pch=15,cex=0.75,bty='n')
  invisible(qp)
}

#' @title tableofmodels
#' @export tableofmodels
tableofmodels <- function(x) {
  ensembles <- names(x)
  ensembles <- ensembles[grep('.DJF',ensembles)]
  for (ensemble in ensembles) { 
    print(sub('.DJF','',ensemble,fix=TRUE))
    models <- sub("\\.r.*$","",attr(x[[ensemble]],"model_id"))
    print(paste(length(models),'different simulations'))
    print(table(models))
  }
}