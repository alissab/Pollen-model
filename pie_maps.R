# plots pie charts for taxa composition on a map; note this is specialized to the STEPPS1 New England domain in the input args defaults, but should work for a different domain
# proportions should be a n x p matrix with columns being different taxa
# centers should be a n x 2 matrix of locations
pieMap=function(proportions, 
                centers, 
                restrict = FALSE,
                inputRestricted = FALSE,
                xlim = c(-52000,940000),
                ylim = c(676000, 1484000),
                radius = NULL,
                scale = 1,
                xlab = 'x',
                ylab = 'y', 
                add_legend, 
                main_title,
                col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
                             "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"),
                cont_shp = cont_shp){
  # plots multiple pie composition charts as a map
  if(!restrict){
    used=1:nrow(centers)
  } else{ # assumes subset of grid
    if(inputRestricted){
      used=1:144
    } else{
      used=49:(256-64)
    }
  }
  centers=as.matrix(centers[used,])
  proportions=as.matrix(proportions[used,])
  if(is.null(xlim)){
    rg=range(centers[,1])
    df=(scale-1)*diff(range(centers[,1]))
    xlim=c(rg[1]-df,rg[2]+df)
  }
  if(is.null(ylim)){
    rg=range(centers[,2])
    df=(scale-1)*diff(range(centers[,2]))
    ylim=c(rg[1]-df,rg[2]+df)
  }
  plot(centers,type='n',xlim=xlim,ylim=ylim,xaxt='n', yaxt='n',ann=FALSE, frame.plot=F, asp=1, main="Pie plots")#, xlab=xlab,ylab=ylab)
  
  # read in North America shape files
  # na_shp <- readOGR("NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
  # na_shp <- sp::spTransform(na_shp, proj_out)
  # cont_shp <- subset(na_shp,
  #                    (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
  # us.shp <- readShapeLines('NA_States_Provinces_Alberts.shp',
  #                          proj4string=CRS('+init=epsg:3175'))
  plot(cont_shp, add=T, lwd=2)
  n=length(centers[,1])
  
  #   col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
  #                "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")
  cols     = col_list[1:ncol(proportions)]
  #   if (ncol(proportions)<13){
  #     cols     = col_list[1:ncol(proportions)]
  #   } else {
  #     cols = rainbow(ncol(proportions))
  #   }
  #if (ncol(proportions) > 12){print('More than 12 categories, aggregate to distinguish!')}
  
  if(is.null(radius)){
    radius=.025*diff(range(centers[,1]))
  }
  minVal=min(proportions)
  # eps=1e-10*max(proportions)
  # if(minVal==0){
  #   warning("Some proportions are zero; proceeding with jittering.")
  #   proportions=proportions+eps
  # }
  if(minVal<0){
    stop("Some proportions are less than zero.")
  }
  if(length(radius)==1){ radius=rep(radius,n)}
  for(i in 1:n){
    if(sum(proportions[i,])>0){
      minVal=min(proportions[i,])
      if(minVal==0){
        warning("Some proportions are zero; proceeding with jittering.")
        eps=1e-10*max(proportions[i,])
        proportions[i,]=proportions[i,]+eps
      }
      pieAdd(as.vector(proportions[i,]),as.vector(centers[i,]),radius=radius[i], col=cols)
    } else{
      #points(centers[i,],pch='X')
    }
  }
  #   par(xpd=NA) 
  #   tmp <- cnvrt.coords(.9,.7, 'tdev')$usr 
  if (add_legend){
    legend.col=c(0,1,1,1,0,1,1,1) 
    #     900000,1500000
    legend('topright', colnames(proportions), pch=rep(22,n), pt.cex=1.6, cex=1.2, pt.bg=cols, col=rep('black', n),
           bg='white', ncol=2)
    title(main=main_title, cex.main=2)
  }
}

# auxiliary fxn needed by pieMap
pieAdd= function (x, center, labels = names(x), edges = 200, radius = 0.8, density = NULL, 
                  angle = 45, col = NULL, border = NULL, lty = NULL,
                  col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
                               "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")) # modified from the pie() function in R
{
  if (!is.numeric(x) || any(is.na(x) | x <= 0)) 
    stop("pie: `x' values must be positive.")
  if (is.null(labels)) 
    labels <- rep("",length(x))
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  
  pin <- par("pin")
  nx <- length(dx)
  if (is.null(col)){ 
    col <- if (is.null(density)) 
      #     col_list = c("#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A", "#B15928", "#FF7F00", 
      #                  "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")
      cols     = col_list[1:nx]
    #     if (ncol(proportions)<13){
    #       cols     = col_list[1:nx]
    #     } else {
    #       cols = rainbow(ncol(proportions))
    #     }
    
    #topo.colors(nx)
    #       cols = brewer.pal(12, 'Paired')
    #       cols = c("white", "black","lightblue", "darkblue", "red","yellow",
    #        "purple","orange","lightgreen","darkgreen", "burlywood", 
    #         "hotpink", "indianred4")
    #       cols = c("#1F78B4", "#33A02C", "#FF7F00", "#E31A1C", "#6A3D9A", "#B15928", 
    #                "#FFFF99", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6")
    #       cols[1:nx]
    #       if (nx > 12){print('More than 12 categories, aggregate to distinguish!')}
  }else par("fg")
  col <- rep(col, length.out = nx)
  border <- rep(border, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(angle, length.out = nx)
  density <- rep(density, length.out = nx)
  for (i in 1:nx) {
    n <- max(2, floor(edges * dx[i]))
    t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
    xc <- c(cos(t2p), 0) * radius + center[1]
    yc <- c(sin(t2p), 0) * radius + center[2]
    polygon(xc, yc, density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    t2p <- 2 * pi * mean(x[i + 0:1])
    xc <- cos(t2p) * radius + center[1]
    yc <- sin(t2p) * radius + center[2]
    if (!is.na(lab <- labels[i]) && lab != "") {
      lines(c(1, 1.05) * xc, c(1, 1.05) * yc)
      text(1.1 * xc, 1.1 * yc, lab, xpd = TRUE, adj = ifelse(xc < 
                                                               0, 1, 0))
    }
  }
  invisible(NULL)
}

