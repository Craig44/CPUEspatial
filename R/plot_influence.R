#' plot_influence
#'
#' @details plotting function to generate CDI plots from The Bentley paper
#'
#' @param influence_info a list that has been outputed by calculate_influence() function
#' @param term character that represents the covariate term you want the plot to generate it for.
#' @return gg plot
#' @export
#'
#'
plot_influence = function(influence_info, term) {
  par(oma=c(1,1,1,1),cex.lab=1.25,cex.axis=1.25)
  layout(matrix(c(1,1,0,2,2,3,2,2,3),3,3,byrow=TRUE))
  # PLot coeffecients
  par(mar=c(0,5,3,0),las=1)
  with(influence_info$terms_ls[[term]],{
    xs = 1:max(as.integer(labs))
    ylim = c(min((lower)),max((upper)))
    if(ylim[1]<0.5*min((MLE))) ylim[1] = 0.5*min((MLE))
    if(ylim[2]>2*max((MLE))) ylim[2] = 2*max((MLE))
    plot(as.integer(labs),(MLE),xlim=range(xs),ylim=ylim,pch=16,cex=1.5,xaxt='n',ylab='', lwd = 2)
    mtext('Coefficient (+- 2SE)',side=2,line=4,las=0,cex=0.8)
    abline(h=1,lty=2)
    abline(v=xs,lty=1,col='grey')
    segments(as.integer(labs),(lower),as.integer(labs),(upper), lwd = 2)
    arrows(as.integer(labs),(lower),as.integer(labs),(upper),angle=90,code=3,length=0.05, lwd = 2)
    axis(3,at=xs,labels=levels(labs)[xs])
  })
  ## distribution of events (observatons) not catch
  par(mar=c(5,5,0,0),las=1)
  xlab = term
  ylab= non_spatial$Call$func_call$time_variable_label
  with(influence_info$all_distr[[term]],{
    xs = 1:max(as.integer(term))
    ys = 1:max(as.integer(get(ylab)))
    plot(NA,xlim=range(xs),ylim=range(ys),xaxt='n',yaxt='n',xlab=xlab,ylab='')
    abline(v=xs,lty=1,col='grey')
    axis(1,at=xs,labels=levels(term)[xs])
    abline(h=ys,lty=1,col='grey')
    axis(2,at=ys,labels=levels(get(ylab))[ys])
    mtext(ylab, side = 2, line=4, las=0, cex=1)
    points(as.integer(term),as.integer(get(ylab)),cex=sqrt(n_prop)*12)
  })
  ## plot influence plots
  #Influence
  par(mar=c(5,0,0,3),las=1)
  ys = 1:nrow(.$influences)
  with(influence_info$influence_df,{
    plot(NA,xlim=range(exp(get(term))),ylim=range(ys),yaxt='n',xlab='Influence')
    abline(h = ys,lty = 1,col='grey')
    abline(v = 1 ,lty = 2)
    points(exp(get(term)),ys,type='o',pch=16,cex=1.8)
    axis(4, at = ys, labels = time)
  })
}
