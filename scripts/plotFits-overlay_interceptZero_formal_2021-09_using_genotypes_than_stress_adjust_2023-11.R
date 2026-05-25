test_plotFITS <- function(pedigree.names, input.dir, output.dir, out.name, 
                          alpha, geom.point.size, geom.line.size, plot.height, plot.width, plot.type, lsq.line, facets.ncols) {
  
  #initialize vectors
  df_div <- c()
  pred.df_fits <- c()
  theory.df_fits <- c()
  gg <- ggplot()
  all.miss<-c()
  #NCOL and NROW do the same, treating a vector as 1-column matrix.
  for(i in 1:NROW(pedigree.names)) {
    for(j in 1:NCOL(pedigree.names)) {
      
      datain <- dget(paste(input.data.dir, pedigree.names[i,j], sep=""))
      
      name <- gsub("(_).*_ABneutral_CG_estimates.Rdata$", "", basename(as.character(pedigree.names[i,j])))
      print(name)
      Genotypes <- gsub(pattern=paste0(name, "_|\\_ABneutral_CG_estimates.Rdata$"), "", 
                     basename(as.character(pedigree.names[i,j])))
      print(Genotypes)
      
      intercept <- datain$estimates[1, "intercept"]
      delta.t <- datain$pedigree[,"delta.t"]
      div.obs <- datain$pedigree[,"div.obs"]-intercept
      n.df <- data.frame(delta.t, div.obs)
      name.div <- cbind(Genotypes, name, n.df)
      df_div <- rbind(df_div, name.div)
      miss<-cbind(paste0(name,"_",Genotypes),length(div.obs[div.obs<0]))
      all.miss<-rbind(all.miss,miss)
      
      #predictive fit
      pred.fit.data <- aggregate(datain$pedigree[,"div.pred"], by=list(datain$pedigree[,"delta.t"]), median)
      pred.fits <- c(intercept, pred.fit.data[,2])
      pred.fits <- pred.fits-intercept
      pred.fit.t <- c(0, pred.fit.data[,1])
      
      pred.df_fits <- rbind(pred.df_fits, data.frame(Genotypes, name, pred.fit.t, pred.fits))
      
      #theoritical fit
      theory.fit.data <- datain$for.fit.plot
      theory.fits <- c(intercept, theory.fit.data[,"div.sim"])
      theory.fits <- theory.fits-intercept
      theory.fit.t <- c(0, theory.fit.data[,"delta.t"])
      theory.df_fits <- rbind(theory.df_fits, data.frame(Genotypes, name, theory.fit.t, theory.fits))
    }
  }
  write.table(all.miss,"output_miss_point_number_below_zero_during_divergence_plot.txt",append = T,row.names=F,col.names = F,quote = F)
  
  if ((plot.type == "data.only")) {
    gg <- gg + geom_point(data=df_div, aes(x=delta.t, y=div.obs, colour=name), alpha = alpha, size=geom.point.size)
  }
  
  if ((plot.type == "fit.only") && (lsq.line == "pred")) {
    gg <- gg + geom_line(data=pred.df_fits, aes(x=pred.df_fits$pred.fit.t, y=pred.df_fits$pred.fits, color=pred.df_fits$name), size=geom.line.size)
  }
  
  if ((plot.type == "fit.only") && (lsq.line == "theory")) {
    gg <- gg + geom_line(data=theory.df_fits, aes(x=theory.fit.t, y=theory.fits, color=name), size=geom.line.size)
  }
  
  if ((plot.type == "both") && (lsq.line == "pred")) {
    gg <- gg + geom_point(data=df_div, aes(x=delta.t, y=div.obs, colour=name), alpha = alpha, size=geom.point.size) + 
      geom_line(data=pred.df_fits, aes(x=pred.fit.t, y=pred.fits, color=name), size=geom.line.size)
  }
  
  if ((plot.type == "both") && (lsq.line == "theory")) {
    gg <- gg + geom_point(data=df_div, aes(x=delta.t, y=div.obs, colour=name), alpha = alpha, size=geom.point.size) + 
      geom_line(data=theory.df_fits, aes(x=theory.fit.t, y=theory.fits, color=name), size=geom.line.size)
  }
  
  
  #pdf(paste(output.dir, out.name, ".pdf", sep=""), colormodel = 'cmyk', width = plot.width, height = plot.height)
  
  gg <- gg + labs(x="delta t (generations)", y="mCG divergence", color="Genotypes") 
  
  ymax <- 1.1*max(df_div$div.obs)
  xmax <- 1.1*max(df_div$delta.t)
  
  Genotypes_1<-str_split_fixed(Genotypes,"_",4)
  Genotypes_2<-paste0(Genotypes_1[1,1]," ",Genotypes_1[1,2]," ",Genotypes_1[1,3]," ",Genotypes_1[1,4])
  
  pdf(paste(output.dir, out.name, ".pdf", sep=""),  width = plot.width, height = plot.height)
  #colormodel = 'cmyk'
  gg <- gg + labs(x=expression(paste(Delta,"t (generations)")), y="mCG divergence", color="Genotypes",title=Genotypes_2)
  
  gg <- gg +
    theme_classic()+
    theme(text = element_text(size = 15),plot.title=element_text(hjust=0.5),axis.text.x=element_text())+
    scale_x_continuous(expand = c(0, 0), limits=c(0, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) 
  #facet_wrap(~Genotypes, ncol=facets.ncols) 
  print(gg)
  dev.off()
}
