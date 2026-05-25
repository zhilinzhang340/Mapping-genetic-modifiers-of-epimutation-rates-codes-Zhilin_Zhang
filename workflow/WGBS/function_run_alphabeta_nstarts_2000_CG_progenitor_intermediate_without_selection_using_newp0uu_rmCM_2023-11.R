
run.alphabeta.new<-function(nodelist,edelist,name,input.dir,output.dir){
    #rm(list=ls())
    library("AlphaBeta")
    source("/mnt/nas/zhilin/RIL/run_AB/rscripts/newbuildPedigree_stand_alone_rmCM_2023_11.R")
    
    plotPedigree(nodelist =nodelist,
                 edgelist =edelist,
                 sampling.design = "progenitor.intermediate", 
                 output.dir = output.dir, 
                 plot.width = 8, plot.height = 8, aspect.ratio = 0.5, 
                 vertex.size = 4, vertex.label = T, out.pdf = paste0("pedigree_",name, sep=""))
    
    ##get full path of the input files
    path_node<-list.files(input.dir, pattern=nodelist, full.names=TRUE)
    path_ede<-list.files(input.dir, pattern=edelist, full.names=TRUE)
    
    ##read files
    node<-read.table(path_node,header = T)
    ede<-read.table(path_ede,header = T)
    
    ##run alphabeta
    output <- newbuildPedigree(path_node, path_ede, cytosine = "CG", posteriorMaxFilter = 0.99)
    p0uu_in <- output$tmpp0;p0uu_in
    write.table(p0uu_in, paste0("p0uu_in_",name,".txt"),row.names=F,quote = F,col.names = F)
    
    ##plot divergence vs delta_t 
    pedigree <- output$Pdata
    write.table(pedigree,paste0("pedigree-pdata_",name,".txt"),row.names=F,quote = F) 
    
    output_ABneutral<- ABneutral(pedigree.data = pedigree, p0uu = p0uu_in,eqp = p0uu_in, 
                                 eqp.weight = 1, Nstarts = 2000, out.dir = output.dir,
                                 out.name = paste0("ABneutral_CG_estimates_",name,sep=""))
    #summary(output_ABneutral)
    ABplot(pedigree.names = paste0("ABneutral_CG_estimates_",name,".Rdata",sep=""), output.dir = output.dir, 
           out.name = paste0("ABneutral_CG_",name,sep=""), plot.height = 8,plot.width = 11)
    
    ABneutral_estimatats<-output_ABneutral$estimates[1:10,1:5]
    library(dplyr)
    ABneutral_estimatats<-tibble::rownames_to_column(ABneutral_estimatats,name)
    write.table(ABneutral_estimatats,paste0("ABneutral_estimatats_",name,".txt"),quote = F,row.names = F)
    
    pdf(paste0(output.dir,"ABneutral_pedigree_",name,".pdf", sep=""),width = 11,height =8)
    plot(output_ABneutral$pedigree[,"delta.t"],  output_ABneutral$pedigree[,"div.obs"])
    points( output_ABneutral$for.fit.plot[,"delta.t"], output_ABneutral$for.fit.plot[,"div.sim"], type="l")
    dev.off()
    
    #run ABnull
    output_ABsUU <- ABnull(pedigree.data = pedigree,out.dir =  output.dir, out.name =  paste0("ABnull_CG_estimates_",name,sep=""))
    summary(output_ABsUU)
    
    
    ##compare two models  ###use "<<-"" for define the global variables
    
    out_ABneutral_ABnull<- FtestRSS(pedigree.select =  paste0("ABneutral_CG_estimates_",name,".Rdata",sep=""), pedigree.null =  paste0("ABnull_CG_estimates_",name,".Rdata",sep=""))
    out_ABneutral_ABnull_Ftest_name<-paste0("out_ABneutral_ABnull_Ftest_",name, sep="")
    out_ABneutral_ABnull_cb<-cbind(out_ABneutral_ABnull_Ftest_name,out_ABneutral_ABnull$Ftest[6])
    write.table(out_ABneutral_ABnull_cb,"out_ABneutral_ABnull_Ftest_p_value.txt",append = T,row.names=F,col.names = F,quote = F)
    # 
    # # ##run ABselectMM
    # output_ABselectMM<- ABselectMM(pedigree.data = pedigree, p0uu = p0uu_in, eqp = p0uu_in,
    #                                eqp.weight = 1, Nstarts = 2000, out.dir = output.dir, out.name = paste0("ABselectMM_CG_estimates_",name,sep=""))
    # #summary(output_ABselectMM)
    # ABplot(pedigree.names = paste0("ABselectMM_CG_estimates_",name,".Rdata",sep=""), output.dir = output.dir,
    #        out.name = paste0("ABselectMM_CG_",name,sep=""), plot.height = 8, plot.width = 11)
    # 
    # ##run ABselectUU
    # output_ABselectUU <- ABselectUU(pedigree.data = pedigree, p0uu = p0uu_in, eqp = p0uu_in,
    #                                 eqp.weight = 1, Nstarts = 2000, out.dir =  output.dir, out.name =paste0("ABselectUU_CG_estimates_",name,sep=""))
    # #summary(output_ABselectUU)
    # ABplot(pedigree.names = paste0("ABselectUU_CG_estimates_",name,".Rdata",sep=""), output.dir = output.dir,
    #        out.name = paste0("ABselectUU_CG_",name,sep=""), plot.height = 8,plot.width = 11)
    # 
    # 
    # 
    # out_ABselectMM_ABneutral<- FtestRSS(pedigree.select =  paste0("ABselectMM_CG_estimates_",name,".Rdata",sep=""), pedigree.null =  paste0("ABneutral_CG_estimates_",name,".Rdata",sep=""))
    # out_ABselectMM_ABneutral_Ftest_name<-paste0("out_ABselectMM_ABneutral_Ftest_",name, sep="")
    # out_ABselectMM_ABneutral_cb<-cbind(out_ABselectMM_ABneutral_Ftest_name,out_ABselectMM_ABneutral$Ftest[6])
    # write.table(out_ABselectMM_ABneutral_cb,"out_ABselectMM_ABneutral_Ftest_p_value.txt",append = T,row.names=F,col.names = F,quote = F)
    # #assign(out_ABselectMM_ABneutral_Ftest_name,out_ABselectMM_ABneutral$Ftest,envir = .GlobalEnv)
    # #out_ABselectMM_ABneutral_Ftest_name<<-out_ABselectMM_ABneutral$Ftest
    # #return(out_ABselectMM_ABneutral_name)
    # # 
    # # 
    # out_ABselectUU_ABneutral<- FtestRSS(pedigree.select =  paste0("ABselectUU_CG_estimates_",name,".Rdata",sep=""), pedigree.null = paste0("ABneutral_CG_estimates_",name,".Rdata",sep=""))
    # out_ABselectUU_ABneutral_Ftest_name<-paste0("out_ABselectUU_ABneutral_Ftest_",name, sep="")
    # out_ABselectUU_ABneutral_cb<-cbind(out_ABselectUU_ABneutral_Ftest_name,out_ABselectUU_ABneutral$Ftest[6])
    # write.table(out_ABselectUU_ABneutral_cb,"out_ABselectUU_ABneutral_Ftest_p_value.txt",append = T,row.names=F,col.names = F,quote = F)
    # #assign(out_ABselectUU_ABneutral_Ftest_name,out_ABselectUU_ABneutral$Ftest,envir = .GlobalEnv)
    # #out_ABselectUU_ABneutral_Ftest_name<<-out_ABselectUU_ABneutral$Ftest
    # #return(out_ABselectUU_ABneutral_name)
    # 
    # 
    ##run BOOTmodel Nboot = 1000
    
    Boutput<- BOOTmodel(pedigree.data = paste0("ABneutral_CG_estimates_",name,".Rdata",sep=""), Nboot = 1000, out.dir = output.dir,
                        out.name = paste0("ABneutral_Boot_CG_estimates_",name,sep=""))
    #summary(Boutput)
    Boutput_boot_base_name<-paste0("Boutput_boot_base_",name, sep="")
    #assign(Boutput_boot_base_name,Boutput$boot.base,envir = .GlobalEnv)
    Boutput_boot_base_alpha_name<-paste0("Boutput_boot_base_alpha_",name, sep="")
    Boutput_boot_base_beta_name<-paste0("Boutput_boot_base_beta_",name, sep="")
    Boutput_boot_base_alpha_cb<-cbind(Boutput_boot_base_alpha_name,Boutput$boot.base$alpha)
    Boutput_boot_base_beta_cb<-cbind(Boutput_boot_base_beta_name,Boutput$boot.base$beta)
    write.table(Boutput_boot_base_alpha_cb,"Boutput_boot_base_alpha.txt",append = T,row.names=F,col.names = F,quote = F)
    write.table(Boutput_boot_base_beta_cb,"Boutput_boot_base_beta.txt",append = T,row.names=F,col.names = F,quote = F)
    
    #Boutput_boot_base_name<<-Boutput$boot.base
    
    Boutput_standard_errors_name<-paste0("Boutput_standard_errors_",name, sep="")
    #assign(Boutput_standard_errors_name,Boutput$standard.errors,envir = .GlobalEnv)
    Boutput_standard_errors_alpha_name<-paste0("Boutput_standard_errors_alpha_",name, sep="")
    Boutput_standard_errors_beta_name<-paste0("Boutput_standard_errors_beta_",name, sep="")
    Boutput_standard_errors_beta_alpha_name<-paste0("Boutput_standard_errors_beta_alpha_",name, sep="")
    #assign: use var name var ; use envir = .GlobalEnv set global environment
    
    Boutput_standard_errors_alpha_cb<-cbind(Boutput_standard_errors_alpha_name,Boutput$standard.errors[1,1])
    Boutput_standard_errors_beta_cb<-cbind(Boutput_standard_errors_beta_name,Boutput$standard.errors[2,1])
    Boutput_standard_errors_beta_alpha_cb<-cbind(Boutput_standard_errors_beta_alpha_name,Boutput$standard.errors[3,1])
    write.table(Boutput_standard_errors_alpha_cb,"Boutput_standard_errors_alpha.txt",append = T,row.names=F,col.names = F,quote = F)
    write.table(Boutput_standard_errors_beta_cb,"Boutput_standard_errors_beta.txt",append = T,row.names=F,col.names = F,quote = F)
    write.table(Boutput_standard_errors_beta_alpha_cb,"Boutput_standard_errors_beta_alpha.txt",append = T,row.names=F,col.names = F,quote = F)
}

