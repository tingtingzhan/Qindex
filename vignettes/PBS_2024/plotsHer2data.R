 
 
##################### ##################### ##################### ##################### ##################### 
############# plots for PSB talk
##################### ##################### ##################### ##################### ##################### 

plots_PSB_2024<- "/Users/ixc107/Documents/ICH_projects/NAHormTxR03data/plots/FigsPSB_2024"
   setwd(plots_PSB_2024)

require(spatstat)
require(spatstat.geom)
 
 load(file="/Users/ixc107/Documents/ICH_projects/NAHormTxR03data/Rdata/trimmed_all/Surger_case36R.RData")
Surger_case=dataout
 head(Surger_case)
 dim(Surger_case)
 
################# remove duplicated points, IF ANY

Surger_case$coord=paste0(Surger_case$Cell.X.Position, "_", Surger_case$Cell.Y.Position )
Surger_case$duplicated=as.numeric(duplicated(Surger_case$coord))
head(Surger_case)
table(Surger_case$duplicated)

data=Surger_case[(Surger_case$duplicated<0.5)&(((Surger_case$Phenotype.ER=="ER+")|(Surger_case$Phenotype.ER=="Others"))),]
head(data)

xrange=c(min(data$Cell.X.Position),max(data$Cell.X.Position))
  yrange=c(min(data$Cell.Y.Position),max(data$Cell.Y.Position)) 
  
  table(data$Phenotype.Ki67)
data$Ki67pos=as.numeric(data$Phenotype.Ki67=="Ki67+")

data$Ki67pos[data$Ki67pos>0.5]="Ki67+"
data$Ki67pos[data$Ki67pos<0.5]="Ki67-"

  table(data$Phenotype.Ki67, data$Ki67pos)
     table(data$Tissue.Category, data$Ki67pos)
     
table(data$Phenotype.ER)     
data$ERpos=as.numeric(data$Phenotype.ER=="ER+")     
data$ERpos[data$ERpos>0.5]="ER+"
data$ERpos[data$ERpos<0.5]="ER-"
  table(data$Phenotype.ER, data$ERpos)
  
boxplot(Nucleus.ER..Opal~Phenotype.ER, data=data)
  
###### ###### ###### MMPP of Ki67 hsps  
   Ki67cont.ppp=ppp(data$Cell.X.Position, data$Cell.Y.Position, xrange, yrange, marks=data$Nucleus.Ki67..Op)  
   plot(Ki67cont.ppp) 

################# trim the data into smaller window
ptrimR=0.55
ptrimL=0.4
  ptrimB=0.40
    ptrimT=0.5
    
 horis_T=(1-ptrimT)*max(data$Cell.Y.Position)+ptrimT*min(data$Cell.Y.Position)    
 horis_B=(1-ptrimB)*min(data$Cell.Y.Position)+ptrimB*max(data$Cell.Y.Position)
 
 vert_L=(1-ptrimL)*min(data$Cell.X.Position)+ptrimL*max(data$Cell.X.Position) 
 vert_R=(1-ptrimR)*max(data$Cell.X.Position)+ptrimR*min(data$Cell.X.Position)
 
 data1=data[(vert_L<data$Cell.X.Position)&(data$Cell.X.Position<vert_R)&(horis_B<data$Cell.Y.Position)&(data$Cell.Y.Position<horis_T), ]  
   dim(data1)
   xrange1=c(min(data1$Cell.X.Position),max(data1$Cell.X.Position))
  yrange1=c(min(data1$Cell.Y.Position),max(data1$Cell.Y.Position))
  
################# ################# ################# ################# ER
  summary(data1$Nucleus.ER..Opal)
ERcont1.ppp=ppp(data1$Cell.X.Position, data1$Cell.Y.Position, xrange1, yrange1, marks=0.001+data1$Nucleus.ER..Opal)  
data1$ERpos=relevel(as.factor(data1$ERpos), ref="ER+")     
ERpos.ppp=ppp(data1$Cell.X.Position, data1$Cell.Y.Position, xrange1, yrange1, marks=factor(data1$ERpos))      

   pdf(file = paste0("ERboth.pdf"),  width = 8, height = 4, pointsize = 14)
  par(mar=c(0,0,0,2), mgp=c(0,0.5,0), mfrow=c(1,2))
     plot(ERpos.ppp, pch=c(19,20), cols=c("red","grey"), cex=0.5, leg.side="right", main=" ")
      plot(ERcont1.ppp, leg.side="right",  main=" ")
  dev.off()    
   
   
################# ################# ################# ################# Ki67

Ki67cont1.ppp=ppp(data1$Cell.X.Position, data1$Cell.Y.Position, xrange1, yrange1, marks=0.001+data1$Nucleus.Ki67..Op)  
data$Ki67pos=relevel(as.factor(data$Ki67pos), ref="Ki67+")     
Ki67pos.ppp=ppp(data1$Cell.X.Position, data1$Cell.Y.Position, xrange1, yrange1, marks=factor(data1$Ki67pos)) 

#par(mfrow=c(1,2))

 pdf(file = paste0("Ki67posMMPP.pdf"),  width = 5, height = 5, pointsize = 12)
  par(mar=c(0,0,0,2), mgp=c(0,0.5,0), mfrow=c(1,1))
   plot(Ki67pos.ppp, pch=c(19,20), cols=c("red","grey"), cex=0.5, leg.side="right", main=" ")
  dev.off() 
  
   pdf(file = paste0("Ki67contMPP.pdf"),  width = 5, height = 5, pointsize = 12)
  par(mar=c(0,0,0,2), mgp=c(0,0.5,0), mfrow=c(1,1))
      plot(Ki67cont1.ppp, leg.side="right", markscale=2, main=" ")
  dev.off() 
  
  
   pdf(file = paste0("Ki67both.pdf"),  width = 8, height = 4, pointsize = 14)
  par(mar=c(0,0,0,2), mgp=c(0,0.5,0), mfrow=c(1,2))
     plot(Ki67pos.ppp, pch=c(19,20), cols=c("red","grey"), cex=0.5, leg.side="right", main=" ")
      plot(Ki67cont1.ppp, leg.side="right", markscale=2, main=" ")
  dev.off() 




  
    
    