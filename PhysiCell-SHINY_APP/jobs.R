source("app.R")
LoadSim("/stock/Physicell_data_backup/22-07-19bis/")

#pheno3<-plotCellPhenotypeAcrossTime(seq(0,300,1),"/stock/Physicell_data_backup/22-07-19bis/","PHENOTYPE")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))

#warburg1<-plotCellPhenotypeAcrossTime(seq(0,300,1),"/stock/Physicell_data_backup/22-07-18/","PHENOTYPE2")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))

#warburg2<-plotCellPhenotypeAcrossTime(seq(0,300,1),"/stock/Physicell_data_backup/22-07-19/","PHENOTYPE2")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))

#warburg3<-plotCellPhenotypeAcrossTime(seq(0,300,1),"/stock/Physicell_data_backup/22-07-19bis/","PHENOTYPE2")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))

#cellpopref<-plotCellPopAcrossTime(seq(0,430,10),"/stock/Physicell_data_backup/22-02-14/")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))

#pHpheno1<-plotCellPhenotypeAcrossTime(seq(0,13,1),"/stock/Physicell_data_backup/22-07-20/","PHENOTYPE")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))
#pHpheno2<-plotCellPhenotypeAcrossTime(seq(0,13,1),"/stock/Physicell_data_backup/22-07-21/","PHENOTYPE")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))
#pHpheno3<-plotCellPhenotypeAcrossTime(seq(0,13,1),"/stock/Physicell_data_backup/22-07-21bis/","PHENOTYPE")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))
#pHpheno4<-plotCellPhenotypeAcrossTime(seq(0,13,1),"/stock/Physicell_data_backup/22-07-22/","PHENOTYPE")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))
#pHwarburg1<-plotCellPhenotypeAcrossTime(seq(0,13,1),"/stock/Physicell_data_backup/22-07-20/","PHENOTYPE2")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))
#pHwarburg2<-plotCellPhenotypeAcrossTime(seq(0,13,1),"/stock/Physicell_data_backup/22-07-21/","PHENOTYPE2")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))
#pHwarburg3<-plotCellPhenotypeAcrossTime(seq(0,13,1),"/stock/Physicell_data_backup/22-07-21bis/","PHENOTYPE2")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))
#pHwarburg4<-plotCellPhenotypeAcrossTime(seq(0,13,1),"/stock/Physicell_data_backup/22-07-22/","PHENOTYPE2")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))

phenoref<-plotCellPhenotypeAcrossTime(seq(0,300,5),"/stock/Physicell_data_backup/22-02-14/","PHENOTYPE")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))

glucosepheno1<-plotCellPhenotypeAcrossTime(seq(0,300,5),"/stock/Physicell_data_backup/22-07-26/","PHENOTYPE")%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)))
