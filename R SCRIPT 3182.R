install.packages("quantmod")
install.packages("tidyquant")
install.packages("imputeTS")
install.packages("dplyr")
library("quantmod")
library("tidyquant")
library("imputeTS")
library("dplyr")


Tickers=c("^GSPC", "^DJI","^GSPE","^SP500-35","^SP500-45","^SP500-60",
          "SHY","^SPGSCI" )



getSymbols(Tickers, from ='2002-07-30', to ='2022-11-02')

OG_GSPC<-c(GSPC$GSPC.Adjusted)
OG_DJI<-c(DJI$DJI.Adjusted)
OG_GSPE<-c(GSPE$GSPE.Adjusted)
OG_SP35<-c(`SP500-35`$`SP500-35.Adjusted`)
OG_SP45<-c(`SP500-45`$`SP500-45.Adjusted`)
OG_SP60<-c(`SP500-60`$`SP500-60.Adjusted`)
OG_SHY<-c(SHY$SHY.Adjusted)
OG_SPGSCI<-c(SPGSCI$SPGSCI.Adjusted)


###############

#For log returns

Ret_Log<- function(x){

seq<-c()  
    
  for (i in 2:length(x)){
    
    d<- log(as.numeric(x[i])/as.numeric(x[i-1]))
    seq<-c(seq,d)
    
    seq[is.na(seq)]<-1
    
  }
  
  return(seq)
  
}

LOG_GSPC<-Ret_Log(OG_GSPC)
LOG_DJI<-Ret_Log(OG_DJI)
LOG_GSPE<-Ret_Log(OG_GSPE)
LOG_SP35<-Ret_Log(OG_SP35)
LOG_SP45<-Ret_Log(OG_SP45)
LOG_SP60<-Ret_Log(OG_SP60)
LOG_SHY<-Ret_Log(OG_SHY)
LOG_SPGSCI<-Ret_Log(OG_SPGSCI)


