### UTILITIES FUNCTIONS ###



# TRIGGERS
# create new trigger with: myTrigger <- makeReactiveTrigger()
# put myTrigger$trigger() in code that should trigger other code
# put myTrigger$depend() in code to run when myTrigger has been triggered
# CANNOT BE USED AS AN EVENT in observeEvent(myTrigger$depend(),{...})
# USE observe({myTrigger$depend() ...}) instead

#####################################################################
### FUNCTION TO MANIPULATE REACTIVE TRIGGERS
#####################################################################
makeReactiveTrigger <- function() {
  rv <- reactiveValues(a = 0)
  list(
    depend = function() {
      rv$a
      invisible()
    },
    trigger = function() {
      rv$a <- isolate(rv$a + 1)
    }
  )
}


#####################################################################
### quantification
#####################################################################
quantif<-function(data,ranges=NULL,nbins=10,rounding=FALSE)
{
  # use round=TRUE to have integer bins (assume (but not checked) x is a set of integers)
  
  if (is.null(ranges)) ranges=range(data)
  
  mn<-ranges[1]
  mx<-ranges[2]

  delta<-(mx-mn)/nbins
  
  
  if (delta>0)
  {
    b<-seq(mn,mx,delta)
    lb<-length(b)
    
    if (rounding==TRUE) cgb<-round(0.5*(b[1:(lb-1)]+b[2:lb]))
    else cgb<-0.5*(b[1:(lb-1)]+b[2:lb])
    
    res<-cgb[.bincode(data,b,FALSE,TRUE)]
  }else{
    res<-mn
  }
  return(res)
  
}

#####################################################################
### CLEAN FORMAT OF LIST GENES TO TRACK
#####################################################################
cleanStrGenesToTrack<-function(strGenesToTrack)
{
strTmp<-strGenesToTrack
strTmp1=gsub(";",",",strTmp)
strTmp2=gsub(" ",",",strTmp1)
strTmp3=gsub(",+",",",strTmp2)

return(sort(strsplit(strTmp3,",")[[1]]))
}

#####################################################################
### CONCATENATE VECTORS OF STRINGS
#####################################################################
implode<-function(strvec,sep=",")
{
  return(paste(strvec,sep="",collapse=sep))
}


