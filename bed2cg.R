bed2cg<-function(bed1){
  ref<-read.table("~/oasis/db/hg19/GPL13534.sort.bed",head=F,sep="\t")
  rbedintersect<-function(bed1,ref){
    Rbedtools<-function(functionstring="intersectBed",bed1,bed2,opt.string=""){
      #create temp files
      a.file=tempfile()
      b.file=tempfile()
      out   =tempfile()
      options(scipen =99) # not to use scientific notation when writing out
      #write bed formatted dataframes to tempfile
      write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
      write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
      # create the command string and call the command using system()
      command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
      cat(command,"\n")
      try(system(command))
      res=read.table(out,header=F)
      unlink(a.file);unlink(b.file);unlink(out)
      res=subset(res,V5!=".")
      return(res)
    }
    merge<-Rbedtools(functionstring="intersectBed",bed1,ref,opt.string="-wao")
    return(merge)
  }
  merge<-rbedintersect(bed1,ref)
  return(merge)
}

