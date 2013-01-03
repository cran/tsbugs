nodes <-
function(bug, part=NULL){
  if(class(bug)!="tsbugs")
    stop("bug must be object with tsbugs")
  if(!is.null(part))
    if(is.na(pmatch(part, names(bug$info))))
      stop("part must be a component of tsbug$info such as likelihood, priors, etc...")
  if(!is.null(part))
    bug<-bug$bug[bug$info[[pmatch(part, names(bug$info))]]]
  else   bug<-bug$bug
  
  nds<-gsub("\t","",bug)
  nds<-gsub("\\[t\\]","",nds)

  st.nds<-grep("~",nds)
  st<-nds[st.nds]
  st.dist<-gsub("(.*)~","",st)
  st.dist<-gsub("\\s","",st.dist)
  st<-gsub("~(.*)","",st)
  st<-gsub("\\s","",st)
  if(length(st)==0)  out<-NULL
  if(length(st)>0)  out<-data.frame(name=st,type="~",dist.trans=st.dist,beg=NA,end=NA,stoc=1,id=st.nds)

  det.nds<-grep("<-",nds)
  det<-nds[det.nds]
  det.trans<-gsub("(.*)<-","",det)
  det.trans<-gsub("\\s","",det.trans)
  det<-gsub("<-(.*)","",det)
  det<-gsub("\\s","",det)
  if(length(det)>0)  out<-rbind(out,data.frame(name=det,type="<-",dist.trans=det.trans,beg=NA,end=NA,stoc=0,id=det.nds))

  st.loops<-grep("[0-9]:" ,nds)
  end.loops<-grep("}" ,nds)
  rg<-nds[st.loops]
  rg<-strsplit(gsub("[^[:digit:]:[:digit:]]","",rg),":")
  if(length(st.loops)>0){
    for(i in 1:dim(out)[1]){
      sel.node<-out$id[i]
      for(j in 1:length(st.loops)){
        if(sel.node>st.loops[j] & sel.node<end.loops[j]){
          out$beg[i]<-as.numeric(rg[[j]][1])
          out$end[i]<-as.numeric(rg[[j]][2])
        }
      }
    } 
  }
  out$name<-as.character(out$name)
  out
}
