adj_prem_rawInput=function(dt,dtx){

ta=(as.matrix(dt))
taex=(as.matrix(dtx))
nid=length(ta[,1])
npid=length(ta[1,])

gexp=rep(0, nid*npid)
id_x=rep("", nid*npid)
prb_x=rep("", nid*npid)

for(t in 1:nid){

for (tt in 1:npid){
gexp[(t-1)*npid+tt]=ta[t,tt]
id_x[(t-1)*npid+tt]=rownames(ta)[t]
prb_x[(t-1)*npid+tt]=colnames(ta)[tt]
}

}

dd=data.frame(gexp,id_x,prb_x)
dresult=lmer(gexp~1+(1|id_x)+(1|prb_x),data=dd)
prem_ef=fixef(dresult)
taex=taex[order(rownames(taex)),]

if(length(dtx[,1])>1){
rd=taex-prem_ef
}else{
rd=taex-prem_ef
}
#rd[rd<0]=0


nid=length(rd[,1])
npid=length(rd[1,])

gexp=rep(0, nid*npid)
id_x=rep("", nid*npid)
prb_x=rep("", nid*npid)

for(t in 1:nid){

for (tt in 1:npid){
gexp[(t-1)*npid+tt]=rd[t,tt]
id_x[(t-1)*npid+tt]=rownames(rd)[t]
prb_x[(t-1)*npid+tt]=colnames(rd)[tt]
}

}

dd=data.frame(gexp,id_x,prb_x)

ddresult=lmer(gexp~1+(1|id_x)+(1|prb_x),data=dd)

return(fixef(ddresult))

}
