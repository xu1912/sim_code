rpath <- "~/test";
setwd (rpath);
#set.seed(5*as.numeric($SLURM_ARRAY_TASK_ID));

library(lme4)
library(MASS)
source("adjPrem.R")
f_O <- paste(rpath, "/result/result_n10_$SLURM_ARRAY_TASK_ID.txt", sep="")


n_id=10
n_prob_intron=5
n_prob_exon=3

prem=1
truem=0.5

v_id=5
v_prob=3
v_err=2

est_1=c()
est_2=c()
est_3=c()

for(jk in 1:500){

re_id=rnorm(n_id,mean=0,sd=sqrt(v_id))
re_prob_intron=rnorm(n_prob_intron,mean=0,sd=sqrt(v_prob))
re_err_intron=rnorm(n_id*n_prob_intron, mean=0, sd=sqrt(v_err))

idx=0
prb_intron=matrix(0,nrow=n_id,ncol=n_prob_intron)
for (i in 1:n_id){
        for (j in 1:n_prob_intron){
                idx=idx+1
                prb_intron[i,j]=prem+re_id[i]+re_prob_intron[j]+re_err_intron[idx]
        }
}

dimnames(prb_intron) <- list(rownames(prb_intron, do.NULL = FALSE, prefix = "id_"),
                          colnames(prb_intron, do.NULL = FALSE, prefix = "prob_intron"))

d_intron=data.frame(prb_intron)

re_prob_exon=rnorm(n_prob_exon,mean=0,sd=sqrt(v_prob))
re_err_exon=rnorm(n_id*n_prob_exon, mean=0, sd=sqrt(v_err))
idx=0
prb_exon=matrix(0,nrow=n_id,ncol=n_prob_exon)
for (i in 1:n_id){
        for (j in 1:n_prob_exon){
                idx=idx+1
                prb_exon[i,j]=truem+prem+re_id[i]+re_prob_exon[j]+re_err_exon[idx]
        }
}

dimnames(prb_exon) <- list(rownames(prb_exon, do.NULL = FALSE, prefix = "id_"),
                          colnames(prb_exon, do.NULL = FALSE, prefix = "prob_exon"))



gexp=rep(0, n_id*n_prob_intron)
id_x=rep("", n_id*n_prob_intron)
prb_x=rep("", n_id*n_prob_intron)
exon_t=rep(0, n_id*n_prob_intron)


for(t in 1:n_id){

        for (tt in 1:n_prob_intron){
                gexp[(t-1)*n_prob_intron+tt]=prb_intron[t,tt]
                id_x[(t-1)*n_prob_intron+tt]=rownames(prb_intron)[t]
                prb_x[(t-1)*n_prob_intron+tt]=colnames(prb_intron)[tt]
        }

}

dd=data.frame(gexp,id_x,prb_x,exon_t)

gexp=rep(0, n_id*n_prob_exon)
id_x=rep("", n_id*n_prob_exon)
prb_x=rep("", n_id*n_prob_exon)
exon_t=rep(1, n_id*n_prob_exon)


for(t in 1:n_id){

        for (tt in 1:n_prob_exon){
                gexp[(t-1)*n_prob_exon+tt]=prb_exon[t,tt]
                id_x[(t-1)*n_prob_exon+tt]=rownames(prb_exon)[t]
                prb_x[(t-1)*n_prob_exon+tt]=colnames(prb_exon)[tt]
        }

}

dd_exon=data.frame(gexp,id_x,prb_x,exon_t)

dd=rbind(dd,dd_exon)

ddresult=lmer(gexp~1+(1|id_x)+(1|prb_x),data=dd)
est_1[jk]=fixef(ddresult)-truem


dresult=lmer(gexp~1+exon_t+(1|id_x)+(1|prb_x),data=dd)
est_2[jk]=fixef(dresult)[2]-truem


est_3[jk]=adj_prem_rawInput(prb_intron,prb_exon)-truem

}


write.table(data.frame(est_1, est_2, est_3), f_O, append=FALSE, sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE);
