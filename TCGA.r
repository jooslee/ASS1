##### Custom code for Keshet, Lee et al. Targeting purine synthesis in ASS1 expressing cancers promotes response to immune checkpoint inhibitors

##### Environment setting 
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

library(parallel)
library(survival)
library("survminer")
library(data.table)

qnorm.array <- function(mat)     
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}

rank.array <- function(mat)     
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average")/(length(mat)+1);
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}

##### Load TCGA data
#setwd("/data/leej55/umiacs/project-scratch/srescues/newmol")
load("./data/prob.TCGA.RData")

##### Association between ASS1 expression and glycolysis/gluconeogenesis 
it=which(prob$types %in% c("THCA","STAD","LUAD","LUSC","DLBC","COAD","OV","CESC","PAAD","ESCA","BRCA"))
pval=NA;iass1=which(prob$genes=="ASS1")
gg=c("PC","PCK1","FBP1","HK1")
for (i in seq(length(gg)))
pval[i]=wilcox.test(prob$mRNA[iass1,prob$mRNAq2[which(prob$genes==gg[i]),]==2],
					prob$mRNA[iass1,prob$mRNAq2[which(prob$genes==gg[i]),]==0],alternative="greater")$p.value
names(pval)=gg

##### TCGA survival analysis
sig.genes=c("PHGDH","PSAT1","PSPH","SHMT1","SHMT2")
nds.score=apply(prob$mRNA.norm2[which(prob$genes %in% c(sig.genes)),],2,mean,na.rm=T)
surv.dt=prob$surv.dt
types=prob$types

ndsq2=rep(NA,length(prob$samples))
for (i in seq(length(utypes))){
	it=which(utypes[i]==prob$types)
	ndsq2[it]=(nds.score[it]>quantile(nds.score[it],1/2,na.rm=T))*1+(nds.score[it]>quantile(nds.score[it],1/2,na.rm=T))*1
}


##### ASS1 
itx=match(c("BRCA","LUAD","COAD"),utypes)
q2=rep(NA,length(prob$samples))
for (i in seq(length(utypes))){
	it=which(utypes[i]==prob$types)
	q2[it]=(prob$mRNA[which(prob$genes=="ASS1"),it]>quantile(prob$mRNA[which(prob$genes=="ASS1"),it],2/3,na.rm=T))*1+
    		(prob$mRNA[which(prob$genes=="ASS1"),it]>quantile(prob$mRNA[which(prob$genes=="ASS1"),it],2/3,na.rm=T))*1
}
dt1=data.frame(prob$surv.dt,ass1=q2,ASS1=prob$mRNA[which(prob$genes=="ASS1"),],types=prob$types,
				samples=prob$samples,nds=ndsq2)
dt1=dt1[dt1$ass1!=1,]
dt0=dt1

pval=msd=hr=NA
for (i in seq(length(utypes))){
	dt1=dt0[dt0$types %in% utypes[i],]
	fit <- survfit(Surv(time, status) ~ ass1, data = dt1)
	tst=survdiff(Surv(time, status) ~ ass1,data=dt1)
	logrank.p <- 1 - pchisq(tst$chisq, length(tst$n) - 1)
	pval[i]=logrank.p
	msd[i]=summary(fit)$table[2,7]-summary(fit)$table[1,7]
	hr[i]=tst$obs[1]/tst$exp[1]-tst$obs[2]/tst$exp[2]
}
names(pval)=names(msd)=names(hr)=utypes

##### ASS1 in high ND score samples
pval=msd=hr=rep(NA,length(utypes))
for (i in seq(length(utypes))){
	dt1=data.frame(prob$surv.dt,ass1=q2,ASS1=prob$mRNA[which(prob$genes=="ASS1"),],types=prob$types,samples=prob$samples,ndsq2=ndsq2)
	dt1=dt0[dt0$types %in% utypes[i],]
	dt1$ndsq2=(dt1$nds>quantile(dt1$nds,1/2,na.rm=T))*1+(dt1$nds>quantile(dt1$nds,1/2,na.rm=T))*1
	dt1=dt1[dt1$ndsq2==2,]
	if (utypes[i] %in% c("BRCA","COAD")) dt1$ass1=(dt1$ASS1>quantile(dt1$ASS1,1/2,na.rm=T))*1+(dt1$ASS1>quantile(dt1$ASS1,1/2,na.rm=T))*1
	if (utypes[i] %in% c("LUAD")) dt1$ass1=(dt1$ASS1>quantile(dt1$ASS1,2/3,na.rm=T))*1+(dt1$ASS1>quantile(dt1$ASS1,2/3,na.rm=T))*1
	dt1=dt1[dt1$ass1!=1,]

	fit <- survfit(Surv(time, status) ~ ass1, data = dt1)
	tst=survdiff(Surv(time, status) ~ ass1,data=dt1)
	logrank.p <- 1 - pchisq(tst$chisq, length(tst$n) - 1)
	if (i %in% itx) {
		pval[i]=logrank.p
		msd[i]=summary(fit)$table[2,7]-summary(fit)$table[1,7]
		hr[i]=tst$obs[1]/tst$exp[1]-tst$obs[2]/tst$exp[2]
	}
}
names(pval)=names(msd)=names(hr)=utypes

##### Association between ND signature and ASS1
sig.genes=c("PHGDH","PSAT1","PSPH","SHMT1","SHMT2")
alpha=1/5

pval=rep(NA,length(utypes))
for (i in seq(length(utypes))) if (utypes[i] %in% c("BRCA","LUAD","LUSC","COAD")){
	it=which(prob$types %in% utypes[i])
	sig.score=apply(prob$mRNA[which(prob$genes %in% sig.genes),it],2,sum,na.rm=T)
	ass1=prob$mRNA[which(prob$genes=="ASS1"),it]

	q1=quantile(sig.score,alpha,na.rm=T)
	q2=quantile(sig.score,1-alpha,na.rm=T)
	pval[i]=wilcox.test(ass1[sig.score>q2],ass1[sig.score<q1],alternative="greater")$p.value
}
names(pval)=utypes

##### Association between gluconeogenetic enzymes and ASS1+ASL+NOS1
neo.genes=c("PC", "PCK1", "PCK2", "FBP1", "FBP2")
tps=list()
tps[1]="BRCA"
tps[2]=c("LUSC","LUAD")
tps[3]="COAD"
pval=NA
for (j in c(1:3)){
	it.lung=which(prob$types %in% tps[[j]])
	score=apply(prob$mRNA.rank2[which(prob$genes %in% sel.genes),it.lung],2,sum,na.rm=T)
	neo.score=apply(prob$mRNA.rank2[which(prob$genes %in% neo.genes),it.lung],2,sum,na.rm=T)
	q1=quantile(score,0.33,na.rm=T)
	q2=quantile(score,0.66,na.rm=T)
	pval[j]=wilcox.test(neo.score[score<q1],neo.score[score>q2],alternative="less")$p.value
}






