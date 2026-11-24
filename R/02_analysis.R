# 02 analysis: multiple imputation (MICE, m = 20) as the main analysis; primary, secondary, confounding control,
# sensitivity, subgroup and diagnostic analyses
suppressMessages({library(data.table);library(survival);library(splines);library(WeightIt);library(MatchIt);library(mice);library(pROC);library(EValue);library(sandwich);library(lmtest);library(cmprsk);library(metafor);library(logistf);library(car)})
set.seed(2026); M<-20; LM<-96
D<-as.data.table(readRDS("D.rds")); S<-as.data.table(readRDS("S.rds")); R<-list()
D[db=="eICU-CRD",bsi:=NA]  # eICU-CRD microbiology table too sparse for bloodstream-infection ascertainment
prep<-function(d){d<-copy(d); d[,sevz:=(sev-mean(sev,na.rm=TRUE))/sd(sev,na.rm=TRUE),by=db]
 d[,`:=`(lalc0=log(pmax(alc0,0.05)),lpf=log(pmin(pmax(pf,20),700)),lanc=log(pmax(anc0,0.1)),old=as.integer(age>=65),hisev=as.integer(sevz>0),sevlym=as.integer(alc0<0.5),DB=as.integer(db=="eICU-CRD"))]; d}
D<-prep(D)
num<-c("sevz","lact","creat","bun","hb","wbc","plt","map","hr","rr","lpf")
M2<-"age+male+chf+ckd+copd+dm+liver+cancer+sepsis"; M3<-paste(M2,"+sevz+lalc0+lact+creat+bun+hb+wbc+plt+map+hr+rr+lpf+vaso+imv+ster+ards")
cvp<-function(cv,d) if(length(unique(d$db))>1) paste(cv,"+DB") else cv
# ---------------- multiple imputation, separately in each database ----------------
iv<-c("d28","x","age","male","chf","ckd","copd","dm","liver","cancer","sepsis","ards","vaso","imv","ster","lalc0",num,"lanc","alb","bili")
mi1<-function(d){z<-d[,..iv]; ini<-mice(z,maxit=0,printFlag=FALSE); pm<-ini$predictorMatrix; meth<-ini$method; meth[names(meth)%in%c("d28","x")]<-""
 mice(z,m=M,maxit=10,method=meth,predictorMatrix=pm,printFlag=FALSE,seed=2026)}
dm<-D[db=="MIMIC-IV"]; de<-D[db=="eICU-CRD"]; IM<-mi1(dm); IE<-mi1(de)
fill<-function(d,imp,i){d<-copy(d); cd<-as.data.table(complete(imp,i)); for(v in setdiff(names(cd),c("d28","x"))) set(d,j=v,value=cd[[v]]); d}
IMP<-list("MIMIC-IV"=lapply(1:M,function(i) fill(dm,IM,i)),"eICU-CRD"=lapply(1:M,function(i) fill(de,IE,i)))
IMP$Pooled<-lapply(1:M,function(i) rbind(IMP[["MIMIC-IV"]][[i]],IMP[["eICU-CRD"]][[i]]))
for(k in names(IMP)) for(i in 1:M){d<-IMP[[k]][[i]]; d[,hisev:=as.integer(sevz>0)]; IMP[[k]][[i]]<-d}
saveRDS(IMP,"IMP.rds")
# ---------------- helpers: Rubin's rules ----------------
rub<-function(est,se){q<-mean(est); u<-mean(se^2); b<-var(est); t<-u+(1+1/length(est))*b; dfm<-(length(est)-1)*(1+u/((1+1/length(est))*b))^2; if(!is.finite(dfm)) dfm<-1e6; c(est=q,se=sqrt(t),p=2*pt(-abs(q/sqrt(t)),dfm))}
ex<-function(r) c(exp(r[1]),exp(r[1]-1.96*r[2]),exp(r[1]+1.96*r[2]),r[3])
fm<-function(v) sprintf("%.2f (%.2f–%.2f)",v[1],v[2],v[3])
cf<-function(fit,tr,V=NULL){if(is.null(V)) V<-vcov(fit); b<-coef(fit); i<-grep(tr,names(b))[1]; c(b[i],sqrt(V[i,i]))}
plg<-function(k,f,tr="^x$",sub=NULL){r<-sapply(IMP[[k]],function(d){if(!is.null(sub)) d<-d[eval(parse(text=sub))]; cf(glm(as.formula(f),binomial,d),tr)}); ex(rub(r[1,],r[2,]))}
pcx<-function(k,cv,tr="^x$",sub=NULL){st<-if(k=="Pooled") "+strata(DB)" else ""; r<-sapply(IMP[[k]],function(d){if(!is.null(sub)) d<-d[eval(parse(text=sub))]; cf(coxph(as.formula(paste("Surv(t28,d28)~x+",cv,st)),d),tr)}); ex(rub(r[1,],r[2,]))}
# ---------------- primary estimate by database ----------------
p0<-plg("MIMIC-IV",paste("d28~x+",M3)); p1<-plg("eICU-CRD",paste("d28~x+",M3))
cat("Model 3 OR (MICE m = 20)  MIMIC-IV",fm(p0),"P =",signif(p0[4],3)," | eICU-CRD",fm(p1),"P =",signif(p1[4],3),"\n")
dbs<-c("MIMIC-IV","eICU-CRD","Pooled")
# ---------------- Table 1 (observed data) ----------------
vb<-c("age","male","chf","copd","ckd","dm","liver","cancer","sepsis","ards","asp","alc0","alc3","anc0","wbc","lact","creat","bun","hb","plt","alb","pf","hr","map","rr","sev","vaso","imv","ster","nalc","d28","hd")
raw<-as.data.table(readRDS("D.rds"))
R$t0<-rbindlist(lapply(dbs[1:2],function(k) rbindlist(lapply(0:1,function(g){z<-raw[db==k&x==g]; data.table(db=k,grp=g,n=nrow(z),t(sapply(vb,function(v){u<-z[[v]]; if(all(u%in%c(0,1,NA))) sprintf("%.1f%%",100*mean(u,na.rm=TRUE)) else sprintf("%.2f (%.2f–%.2f)",median(u,na.rm=TRUE),quantile(u,.25,na.rm=TRUE),quantile(u,.75,na.rm=TRUE))})))}))))
R$t0p<-sapply(vb,function(v) sapply(dbs[1:2],function(k){z<-raw[db==k]; if(all(is.na(z[[v]]))) NA else if(all(z[[v]]%in%c(0,1,NA))) chisq.test(table(z[[v]],z$x))$p.value else wilcox.test(z[[v]]~z$x)$p.value}))
R$miss<-raw[,c(list(Database=db[1]),lapply(.SD,function(z) sprintf("%.1f",100*mean(is.na(z))))),by=db,.SDcols=c("lact","creat","bun","hb","wbc","plt","pf","hr","map","rr","sev","anc0","alb","bili")][,-1]
# ---------------- Primary outcome ----------------
R$t1a<-rbindlist(lapply(dbs,function(k){d<-IMP[[k]][[1]]; o3<-plg(k,paste("d28~x+",cvp(M3,d))); h3<-pcx(k,M3)
 data.table(db=k,grp=0:1,n=c(sum(d$x==0),sum(d$x==1)),ev=c(sum(d$d28[d$x==0]),sum(d$d28[d$x==1])),crude=c("1",fm(plg(k,"d28~x"))),m2=c("1",fm(plg(k,paste("d28~x+",cvp(M2,d))))),m3=c("1",fm(o3)),p3=c(NA,o3[4]),hr=c("1",fm(h3)),ph=c(NA,h3[4]))}))
# adjusted (marginally standardized) 28-day mortality: 100 bootstrap samples within each imputed data set, combined by Rubin's rules
stdz<-function(d,cv){f<-glm(as.formula(paste("d28~x+",cv)),binomial,d); d1<-copy(d); d1[,x:=1L]; d0<-copy(d); d0[,x:=0L]; c(mean(predict(f,d1,type="response")),mean(predict(f,d0,type="response")))}
R$std<-rbindlist(lapply(dbs,function(k){cv<-cvp(M3,IMP[[k]][[1]]); r<-lapply(IMP[[k]],function(d){e<-stdz(d,cv); b<-t(replicate(100,stdz(d[sample(.N,replace=TRUE)],cv))); c(e,e[1]-e[2],apply(cbind(b,b[,1]-b[,2]),2,sd))})
 r<-do.call(rbind,r); g<-function(j){q<-rub(r[,j],r[,j+3]); c(q[1],q[1]-1.96*q[2],q[1]+1.96*q[2])}; a<-g(1); b<-g(2); c3<-g(3)
 data.table(db=k,p1=a[1],p1lo=a[2],p1hi=a[3],p0=b[1],p0lo=b[2],p0hi=b[3],rd=c3[1],rdlo=c3[2],rdhi=c3[3])}))
R$t1b<-rbindlist(lapply(dbs,function(k){cv<-cvp(M3,IMP[[k]][[1]]); for(i in 1:M){d<-IMP[[k]][[i]]; d[,xd:=-pmin(alc3,3)/0.1]; d[,xc:=relevel(factor(cut(alc3,c(-Inf,0.5,0.8,1,Inf),right=FALSE,labels=c("lt05","c05","c08","ge10"))),"ge10")]; d[,xt:=as.integer(cut(alc3,c(-Inf,0.5,0.8,1,Inf),right=FALSE))]; IMP[[k]][[i]]<<-d}
 a<-plg(k,paste("d28~xd+",cv),"^xd$"); h<-{st<-if(k=="Pooled") "+strata(DB)" else ""; r<-sapply(IMP[[k]],function(d) cf(coxph(as.formula(paste("Surv(t28,d28)~xd+",cv,st)),d),"^xd$")); ex(rub(r[1,],r[2,]))}
 cs<-lapply(c("xclt05","xcc05","xcc08"),function(tr) plg(k,paste("d28~xc+",cv),paste0("^",tr,"$")))
 t<-plg(k,paste("d28~xt+",cv),"^xt$")
 rbind(data.table(db=k,a="Per 0.1 ×10⁹/L lower ALC at 48–96 h, OR (Model 3)",v=fm(a),p=a[4]),data.table(db=k,a="Per 0.1 ×10⁹/L lower ALC at 48–96 h, HR (Cox Model 3)",v=fm(h),p=h[4]),
  data.table(db=k,a=paste0("ALC 48–96 h ",c("<0.5","0.5-0.79","0.8-0.99")," vs ≥1.0 ×10⁹/L, OR"),v=sapply(cs,fm),p=sapply(cs,`[`,4)),data.table(db=k,a="P for trend (4 categories)",v="—",p=t[4]))}))
# spline (Wald tests pooled with D1), Kaplan-Meier, incremental discrimination
# restricted cubic spline basis (Harrell parameterisation): linear term plus k-2 non-linear terms
rcsb<-function(x,k){nk<-length(k); nr<-(k[nk]-k[1])^2; cbind(x,matrix(sapply(1:(nk-2),function(j) (pmax(x-k[j],0)^3-pmax(x-k[nk-1],0)^3*(k[nk]-k[j])/(k[nk]-k[nk-1])+pmax(x-k[nk],0)^3*(k[nk-1]-k[j])/(k[nk]-k[nk-1]))/nr),nrow=length(x)))}
R$t1c<-rbindlist(lapply(dbs,function(k){cv<-cvp(M3,IMP[[k]][[1]]); st<-if(k=="Pooled") "+strata(DB)" else ""; d<-IMP[[k]][[1]]
 lr<-survdiff(as.formula(paste("Surv(t28,d28)~x",st)),d); plr<-pchisq(lr$chisq,1,lower.tail=FALSE)
 lo<-quantile(d$alc3,.01); hi<-quantile(d$alc3,.99); k4<-quantile(pmin(pmax(d$alc3,lo),hi),c(.05,.35,.65,.95))
 rb<-function(d){d<-copy(d); B<-rcsb(pmin(pmax(d$alc3,lo),hi),k4); d[,`:=`(r1=B[,1],r2=B[,2],r3=B[,3])]; d}
 f1<-as.mira(lapply(IMP[[k]],function(d) glm(as.formula(paste("d28~r1+r2+r3+",cv)),binomial,rb(d))))
 fl<-as.mira(lapply(IMP[[k]],function(d) glm(as.formula(paste("d28~r1+",cv)),binomial,rb(d))))
 f0<-as.mira(lapply(IMP[[k]],function(d) glm(as.formula(paste("d28~",cv)),binomial,d)))
 po<-D1(f1,f0)$result[1,4]; pnl<-D1(f1,fl)$result[1,4]
 au<-sapply(IMP[[k]],function(d) c(auc(d$d28,fitted(glm(as.formula(paste("d28~",cv)),binomial,d)),quiet=TRUE),auc(d$d28,fitted(glm(as.formula(paste("d28~x+",cv)),binomial,d)),quiet=TRUE)))
 data.table(db=k,a=c("Kaplan–Meier log-rank (from 96-h landmark)","RCS of ALC at 48–96 h (4 knots), overall P (pooled Wald test)","RCS, P for non-linearity","Adding persistent lymphopenia to Model 3, mean AUC across imputations"),v=c("—","—","—",sprintf("AUC %.3f → %.3f",mean(au[1,]),mean(au[2,]))),p=c(plr,po,pnl,NA))}))
# ---------------- Secondary outcomes ----------------
so<-function(k,lab,y,type){d0<-IMP[[k]][[1]]; if(all(is.na(d0[[y]]))) return(NULL); cv<-cvp(M3,d0); fr<-sprintf("%.1f%%",100*mean(d0[[y]],na.rm=TRUE))
 r<-sapply(IMP[[k]],function(d){d<-d[!is.na(get(y))]
  if(type=="lg") cf(glm(as.formula(paste(y,"~x+",cv)),binomial,d),"^x$")
  else if(type=="cox"){d[,tt:=fifelse(get(y)==1,death_h,ifelse(y=="d90",2160,8760))-LM]; cf(coxph(as.formula(paste("Surv(tt,",y,")~x+",cv)),d),"^x$")}
  else if(type=="fg"){ev<-fifelse(d$bsi==1,1L,fifelse(d$d28==1,2L,0L)); tt<-fifelse(ev==2,d$death_h-LM,672-LM); X<-model.matrix(as.formula(paste("~x+",cv)),d)[,-1]; f<-crr(tt,ev,X,failcode=1,cencode=0); c(f$coef["x"],sqrt(f$var[1,1]))}
  else cf(lm(as.formula(paste(y,"~x+",cv)),d),"^x$")})
 q<-rub(r[1,],r[2,])
 if(type=="lin") data.table(db=k,o=lab,f=sprintf("%.1f (%.1f–%.1f)",median(d0[[y]],na.rm=TRUE),quantile(d0[[y]],.25,na.rm=TRUE),quantile(d0[[y]],.75,na.rm=TRUE)),m="Linear Model 3, β (days)",v=sprintf("%.2f (%.2f to %.2f)",q[1],q[1]-1.96*q[2],q[1]+1.96*q[2]),p=q[3])
 else data.table(db=k,o=lab,f=fr,m=c(lg="Logistic Model 3, OR",cox="Cox Model 3, HR",fg="Fine–Gray (death = competing), SHR Model 3")[type],v=fm(ex(q)),p=q[3])}
R$t2<-rbindlist(list(so("MIMIC-IV","In-hospital mortality","hd","lg"),so("eICU-CRD","In-hospital mortality","hd","lg"),so("Pooled","In-hospital mortality","hd","lg"),
 so("MIMIC-IV","90-day mortality","d90","cox"),so("MIMIC-IV","1-year mortality","d365","cox"),so("MIMIC-IV","Bloodstream infection (96 h–day 28)","bsi","lg"),so("MIMIC-IV","Bloodstream infection (96 h–day 28)","bsi","fg"),
 so("MIMIC-IV","Ventilator-free days to day 28","vfd","lin"),so("MIMIC-IV","ICU length of stay","los","lin"),so("eICU-CRD","ICU length of stay","los","lin"),so("Pooled","ICU length of stay","los","lin")),fill=TRUE)
# ---------------- Confounding control (each method applied within every imputed data set; Rubin's rules) ----------------
smd<-function(d,ps,wt){X<-model.matrix(ps,d)[,-1]; max(abs(apply(X,2,function(z){m1<-weighted.mean(z[d$x==1],wt[d$x==1]);m0<-weighted.mean(z[d$x==0],wt[d$x==0]); (m1-m0)/sqrt((var(z[d$x==1])+var(z[d$x==0]))/2)})))}
R$t3<-rbindlist(lapply(dbs,function(k){cv<-cvp(M3,IMP[[k]][[1]]); ps<-as.formula(paste("x~",cv)); st<-if(k=="Pooled") "+strata(DB)" else ""
 r<-lapply(IMP[[k]],function(d){w1<-weightit(ps,d,method="glm",estimand="ATE",stabilize=TRUE); d[,iw:=pmin(w1$weights,quantile(w1$weights,.99))]; w2<-weightit(ps,d,method="glm",estimand="ATO"); d[,ow:=w2$weights]
  fw<-function(f,w){m<-suppressWarnings(glm(as.formula(f),quasibinomial,d,weights=w)); cf(m,"^x$",vcovHC(m,type="HC0"))}
  ch<-coxph(as.formula(paste("Surv(t28,d28)~x",st)),d,weights=iw,robust=TRUE)
  mt<-matchit(ps,d,method="nearest",distance="glm",link="linear.logit",caliper=0.2,ratio=1); md<-as.data.table(match.data(mt))
  g1<-glm(d28~x,binomial,md); g2<-glm(as.formula(paste("d28~x+",cv)),binomial,md)
  out<-rbind(iptw=fw("d28~x",d$iw),iptwh=c(coef(ch)["x"],sqrt(vcov(ch)["x","x"])),dr=fw(paste("d28~x+",cv),d$iw),ow=fw("d28~x",d$ow),psm=cf(g1,"^x$",vcovCL(g1,cluster=~subclass)),psmdr=cf(g2,"^x$",vcovCL(g2,cluster=~subclass)))
  if(k=="eICU-CRD"){g3<-glm(as.formula(paste("d28~x+",cv)),binomial,d); out<-rbind(out,hosp=cf(g3,"^x$",vcovCL(g3,cluster=~hosp)))}
  list(out=out,np=sum(md$x==0),smd=c(smd(d,ps,rep(1,nrow(d))),smd(d,ps,d$iw),smd(d,ps,d$ow)))})
 nm<-rownames(r[[1]]$out); q<-lapply(nm,function(j) ex(rub(sapply(r,function(z) z$out[j,1]),sapply(r,function(z) z$out[j,2])))); names(q)<-nm
 o2<-plg(k,paste("d28~x+",cvp(M2,IMP[[k]][[1]]))); o3<-plg(k,paste("d28~x+",cv)); ev<-evalues.OR(o3[1],o3[2],o3[3],rare=mean(IMP[[k]][[1]]$d28)<0.15)
 lab<-c(iptw="IPTW (ATE, stabilized, truncated 99th pct), OR",iptwh="IPTW, HR (Cox, robust SE)",dr="Doubly robust (IPTW + Model 3), OR",ow="Overlap weights (ATO), OR",psm=sprintf("PSM 1:1, caliper 0.2 SD of logit PS (median %d pairs), OR",as.integer(median(sapply(r,`[[`,"np")))),psmdr="PSM + Model 3 covariate adjustment in matched sample, OR",hosp="Model 3 OR, hospital-clustered robust SE")
 sm<-colMeans(t(sapply(r,`[[`,"smd")))
 rbind(data.table(db=k,a=c("Model 2, OR","Model 3, OR"),v=c(fm(o2),fm(o3)),p=c(o2[4],o3[4])),data.table(db=k,a=lab[nm],v=sapply(q,fm),p=sapply(q,`[`,4)),
  data.table(db=k,a="E-value (Model 3 OR; CI bound)",v=sprintf("%.2f (%.2f)",ev[2,1],ev[2,2]),p=NA),data.table(db=k,a=sprintf("Max absolute SMD: unweighted %.3f; IPTW %.3f; overlap %.3f",sm[1],sm[2],sm[3]),v="—",p=NA))}))
# two-stage (database-level) meta-analysis of the model 3 log ORs
mm<-rbindlist(lapply(dbs[1:2],function(k){r<-sapply(IMP[[k]],function(d) cf(glm(as.formula(paste("d28~x+",M3)),binomial,d),"^x$")); q<-rub(r[1,],r[2,]); data.table(db=k,yi=q[1],sei=q[2])}))
fe<-rma(yi=mm$yi,sei=mm$sei,method="FE"); re<-rma(yi=mm$yi,sei=mm$sei,method="REML")
R$meta<-data.table(a=c("Two-stage meta-analysis, fixed effect, OR","Two-stage meta-analysis, random effects (REML), OR","Heterogeneity between databases, I² (%)"),v=c(fm(exp(c(fe$b,fe$ci.lb,fe$ci.ub))),fm(exp(c(re$b,re$ci.lb,re$ci.ub))),sprintf("%.1f",re$I2)),p=c(fe$pval,re$pval,re$QEp))
# ---------------- Sensitivity analyses ----------------
SA<-list(); add<-function(k,lab,sub=NULL,f="x",cv=M3,mk=NULL){if(!is.null(mk)) for(i in 1:M){d<-IMP[[k]][[i]]; eval(parse(text=mk)); IMP[[k]][[i]]<<-d}
 d1<-IMP[[k]][[1]]; if(!is.null(sub)) d1<-d1[eval(parse(text=sub))]; r<-plg(k,paste("d28~",f,"+",cvp(cv,d1)),paste0("^",f,"$"),sub)
 SA[[length(SA)+1]]<<-data.table(db=k,a=lab,n=nrow(d1),ev=sum(d1$d28),v=fm(r),p=r[4])}
for(k in dbs){
 add(k,"Persistent threshold <0.8 ×10⁹/L",f="x2",mk="d[,x2:=as.integer(alc3<0.8)]"); add(k,"Exposure = lowest (not last) ALC 48–96 h <1.0",f="x3",mk="d[,x3:=as.integer(alcmin3<1)]")
 add(k,"Exposure window 48–72 h",sub="!is.na(alc72)",f="x4",mk="d[,x4:=as.integer(alc72<1)]")
 add(k,"Model 3 + number of ALC measurements 0–96 h (measurement intensity)",cv=paste(M3,"+log(nalc)")); add(k,"Model 3 + baseline neutrophil count",cv=paste(M3,"+lanc")); add(k,"Model 3 + albumin and bilirubin",cv=paste(M3,"+alb+bili"))
 add(k,"Excluding systemic corticosteroids 0–96 h",sub="ster==0",cv=sub("\\+ster","",M3)); add(k,"Excluding sepsis (non-septic pneumonia)",sub="sepsis==0",cv=sub("\\+sepsis","",M3))
 add(k,"Excluding solid cancer",sub="cancer==0",cv=sub("\\+cancer","",M3)); add(k,"Excluding aspiration pneumonia",sub="asp==0")
 add(k,"Pneumonia as principal diagnosis (MIMIC-IV: top-3 discharge code; eICU-CRD: APACHE admission diagnosis)",sub="pna3==1")
 add(k,"Excluding death before day 7",sub="is.na(death_h)|death_h>168"); add(k,"Invasive ventilation 0–24 h only",sub="imv==1",cv=sub("\\+imv","",M3)); add(k,"No invasive ventilation 0–24 h",sub="imv==0",cv=sub("\\+imv","",M3))}
SA<-rbindlist(SA)
# single median imputation and complete-case analyses
med<-copy(D); for(v in c(num,"lalc0")) med[,(v):=fifelse(is.na(get(v)),median(get(v),na.rm=TRUE),get(v)),by=db]
for(k in dbs){z<-if(k=="Pooled") med else med[db==k]; o<-glm(as.formula(paste("d28~x+",cvp(M3,z))),binomial,z); r<-ex(c(cf(o,"^x$"),coef(summary(o))["x",4])); SA<-rbind(SA,data.table(db=k,a="Single median imputation",n=nrow(z),ev=sum(z$d28),v=fm(r),p=r[4]))
 z<-if(k=="Pooled") D else D[db==k]; z<-z[complete.cases(z[,c(num,"lalc0"),with=FALSE])]; o<-glm(as.formula(paste("d28~x+",cvp(M3,z))),binomial,z); r<-ex(c(cf(o,"^x$"),coef(summary(o))["x",4])); SA<-rbind(SA,data.table(db=k,a="Complete-case analysis",n=nrow(z),ev=sum(z$d28),v=fm(r),p=r[4]))}
# inverse probability of selection weighting (patients without a 48-96 h ALC)
S<-prep(S); for(v in c(num,"lalc0")) S[,(v):=fifelse(is.na(get(v)),median(get(v),na.rm=TRUE),get(v)),by=db]
selv<-"age+male+chf+ckd+copd+dm+liver+cancer+sepsis+ards+sevz+lalc0+lact+creat+bun+hb+wbc+plt+map+hr+rr+lpf+vaso+imv+ster"
for(k in dbs[1:2]){s<-S[db==k]; g<-glm(as.formula(paste("inc~",selv)),binomial,s); s[,ps:=fitted(g)]; w<-1/s[inc==1]$ps; w<-pmin(w,quantile(w,.99)); ids<-s[inc==1]$id
 r<-sapply(IMP[[k]],function(d){d<-copy(d); d[,sw:=w[match(id,ids)]]; m<-suppressWarnings(glm(as.formula(paste("d28~x+",M3)),quasibinomial,d,weights=sw)); cf(m,"^x$",vcovHC(m,type="HC0"))}); q<-ex(rub(r[1,],r[2,]))
 SA<-rbind(SA,data.table(db=k,a="Inverse probability of selection weighting (missing 48–96 h ALC)",n=nrow(IMP[[k]][[1]]),ev=sum(IMP[[k]][[1]]$d28),v=fm(q),p=q[4]))}
# Firth penalized logistic regression where events per variable were low
for(k in dbs[1:2]){for(sb in c("sepsis==0","ster==0","imv==1")){cv<-sub(paste0("\\+",sub("==.*","",sb)),"",M3); r<-sapply(IMP[[k]],function(d){d<-d[eval(parse(text=sb))]; f<-logistf(as.formula(paste("d28~x+",cv)),d,plcontrol=logistpl.control(maxit=200),pl=FALSE); c(coef(f)["x"],sqrt(diag(vcov(f)))["x"])}); q<-ex(rub(r[1,],r[2,])); d1<-IMP[[k]][[1]][eval(parse(text=sb))]
 SA<-rbind(SA,data.table(db=k,a=paste0("Firth penalized regression: ",c("sepsis==0"="excluding sepsis","ster==0"="excluding corticosteroids","imv==1"="invasive ventilation only")[sb]),n=nrow(d1),ev=sum(d1$d28),v=fm(q),p=q[4]))}}
R$t4a<-SA
# events per variable
np<-length(strsplit(M3,"\\+")[[1]])+1
R$epv<-rbindlist(lapply(dbs,function(k) data.table(db=k,events=sum(IMP[[k]][[1]]$d28),params=np+(k=="Pooled"),epv=sum(IMP[[k]][[1]]$d28)/(np+(k=="Pooled")))))
# ---------------- Subgroups (pooled, MICE) ----------------
sg<-list(old="Age ≥65 years",male="Male sex",sepsis="Sepsis",vaso="Vasopressor use (0–24 h)",imv="Invasive ventilation (0–24 h)",ster="Systemic corticosteroids (0–96 h)",copd="COPD",cancer="Solid cancer",sevlym="Baseline ALC <0.5 ×10⁹/L",hisev="Severity above database mean",DB="Database: eICU-CRD")
cvr<-paste(M3,"+DB")
R$t4b<-rbindlist(lapply(names(sg),function(v){cv<-gsub(paste0("\\+",v,"(?![a-z])|^",v,"\\+"),"",cvr,perl=TRUE); if(v=="old") cv<-sub("^age\\+","",cv); if(v=="sevlym") cv<-sub("\\+lalc0","",cv); if(v=="hisev") cv<-sub("\\+sevz","",cv)
 pi<-plg("Pooled",paste("d28~x*",v,"+",cv),paste0("^x:",v,"$"))[4]
 rbindlist(lapply(1:0,function(l){r<-plg("Pooled",paste("d28~x+",cv),"^x$",paste0(v,"==",l)); z<-IMP$Pooled[[1]][get(v)==l]; data.table(sg=sg[[v]],lev=ifelse(l==1,"Yes","No"),n=nrow(z),ev=sum(z$d28),v=fm(r),est=r[1],lo=r[2],hi=r[3],p=r[4],pint=pi)}))}))
# ---------------- Model diagnostics (first imputed data set) ----------------
R$diag<-rbindlist(lapply(dbs,function(k){d<-IMP[[k]][[1]]; cv<-cvp(M3,d); st<-if(k=="Pooled") "+strata(DB)" else ""
 z<-cox.zph(coxph(as.formula(paste("Surv(t28,d28)~x+",cv,st)),d)); g<-glm(as.formula(paste("d28~x+",cv)),binomial,d); vf<-vif(g); g2<-glm(as.formula(paste("d28~x+",cv,"+lanc")),binomial,d); vf2<-vif(g2)
 data.table(db=k,a=c("Schoenfeld residuals: persistent lymphopenia","Schoenfeld residuals: global test","Maximum VIF, model 3","VIF: WBC / log baseline ALC (model 3)","VIF: WBC / log baseline ALC / log neutrophils (model 3 + neutrophils)"),
  v=c("—","—",sprintf("%.2f",max(vf)),sprintf("%.2f / %.2f",vf["wbc"],vf["lalc0"]),sprintf("%.2f / %.2f / %.2f",vf2["wbc"],vf2["lalc0"],vf2["lanc"])),p=c(z$table["x","p"],z$table["GLOBAL","p"],NA,NA,NA))}))
# ---------------- Selection comparison: included vs excluded (no ALC at 48-96 h) ----------------
S0<-as.data.table(readRDS("S.rds")); sv<-c("age","male","sepsis","copd","cancer","ards","alc0","lact","creat","pf","sev","vaso","imv","ster","los","d28","hd")
R$sel<-rbindlist(lapply(dbs[1:2],function(k) rbindlist(lapply(1:0,function(g){z<-S0[db==k&inc==g]; data.table(db=k,inc=g,n=nrow(z),t(sapply(sv,function(v){u<-z[[v]]; if(all(u%in%c(0,1,NA))) sprintf("%.1f%%",100*mean(u,na.rm=TRUE)) else sprintf("%.2f (%.2f–%.2f)",median(u,na.rm=TRUE),quantile(u,.25,na.rm=TRUE),quantile(u,.75,na.rm=TRUE))})))}))))
R$selp<-sapply(sv,function(v) sapply(dbs[1:2],function(k){z<-S0[db==k]; if(all(is.na(z[[v]]))) NA else if(all(z[[v]]%in%c(0,1,NA))) chisq.test(table(z[[v]],z$inc))$p.value else wilcox.test(z[[v]]~z$inc)$p.value}))
saveRDS(R,"R.rds"); saveRDS(IMP,"IMP.rds"); cat("ANALYSIS DONE\n")
print(R$t1a); print(R$std); print(R$t1b); print(R$t1c); print(R$t2); print(R$t3[,.(db,a,v,p=signif(p,2))]); print(R$meta); print(R$t4a[,.(db,a,n,ev,v,p=signif(p,2))]); print(R$epv); print(R$t4b[,.(sg,lev,n,v,p=signif(p,2),pint=signif(pint,2))]); print(R$diag); print(R$sel)
