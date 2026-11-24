# 04 figures: Figure 1 flow chart, Figure 2 Kaplan-Meier, Supplementary Figure S1 restricted cubic spline, Supplementary Figure S2 subgroup forest plot (TIFF 300 dpi)
suppressMessages({library(data.table);library(ggplot2);library(splines);library(survival);library(patchwork)})
R<-readRDS("R.rds"); IMP<-readRDS("IMP.rds"); D<-as.data.table(readRDS("D.rds")); fl<-readRDS("fl.rds"); FT<-"Times New Roman"; dir.create("fig",showWarnings=FALSE)
sv<-function(p,f,w,h){ggsave(file.path("fig",paste0(f,".tiff")),p,width=w,height=h,dpi=300,compression="none",bg="white"); ggsave(file.path("fig",paste0(f,".png")),p,width=w,height=h,dpi=150,bg="white")}
fm<-function(n) format(n,big.mark=","); pv<-function(p) if(p<0.001) "P < 0.001" else sprintf("P = %.3f",p)
## ---- Figure 1: flow chart (arrow heads stop short of the boxes)
bx<-function(x,y,w,h,lab,id) data.table(id=id,x=x,y=y,w=w,h=h,lab=lab)
mk<-function(x0,db,main,exc){B1<-list(); E1<-list(); n<-length(main); ys<-seq(10,10-1.9*(n-1),by=-1.9)
 for(i in 1:n) B1[[i]]<-bx(x0,ys[i],3.3,1.0,main[i],paste0(db,"m",i))
 for(i in seq_along(exc)) if(!is.na(exc[i])) E1[[length(E1)+1]]<-bx(x0+3.6,ys[i]-0.95,2.95,0.85,exc[i],paste0(db,"e",i))
 list(b=rbindlist(B1),e=rbindlist(E1))}
cn<-D[,.N,by=.(db,x)]
mm<-function(f,p,lab) c(sprintf("%s\nFirst ICU stays of adults\nn = %s",lab,fm(f[[paste0(p,"0")]])),sprintf("Pneumonia\nn = %s",fm(f[[paste0(p,"_pna")]])),sprintf("No haematological malignancy,\nHIV/AIDS or transplantation\nn = %s",fm(f[[paste0(p,"_imm")]])),sprintf("Alive and in hospital\nat the 96-h landmark\nn = %s",fm(f[[paste0(p,"_lm")]])),sprintf("Baseline lymphopenia\n(first ALC −6 to 24 h <1.0 × 10⁹/L)\nn = %s",fm(f[[paste0(p,"_lym")]])),sprintf("Analysed cohort\n(≥1 ALC at 48–96 h)\nn = %s",fm(f[[paste0(p,"_b3")]])))
me<-function(f,p) c(NA,sprintf("Haematological malignancy,\nHIV/AIDS or transplantation:\nn = %s",fm(f[[paste0(p,"_pna")]]-f[[paste0(p,"_imm")]])),sprintf("Death or discharge\nbefore 96 h: n = %s",fm(f[[paste0(p,"_imm")]]-f[[paste0(p,"_lm")]])),sprintf("No baseline ALC n = %s;\nbaseline ALC ≥1.0 n = %s",fm(f[[paste0(p,"_lm")]]-f[[paste0(p,"_b0")]]),fm(f[[paste0(p,"_b0")]]-f[[paste0(p,"_lym")]])),sprintf("No ALC at 48–96 h:\nn = %s",fm(f[[paste0(p,"_lym")]]-f[[paste0(p,"_b3")]])),NA)
M<-mk(0,"M",mm(fl,"m","MIMIC-IV v3.1"),me(fl,"m")); E<-mk(8.4,"E",mm(fl,"e","eICU-CRD v2.0"),me(fl,"e"))
yg<-10-1.9*5-1.75
gb<-rbind(bx(-0.95,yg,1.75,0.95,sprintf("Recovered\nn = %s",fm(cn[db=="MIMIC-IV"&x==0]$N)),"Mg0"),bx(0.95,yg,1.75,0.95,sprintf("Persistent\nn = %s",fm(cn[db=="MIMIC-IV"&x==1]$N)),"Mg1"),
 bx(7.45,yg,1.75,0.95,sprintf("Recovered\nn = %s",fm(cn[db=="eICU-CRD"&x==0]$N)),"Eg0"),bx(9.35,yg,1.75,0.95,sprintf("Persistent\nn = %s",fm(cn[db=="eICU-CRD"&x==1]$N)),"Eg1"))
BX<-rbind(M$b,M$e,E$b,E$e,gb); g<-0.13; ar<-list()
vs<-function(b) for(i in 1:(nrow(b)-1)) ar[[length(ar)+1]]<<-data.table(x=b$x[i],y=b$y[i]-b$h[i]/2,xe=b$x[i+1],ye=b$y[i+1]+b$h[i+1]/2+g)
hz<-function(b,e) for(i in 1:nrow(e)) ar[[length(ar)+1]]<<-data.table(x=b$x[1],y=e$y[i],xe=e$x[i]-e$w[i]/2-g,ye=e$y[i])
sp<-function(b,gx){lb<-b[.N]; for(i in 1:nrow(gx)) ar[[length(ar)+1]]<<-data.table(x=lb$x,y=lb$y-lb$h/2,xe=gx$x[i],ye=gx$y[i]+gx$h[i]/2+g)}
vs(M$b); vs(E$b); hz(M$b,M$e); hz(E$b,E$e); sp(M$b,gb[1:2]); sp(E$b,gb[3:4]); AR<-rbindlist(ar)
p1<-ggplot()+geom_rect(data=BX,aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),fill="white",colour="black",linewidth=.4)+geom_text(data=BX,aes(x,y,label=lab),family=FT,size=2.55,lineheight=.92)+
 geom_segment(data=AR,aes(x=x,y=y,xend=xe,yend=ye),arrow=arrow(length=unit(1.5,"mm"),type="closed"),linewidth=.35)+coord_cartesian(xlim=c(-2,13.6),ylim=c(yg-0.7,10.7),expand=FALSE)+theme_void(base_family=FT)+theme(plot.background=element_rect(fill="white",colour=NA))
sv(p1,"Figure1",8.4,7.4)
## ---- Figure 2: Kaplan-Meier from the 96-h landmark with numbers at risk
D[,grp:=factor(ifelse(x==1,"Persistent lymphopenia","Recovered"),levels=c("Recovered","Persistent lymphopenia"))]
kmp<-function(k,tag){d<-D[db==k]; f<-survfit(Surv(t28/24+4,d28)~grp,d); s<-summary(f,times=seq(4,28,0.25),extend=TRUE)
 z<-data.table(t=s$time,cum=1-s$surv,lo=1-s$upper,hi=1-s$lower,g=factor(sub("grp=","",s$strata),levels=c("Recovered","Persistent lymphopenia"))); r<-summary(f,times=c(4,7,14,21,28),extend=TRUE); rk<-data.table(t=r$time,n=r$n.risk,g=factor(sub("grp=","",r$strata),levels=c("Persistent lymphopenia","Recovered")))
 pl<-R$t1c[db==k]$p[1]; hr<-R$t1a[db==k&grp==1]
 a<-ggplot(z,aes(t,cum,colour=g))+geom_step(linewidth=.7)+scale_colour_manual(values=c("#4D4D4D","#B2182B"),name=NULL)+scale_x_continuous(breaks=c(4,7,14,21,28),limits=c(3.5,28.5))+scale_y_continuous(labels=function(v) paste0(round(100*v),"%"),limits=c(0,NA))+
  labs(title=paste0(tag," ",k),x=NULL,y="Cumulative mortality")+annotate("text",x=4.3,y=max(z$cum)*0.97,hjust=0,family=FT,size=3.1,label=sprintf("Log-rank %s\nAdjusted HR %s",pv(pl),hr$hr))+
  theme_classic(base_family=FT,base_size=10.5)+theme(legend.position=c(.75,.18),plot.title=element_text(face="bold",size=11))
 b<-ggplot(rk,aes(t,g,label=fm(n)))+geom_text(family=FT,size=3)+scale_x_continuous(breaks=c(4,7,14,21,28),limits=c(3.5,28.5))+labs(x="Days since ICU admission",y=NULL,subtitle="Number at risk")+theme_classic(base_family=FT,base_size=10)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.subtitle=element_text(size=9))
 a/b+plot_layout(heights=c(4,1))}
p2<-(kmp("MIMIC-IV","A")|kmp("eICU-CRD","B")); sv(p2,"Figure2",10,4.8)
## ---- Figure 3: restricted cubic spline, model 3 (MIMIC-IV, eICU-CRD, pooled)
M3<-"age+male+chf+ckd+copd+dm+liver+cancer+sepsis+sevz+lalc0+lact+creat+bun+hb+wbc+plt+map+hr+rr+lpf+vaso+imv+ster+ards"
# spline curve pooled across imputed data sets (Rubin's rules for the spline coefficients)
rcsb<-function(x,k){nk<-length(k); nr<-(k[nk]-k[1])^2; cbind(x,matrix(sapply(1:(nk-2),function(j) (pmax(x-k[j],0)^3-pmax(x-k[nk-1],0)^3*(k[nk]-k[j])/(k[nk]-k[nk-1])+pmax(x-k[nk],0)^3*(k[nk-1]-k[j])/(k[nk]-k[nk-1]))/nr),nrow=length(x)))}
rcs1<-function(L,lab,cv){d0<-L[[1]]; lo<-quantile(d0$alc3,.01); hi<-quantile(d0$alc3,.99); k4<-quantile(pmin(pmax(d0$alc3,lo),hi),c(.05,.35,.65,.95))
 fs<-lapply(L,function(d){d<-copy(d); B<-rcsb(pmin(pmax(d$alc3,lo),hi),k4); d[,`:=`(r1=B[,1],r2=B[,2],r3=B[,3])]; f<-glm(as.formula(paste("d28~r1+r2+r3+",cv)),binomial,d); list(b=coef(f)[c("r1","r2","r3")],V=vcov(f)[c("r1","r2","r3"),c("r1","r2","r3")])})
 Q<-sapply(fs,`[[`,"b"); qb<-rowMeans(Q); U<-Reduce(`+`,lapply(fs,`[[`,"V"))/length(fs); Tv<-U+(1+1/length(fs))*cov(t(Q))
 g<-seq(quantile(d0$alc3,.02),quantile(d0$alc3,.98),length=200); Lm<-sweep(rcsb(g,k4),2,as.vector(rcsb(1,k4))); lp<-as.vector(Lm%*%qb); se<-sqrt(rowSums((Lm%*%Tv)*Lm))
 data.table(lab=lab,g,or=exp(lp),lo=exp(lp-1.96*se),hi=exp(lp+1.96*se))}
rc<-rbind(rcs1(IMP[["MIMIC-IV"]],"A MIMIC-IV",M3),rcs1(IMP[["eICU-CRD"]],"B eICU-CRD",M3),rcs1(IMP$Pooled,"C Pooled (database indicator)",paste(M3,"+DB")))
pt<-R$t1c[grepl("RCS",a)]; lab<-sapply(c("MIMIC-IV","eICU-CRD","Pooled"),function(k) sprintf("P overall %s\nP non-linear %s",sub("P ","",pv(pt[db==k]$p[1])),sub("P ","",pv(pt[db==k]$p[2])))); names(lab)<-unique(rc$lab)
an<-data.table(lab=names(lab),txt=lab)
p3<-ggplot(rc,aes(g,or))+geom_ribbon(aes(ymin=lo,ymax=hi),fill="grey85")+geom_line(colour="#B2182B",linewidth=.8)+geom_hline(yintercept=1,linetype=2)+geom_vline(xintercept=1,linetype=3,colour="grey40")+facet_wrap(~lab,nrow=1)+
 geom_text(data=an,aes(x=-Inf,y=Inf,label=txt),hjust=-0.08,vjust=1.3,family=FT,size=3)+scale_y_log10()+labs(x="Absolute lymphocyte count at 48–96 h (× 10⁹/L; reference 1.0)",y="Adjusted OR for 28-day mortality (log scale)")+
 theme_classic(base_family=FT,base_size=10.5)+theme(strip.background=element_blank(),strip.text=element_text(face="bold",hjust=0))
sv(p3,"FigureS1",10,3.8)
## ---- Figure 4: subgroup forest plot (pooled)
sg<-copy(R$t4b); sg[,lev2:=ifelse(lev=="Yes","Yes","No")]
nm<-c("Age ≥65 years"="Age ≥65 years","Male sex"="Male sex","Sepsis"="Sepsis","Vasopressor use (0–24 h)"="Vasopressor use, 0–24 h","Invasive ventilation (0–24 h)"="Invasive ventilation, 0–24 h","Systemic corticosteroids (0–96 h)"="Systemic corticosteroids, 0–96 h","COPD"="COPD","Solid cancer"="Solid cancer","Baseline ALC <0.5 ×10⁹/L"="Baseline ALC <0.5 × 10⁹/L","Severity above database mean"="Severity score above database mean","Database: eICU-CRD"="eICU-CRD (vs MIMIC-IV)")
sg[,sg2:=nm[sg]]; rows<-list(); for(s in unique(sg$sg2)){z<-sg[sg2==s]; rows[[length(rows)+1]]<-data.table(lab=s,est=NA,lo=NA,hi=NA,n="",ci="",pint=sprintf("%.3f",z$pint[1]),b=TRUE); for(l in c("Yes","No")){w<-z[lev2==l]; rows[[length(rows)+1]]<-data.table(lab=paste0("    ",l),est=w$est,lo=w$lo,hi=w$hi,n=paste0(fm(w$ev),"/",fm(w$n)),ci=sprintf("%.2f (%.2f–%.2f)",w$est,w$lo,w$hi),pint="",b=FALSE)}}
FR<-rbindlist(rows); FR[,y:=rev(seq_len(.N))]
th<-theme_void(base_family=FT,base_size=10)
ylim<-c(0.5,nrow(FR)+1.2); hd<-function(t,x) annotate("text",x=x,y=nrow(FR)+1,label=t,family=FT,fontface="bold",size=3.1,hjust=0)
pa<-ggplot(FR)+geom_text(aes(0,y,label=lab,fontface=ifelse(b,"bold","plain")),family=FT,size=3,hjust=0)+geom_text(aes(1.45,y,label=n),family=FT,size=3,hjust=0)+hd("Subgroup",0)+hd("Deaths/N",1.45)+coord_cartesian(xlim=c(0,2.1),ylim=ylim,expand=FALSE)+th
pb<-ggplot(FR[!is.na(est)])+geom_vline(xintercept=1,linetype=2)+geom_errorbar(aes(y=y,xmin=lo,xmax=hi),width=0.25,orientation="y")+geom_point(aes(est,y),shape=15,size=2)+scale_x_log10(breaks=c(0.5,1,2,4))+coord_cartesian(xlim=c(0.7,4.5),ylim=ylim,expand=FALSE)+labs(x="Adjusted OR (log scale)")+theme_classic(base_family=FT,base_size=10)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
pc<-ggplot(FR)+geom_text(aes(0,y,label=ci),family=FT,size=3,hjust=0)+geom_text(aes(1.25,y,label=pint),family=FT,size=3,hjust=0)+hd("OR (95% CI)",0)+hd("P for interaction",1.25)+coord_cartesian(xlim=c(0,2.1),ylim=ylim,expand=FALSE)+th
p4<-pa+pb+pc+plot_layout(widths=c(1.25,1,1.1)); sv(p4,"FigureS2",10,7.2)
cat("figures done\n")
