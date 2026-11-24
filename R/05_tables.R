# 05 tables: manuscript Tables 1-3 and supplementary tables (saved as T.rds for the manuscript build)
suppressMessages({library(data.table)})
R<-readRDS("R.rds"); T<-list()
pf<-function(p) ifelse(is.na(p),"",ifelse(p<0.001,"<0.001",sprintf("%.3f",p)))
fi<-function(v,p) ifelse(is.na(p),v,paste0(v,"; P ",ifelse(p<0.001,"< 0.001",paste0("= ",sprintf("%.3f",p)))))
rf<-function(z,dg) sapply(z,function(u){if(grepl("%",u)||grepl("NA",u)) return(u); x<-as.numeric(regmatches(u,gregexpr("-?[0-9]+\\.?[0-9]*",u))[[1]]); sprintf(paste0("%.",dg,"f (%.",dg,"f–%.",dg,"f)"),x[1],x[2],x[3])})
dbs<-c("MIMIC-IV","eICU-CRD","Pooled")
# ---- Table 1
vl<-c(n="Patients, n",age="Age, years",male="Male sex",chf="Heart failure",copd="COPD",ckd="Chronic kidney disease",dm="Diabetes mellitus",liver="Liver disease",cancer="Solid cancer",sepsis="Sepsis",ards="ARDS",asp="Aspiration pneumonia",
 alc0="Baseline ALC, × 10⁹/L",alc3="ALC at 48–96 h, × 10⁹/L",anc0="Baseline neutrophil count, × 10⁹/L",wbc="White blood cell count, × 10⁹/L",lact="Lactate, mmol/L",creat="Creatinine, mg/dL",bun="Blood urea nitrogen, mg/dL",hb="Haemoglobin, g/dL",plt="Platelet count, × 10⁹/L",alb="Albumin, g/dL",
 pf="Lowest PaO₂/FiO₂, mmHg",hr="Heart rate, beats/min",map="Mean arterial pressure, mmHg",rr="Respiratory rate, breaths/min",sev="Severity score*",vaso="Vasopressor use, 0–24 h",imv="Invasive ventilation, 0–24 h",ster="Systemic corticosteroids, 0–96 h",nalc="ALC measurements, 0–96 h, n",d28="28-day mortality",hd="In-hospital mortality")
t0<-copy(R$t0); t0[,n:=format(n,big.mark=",")]
for(v in c("age","plt","sev","nalc","pf","hr","map","rr","bun")) t0[,(v):=rf(get(v),0)]; for(v in c("anc0","wbc","lact","creat","hb","alb")) t0[,(v):=rf(get(v),1)]
w<-melt(t0[,c("db","grp",names(vl)),with=FALSE],id.vars=c("db","grp"),variable.factor=FALSE,value.name="v"); w<-dcast(w,variable~db+grp,value.var="v")
w[,`P (MIMIC-IV)`:=sapply(variable,function(v) if(v=="n") "" else pf(R$t0p["MIMIC-IV",v]))]; w[,`P (eICU-CRD)`:=sapply(variable,function(v) if(v=="n") "" else pf(R$t0p["eICU-CRD",v]))]
w<-w[match(names(vl),variable)]; w[,Characteristic:=vl[variable]]
T$t1<-w[,.(Characteristic,`MIMIC-IV_0`,`MIMIC-IV_1`,`P (MIMIC-IV)`,`eICU-CRD_0`,`eICU-CRD_1`,`P (eICU-CRD)`)]
# ---- Table 2
a<-R$t1a; s<-R$std; b<-R$t1b; OT<-R$t2
g<-function(k,col) a[db==k&grp==1][[col]]
ev<-function(k,gr){z<-a[db==k&grp==gr]; sprintf("%d/%s (%.1f%%)",z$ev,format(z$n,big.mark=","),100*z$ev/z$n)}
row<-function(lab,f) data.table(Outcome=lab,`MIMIC-IV`=f("MIMIC-IV"),`eICU-CRD`=f("eICU-CRD"),Pooled=f("Pooled"))
sec<-function(lab,on,mm=NULL) row(lab,function(k){z<-OT[OT$db==k&OT$o==on]; if(!is.null(mm)) z<-z[grepl(mm,z$m,fixed=TRUE)]; if(nrow(z)==0) "NA" else fi(z$v[1],z$p[1])})
mt<-R$meta
T$t2<-rbind(row("28-day mortality (primary outcome)",function(k) ""),
 row("  Deaths, recovered group",function(k) ev(k,0)),row("  Deaths, persistent lymphopenia group",function(k) ev(k,1)),
 row("  Model 1 (unadjusted), OR",function(k) g(k,"crude")),row("  Model 2, OR",function(k) g(k,"m2")),row("  Model 3, OR (primary analysis)",function(k) fi(g(k,"m3"),g(k,"p3"))),row("  Model 3, Cox HR",function(k) fi(g(k,"hr"),g(k,"ph"))),
 row("  Two-stage meta-analysis of database estimates, OR (fixed / random effects)",function(k) if(k=="Pooled") paste0(fi(mt$v[1],mt$p[1])," / ",sub(" \\(.*","",mt$v[2]),"; I² = ",mt$v[3],"%") else "NA"),
 row("  Adjusted mortality, recovered vs persistent, %",function(k){z<-s[db==k]; sprintf("%.1f vs %.1f",100*z$p0,100*z$p1)}),row("  Adjusted risk difference, percentage points",function(k){z<-s[db==k]; sprintf("%.1f (%.1f–%.1f)",100*z$rd,100*z$rdlo,100*z$rdhi)}),
 row("  Per 0.1 × 10⁹/L lower ALC at 48–96 h, Model 3 OR",function(k){z<-b[db==k][1]; fi(z$v,z$p)}),row("  P for trend across ALC categories at 48–96 h",function(k) pf(b[db==k&grepl("trend",a)]$p)),
 row("Secondary outcomes (Model 3)",function(k) ""),sec("  In-hospital mortality, OR","In-hospital mortality"),sec("  90-day mortality, HR","90-day mortality"),sec("  1-year mortality, HR","1-year mortality"),
 sec("  Ventilator-free days to day 28, β (days)","Ventilator-free days to day 28"),sec("  ICU length of stay, β (days)","ICU length of stay"),sec("  Bloodstream infection, 96 h to day 28, OR","Bloodstream infection (96 h–day 28)","Logistic"),sec("  Bloodstream infection, Fine–Gray SHR","Bloodstream infection (96 h–day 28)","Fine"))
# ---- Table 3
cc<-R$t3; sa<-R$t4a
c3<-function(lab,pat,src) row(lab,function(k){z<-src[db==k&grepl(pat,a,fixed=TRUE)]; if(!nrow(z)) "NA" else fi(z$v[1],z$p[1])})
T$t3<-rbind(row("Confounding control (multiply imputed data, Rubin's rules)",function(k) ""),c3("  IPTW (ATE), OR","IPTW (ATE",cc),c3("  IPTW, Cox HR","IPTW, HR",cc),c3("  Doubly robust (IPTW + Model 3), OR","Doubly robust",cc),c3("  Overlap weights (ATO), OR","Overlap weights",cc),
 row("  PSM 1:1, OR (median matched pairs)",function(k){z<-cc[db==k&grepl("^PSM 1:1",a)]; paste0(fi(z$v,z$p)," (",sub(".*median (\\d+) pairs.*","\\1",z$a),")")}),c3("  PSM + Model 3, OR","PSM + Model 3",cc),c3("  Hospital-clustered robust SE, OR","hospital-clustered",cc),
 row("  E-value (point estimate; CI limit)",function(k) cc[db==k&grepl("E-value",a)]$v),row("  Maximum absolute SMD, unweighted / IPTW",function(k){z<-cc[db==k&grepl("SMD",a)]$a; m<-regmatches(z,gregexpr("[0-9]\\.[0-9]+",z))[[1]]; paste(m[1],"/",m[2])}),
 row("Sensitivity analyses (Model 3 OR)",function(k) ""),
 c3("  Persistence threshold <0.8 × 10⁹/L","Persistent threshold <0.8",sa),c3("  Lowest rather than last ALC at 48–96 h","lowest (not last)",sa),c3("  Exposure window 48–72 h","window 48–72",sa),c3("  Adjusted for number of ALC measurements","measurement intensity",sa),
 c3("  Adjusted for baseline neutrophil count","neutrophil",sa),c3("  Adjusted for albumin and bilirubin","albumin and bilirubin",sa),c3("  Pneumonia as principal diagnosis†","principal diagnosis",sa),c3("  Excluding aspiration pneumonia","aspiration",sa),
 c3("  Excluding sepsis","Excluding sepsis",sa),c3("  Excluding solid cancer","solid cancer",sa),c3("  Excluding corticosteroid use, 0–96 h","Excluding systemic corticosteroids",sa),c3("  Excluding death before day 7","before day 7",sa),
 c3("  Invasive ventilation within 24 h","Invasive ventilation 0–24 h only",sa),c3("  No invasive ventilation within 24 h","No invasive ventilation",sa),c3("  Selection weighting for missing ALC at 48–96 h","selection weighting",sa),
 c3("  Single median imputation","Single median imputation",sa),c3("  Complete-case analysis","Complete-case",sa),
 c3("  Firth regression, excluding sepsis","Firth penalized regression: excluding sepsis",sa),c3("  Firth regression, excluding corticosteroid use","Firth penalized regression: excluding corticosteroids",sa),c3("  Firth regression, invasive ventilation only","Firth penalized regression: invasive ventilation only",sa))
# ---- Supplementary
T$miss<-R$miss
T$s3<-rbind(b[,.(Database=db,Analysis=a,Estimate=v,P=pf(p))],R$t1c[,.(Database=db,Analysis=a,Estimate=v,P=pf(p))])
T$s4<-R$t4b[,.(Subgroup=sg,Level=lev,`Deaths/N`=paste0(ev,"/",format(n,big.mark=",")),`OR (95% CI)`=v,P=pf(p),`P for interaction`=pf(pint))]
sl<-c(n="Patients, n",age="Age, years",male="Male sex",sepsis="Sepsis",copd="COPD",cancer="Solid cancer",ards="ARDS",alc0="Baseline ALC, × 10⁹/L",lact="Lactate, mmol/L",creat="Creatinine, mg/dL",pf="Lowest PaO₂/FiO₂, mmHg",sev="Severity score*",vaso="Vasopressor use, 0–24 h",imv="Invasive ventilation, 0–24 h",ster="Systemic corticosteroids, 0–96 h",los="ICU length of stay, days",d28="28-day mortality",hd="In-hospital mortality")
s0<-copy(R$sel); s0[,n:=format(n,big.mark=",")]; for(v in c("age","sev","pf")) s0[,(v):=rf(get(v),0)]; for(v in c("lact","creat","los")) s0[,(v):=rf(get(v),1)]
w<-melt(s0[,c("db","inc",names(sl)),with=FALSE],id.vars=c("db","inc"),variable.factor=FALSE,value.name="v"); w<-dcast(w,variable~db+inc,value.var="v")
w[,`P (MIMIC-IV)`:=sapply(variable,function(v) if(v=="n") "" else pf(R$selp["MIMIC-IV",v]))]; w[,`P (eICU-CRD)`:=sapply(variable,function(v) if(v=="n") "" else pf(R$selp["eICU-CRD",v]))]
w<-w[match(names(sl),variable)]; w[,Characteristic:=sl[variable]]
T$s5<-w[,.(Characteristic,`MIMIC-IV_1`,`MIMIC-IV_0`,`P (MIMIC-IV)`,`eICU-CRD_1`,`eICU-CRD_0`,`P (eICU-CRD)`)]
T$s6<-rbind(R$diag[,.(Database=db,Check=a,Value=v,P=pf(p))],R$epv[,.(Database=db,Check="Events per model parameter (model 3)",Value=sprintf("%d events / %d parameters = %.1f",events,params,epv),P="")])
saveRDS(T,"T.rds"); print(T$t2); print(T$t3); cat("tables done\n")
