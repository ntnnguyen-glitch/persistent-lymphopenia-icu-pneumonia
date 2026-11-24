# 01 cohort: adults with pneumonia admitted to the ICU, lymphopenic at baseline (ALC <1.0 x10^9/L, first value -6 to 24 h),
# alive and in hospital at the 96-h landmark with >=1 ALC during 48-96 h (patients without a 48-96 h ALC are kept in S.rds
# for the selection-bias comparison).
# Exposure: persistent lymphopenia = last ALC in 48-96 h <1.0 x10^9/L (reference: recovered, >=1.0)
suppressMessages({library(DBI);library(data.table)})
con<-dbConnect(duckdb::duckdb(shared_home=FALSE),":memory:"); invisible(dbExecute(con,sprintf("ATTACH '%s' AS w (READ_ONLY)",Sys.getenv("DUCKDB","w.duckdb"))))
Q<-function(s) as.data.table(dbGetQuery(con,s)); H<-function(a,b) as.numeric(difftime(a,b,units="hours")); fl<-list(); LM<-96
# ================= MIMIC-IV =================
m<-Q("select * from w.st"); fl$m0<-nrow(m)
dx<-Q("select hadm_id,icd_code from w.dx"); fx<-function(p) as.integer(m$hadm_id%in%dx[grepl(p,icd_code)]$hadm_id)
m[,pna:=fx("^(48[0-6]|J1[2-8]|J69|5070)")]; m<-m[pna==1]; fl$m_pna<-nrow(m)
m[,`:=`(hem=fx("^(20[0-8]|C8[1-9]|C9[0-6])"),hiv=fx("^(042|B2[0-4])"),tx=fx("^(V42|Z94|99681|99682|99685|99687|T86)"))]
m<-m[hem==0&hiv==0&tx==0]; fl$m_imm<-nrow(m)
dth<-as.POSIXct(m$deathtime,tz="UTC"); dd<-as.POSIXct(ifelse(is.na(m$dod),NA,paste(as.character(m$dod),"23:59:00")),tz="UTC"); dth[is.na(dth)]<-dd[is.na(dth)]; m[,death_h:=pmax(H(dth,intime),0)]
m[,disch_h:=H(dischtime,intime)]; m<-m[(is.na(death_h)|death_h>LM)&disch_h>LM]; fl$m_lm<-nrow(m)
ml<-as.data.table(readRDS("l_mlab.rds")); ml<-merge(ml[hadm_id%in%m$hadm_id&!is.na(v)],m[,.(hadm_id,stay_id,intime)],by="hadm_id"); ml[,h:=H(charttime,intime)]; ml<-ml[h>=-6&h<=168]
# ALC: reported absolute count (itemid 51133) or WBC x lymphocyte % from the same specimen (specimen_id)
ab<-ml[itemid==51133&v>0][,v:=fifelse(v>50,v/1000,v)][,.(stay_id,sp,h,alc=v)]
pr<-merge(ml[itemid==51244&v>=0&v<=100,.(stay_id,sp,h,lp=v)],unique(ml[itemid==51301&v>0&v<200,.(sp,w=v)],by="sp"),by="sp")[,.(stay_id,sp,h,alc=w*lp/100)]
A<-unique(rbind(ab,pr),by=c("stay_id","sp"))[alc<50]
b0<-A[h>=-6&h<=24][order(h)][,.(alc0=alc[1]),by=stay_id]; b3<-A[h>48&h<=LM][order(-h)][,.(alc3=alc[1],alcmin3=min(alc)),by=stay_id]; b72<-A[h>48&h<=72][order(-h)][,.(alc72=alc[1]),by=stay_id]
nm<-A[h>=-6&h<=LM,.(nalc=.N),by=stay_id]
m<-merge(m,b0,by="stay_id"); fl$m_b0<-nrow(m); m<-m[alc0<1]; fl$m_lym<-nrow(m); m<-merge(m,b3,by="stay_id",all.x=TRUE); fl$m_b3<-sum(!is.na(m$alc3))
m[,`:=`(alc72=b72$alc72[match(stay_id,b72$stay_id)],nalc=nm$nalc[match(stay_id,nm$stay_id)])]
nu<-merge(ml[itemid==51256&v>=0&v<=100&h>=-6&h<=24,.(stay_id,sp,h,np=v)],unique(ml[itemid==51301&v>0&v<200,.(sp,w=v)],by="sp"),by="sp")[order(h)][,.(anc=w[1]*np[1]/100),by=stay_id]
m[,anc0:=nu$anc[match(stay_id,nu$stay_id)]]
mdx<-as.data.table(readRDS("l_mdx.rds")); m[,pna3:=as.integer(hadm_id%in%mdx[!is.na(seq)&seq<=3&grepl("^(48[0-6]|J1[2-8]|J69|5070)",icd_code)]$hadm_id)]
m[,`:=`(chf=fx("^(428|I50)"),copd=fx("^(49[0-2]|496|J4[0-4])"),ckd=fx("^(585|N18)"),dm=fx("^(250|E1[0-4])"),liver=fx("^(571|K70|K7[2-6])"),cancer=fx("^(1[4-9][0-9]|19[0-9]|C[0-7][0-9]|C80)"),sepsis=fx("^(99591|99592|78552|A40|A41|R652)"),ards=fx("^(51882|J80)"),asp=fx("^(5070|J69)"),male=as.integer(gender=="M"))]
rx<-as.data.table(readRDS("l_mrx.rds")); rx<-merge(rx[hadm_id%in%m$hadm_id],m[,.(hadm_id,stay_id,intime)],by="hadm_id"); rx[,h:=H(starttime,intime)]; rx[,r:=toupper(route)]
rx<-rx[!grepl("TP|TOPICAL|OPHTH|INH|NEB|EAR|NU|IN|PR|RECTAL|INTRA-ARTICULAR|IA",r)|r%in%c("IV","IV PUSH","IV BOLUS","PO","PO/NG","NG","IV DRIP","IM")]
m[,`:=`(ster=as.integer(stay_id%in%rx[h>=-24&h<=LM]$stay_id),ster0=as.integer(stay_id%in%rx[h>=-24&h<=24]$stay_id))]
L<-Q("select stay_id,specimen_id,itemid,charttime,v,value from w.lab"); L<-merge(L[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); L[,h:=H(charttime,intime)]
gv<-function(it,lo,hi,fun="first",a=-6,z=24){x<-L[itemid%in%it&v>=lo&v<=hi&h>=a&h<=z][order(h)]; x<-if(fun=="first") x[,.(v=v[1]),by=stay_id] else x[,.(v=if(fun=="max") max(v) else min(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(lact=gv(50813,.1,30),bun=gv(51006,1,250),creat=gv(50912,.1,25),hb=gv(51222,3,20),wbc=gv(51301,.1,200),plt=gv(51265,1,1500),pltmin=gv(51265,1,1500,"min",0,24),alb=gv(50862,1,6),bili=gv(50885,.1,60,"max"),crmax=gv(50912,.1,25,"max",0,24))]
# lowest PaO2/FiO2, 0-24 h: arterial PaO2 with FiO2 from the same blood gas or charted FiO2 within the preceding 4 h
spx<-L[itemid==52033,.(specimen_id,st=value)]; g<-dcast(L[itemid%in%c(50821,50816)&h>=0&h<=24],stay_id+specimen_id+h~itemid,value.var="v",fun.aggregate=function(x)x[1]); setnames(g,c("50821","50816"),c("pao2","fl"),skip_absent=TRUE)
g[,st:=spx$st[match(specimen_id,spx$specimen_id)]]; g<-g[(is.na(st)|st=="ART.")&!is.na(pao2)&pao2>=20&pao2<=700]
ce<-Q("select stay_id,itemid,charttime,v from w.ce"); ce<-merge(ce[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); ce[,h:=H(charttime,intime)]
cm<-function(its,lo,hi,fun=mean,a=0,z=24){x<-ce[itemid%in%its&v>=lo&v<=hi&h>=a&h<=z,.(v=fun(v)),by=stay_id]; x$v[match(m$stay_id,x$stay_id)]}
m[,`:=`(hr=cm(220045,20,250),map=cm(c(220052,220181,225312),20,200),mapmin=cm(c(220052,220181,225312),20,200,min),rr=cm(c(220210,224690),3,70))]
ip<-as.data.table(readRDS("l_minp.rds")); ip<-merge(ip[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); ip[,h:=H(starttime,intime)]
m[,vaso:=as.integer(stay_id%in%ip[h>=-6&h<=24]$stay_id)]
pe<-Q("select stay_id,starttime,endtime from w.pe where itemid=225792"); pe<-merge(pe[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); pe[,`:=`(s=H(starttime,intime),e=H(endtime,intime))]
m[,imv:=as.integer(stay_id%in%pe[s>=-6&s<=24]$stay_id)]
mce<-as.data.table(readRDS("l_mce.rds")); mce<-merge(mce[stay_id%in%m$stay_id],m[,.(stay_id,intime)],by="stay_id"); mce[,h:=H(charttime,intime)]
fc<-mce[itemid==223835&!is.na(v)&v>=21&v<=100,.(stay_id,hf=h,f=v)]; pc2<-merge(g[,.(stay_id,specimen_id,ha=h)],fc,by="stay_id",allow.cartesian=TRUE)[hf<=ha&ha-hf<=4][order(ha-hf)][,.SD[1],by=specimen_id]
g[,fc:=pc2$f[match(specimen_id,pc2$specimen_id)]]; g[,f:=fifelse(!is.na(fl)&fl>=21&fl<=100,fl,fc)]; pfw<-g[!is.na(f),.(pf=min(pao2/(f/100))),by=stay_id]; m[,pf:=pfw$pf[match(stay_id,pfw$stay_id)]]
# invasive ventilation: procedure event 225792 or charted endotracheal/tracheostomy tube as oxygen delivery device, -6 to 24 h
ett<-mce[itemid==226732&value%in%c("Endotracheal tube","Tracheostomy tube")&h>=-6&h<=24]$stay_id
m[,imv:=as.integer(imv==1|stay_id%in%ett)]
gc<-dcast(mce[itemid%in%c(220739,223900,223901)&!is.na(v)&h>=0&h<=24],stay_id+charttime~itemid,value.var="v",fun.aggregate=function(x)x[1]); setnames(gc,c("220739","223900","223901"),c("e","vb","mo"),skip_absent=TRUE)
gcs<-gc[!is.na(e)&!is.na(vb)&!is.na(mo),.(gcs=min(e+vb+mo)),by=stay_id]; m[,gcs:=gcs$gcs[match(stay_id,gcs$stay_id)]]
# first-day SOFA (six organ systems); missing components scored 0; cardiovascular component approximated (no infusion doses)
sc<-function(x,br) fifelse(is.na(x),0L,as.integer(findInterval(x,br)))
m[,resp:=fifelse(is.na(pf),0L,fifelse(pf<100&imv==1,4L,fifelse(pf<200&imv==1,3L,fifelse(pf<300,2L,fifelse(pf<400,1L,0L)))))]
m[,cns:=fifelse(is.na(gcs),0L,fifelse(gcs<6,4L,fifelse(gcs<10,3L,fifelse(gcs<13,2L,fifelse(gcs<15,1L,0L)))))]
m[,sev:=as.numeric(resp+cns+sc(bili,c(1.2,2,6,12))+fifelse(vaso==1,3L,fifelse(!is.na(mapmin)&mapmin<70,1L,0L))+sc(crmax,c(1.2,2,3.5,5))+sc(-pltmin,-c(150,100,50,20)))]
v28<-pe[s<672,.(vh=sum(pmin(e,672)-pmax(s,0))),by=stay_id]; m[,vh:=v28$vh[match(stay_id,v28$stay_id)]]
mi<-as.data.table(readRDS("l_mmicro.rds")); mi<-mi[subject_id%in%m$subject_id&!is.na(org_name)&org_name!=""&!grepl("COAGULASE NEGATIVE|CORYNEBACTERIUM|BACILLUS|PROPIONIBACTERIUM|MICROCOCCUS|CANCELLED",org_name)]
mi[,t:=fifelse(is.na(charttime),chartdate,charttime)]; mi<-merge(mi,m[,.(subject_id,stay_id,intime)],by="subject_id"); mi[,h:=H(t,intime)]
m[,bsi:=as.integer(stay_id%in%mi[h>LM&h<=672]$stay_id)]
m[,`:=`(vfd=fifelse(!is.na(death_h)&death_h<=672,0,pmax(0,28-fifelse(is.na(vh),0,vh)/24)),los=H(outtime,intime)/24,hd=hospital_expire_flag,db="MIMIC-IV",id=stay_id,hosp=NA_character_)]
# ================= eICU-CRD (diagnosis strings in lower case) =================
e<-Q("select cast(patientunitstayid as bigint) s,uniquepid,gender,age,try_cast(hospitaldischargeoffset as int) hdo,hospitaldischargestatus hds,try_cast(unitdischargeoffset as int) udo,try_cast(unitvisitnumber as int) uvn,hospitaldischargeyear yr,hospitalid hosp from w.ept")
e<-e[uvn==1][order(yr,s)][,.SD[1],by=uniquepid]; e[,age:=as.numeric(ifelse(age=="> 89","90",age))]; e<-e[!is.na(age)&age>=18&hds%in%c("Alive","Expired")&!is.na(hdo)]; fl$e0<-nrow(e)
edx<-Q("select s,d from w.edx"); eph<-Q("select s,d from w.eph"); fe<-function(p) as.integer(e$s%in%c(edx[grepl(p,d,perl=TRUE)]$s,eph[grepl(p,d,perl=TRUE)]$s)); fd<-function(p) as.integer(e$s%in%edx[grepl(p,d,perl=TRUE)]$s)
e<-e[s%in%edx[grepl("pneumonia",d)]$s]; fl$e_pna<-nrow(e)
e[,`:=`(hem=fe("leukemia|lymphoma|myeloma|myelodysplas"),hiv=fe("\\bhiv\\b|\\baids\\b"),tx=fe("transplant"))]; e<-e[hem==0&hiv==0&tx==0]; fl$e_imm<-nrow(e)
e[,death_h:=fifelse(hds=="Expired",pmax(hdo/60,0),NA_real_)]; e[,disch_h:=hdo/60]; e<-e[(is.na(death_h)|death_h>LM)&disch_h>LM]; fl$e_lm<-nrow(e)
el<-as.data.table(readRDS("l_elab.rds")); el<-el[s%in%e$s&!is.na(v)&o>=-360&o<=10080]
A<-merge(el[n=="-lymphs"&v>=0&v<=100,.(s,o,lp=v)],unique(el[n=="WBC x 1000"&v>0&v<200,.(s,o,w=v)],by=c("s","o")),by=c("s","o")); A[,alc:=w*lp/100]
b0<-A[o>=-360&o<=1440][order(o)][,.(alc0=alc[1]),by=s]; b3<-A[o>2880&o<=LM*60][order(-o)][,.(alc3=alc[1],alcmin3=min(alc)),by=s]; b72<-A[o>2880&o<=4320][order(-o)][,.(alc72=alc[1]),by=s]; nm<-A[o>=-360&o<=LM*60,.(nalc=.N),by=s]
e<-merge(e,b0,by="s"); fl$e_b0<-nrow(e); e<-e[alc0<1]; fl$e_lym<-nrow(e); e<-merge(e,b3,by="s",all.x=TRUE); fl$e_b3<-sum(!is.na(e$alc3))
e[,`:=`(alc72=b72$alc72[match(s,b72$s)],nalc=nm$nalc[match(s,nm$s)])]
nu<-merge(el[n=="-polys"&v>=0&v<=100&o>=-360&o<=1440,.(s,o,np=v)],unique(el[n=="WBC x 1000"&v>0&v<200,.(s,o,w=v)],by=c("s","o")),by=c("s","o"))[order(o)][,.(anc=w[1]*np[1]/100),by=s]; e[,anc0:=nu$anc[match(s,nu$s)]]
e[,`:=`(chf=fe("heart failure"),copd=fe("copd|chronic obstructive|emphysema"),ckd=fe("chronic kidney|chronic renal|renal insufficiency"),dm=fe("diabetes"),liver=fe("cirrhosis|chronic liver"),cancer=fe("cancer|malignan|carcinoma|metasta"),sepsis=fd("sepsis|septic"),ards=fd("\\bards\\b|acute respiratory distress"),asp=fd("aspiration"),male=as.integer(gender=="Male"))]
em<-as.data.table(readRDS("l_emed.rds")); em[,r:=tolower(routeadmin)]; em<-em[s%in%e$s&!grepl("inhal|neb|topical|ophth|otic|rect|intra-art|cream|oint",r)&!grepl("cream|oint|inhal|drop|enema|foam|lotion",tolower(drugname))]
et<-as.data.table(readRDS("l_etr.rds")); et<-et[s%in%e$s&grepl("steroid|hydrocortisone|methylpred",d)]
st<-rbind(em[,.(s,o)],et[,.(s,o)]); e[,`:=`(ster=as.integer(s%in%st[o>=-1440&o<=LM*60]$s),ster0=as.integer(s%in%st[o>=-1440&o<=1440]$s))]
EL<-Q("select s,o,n,v from w.elab where o between -360 and 1440"); EL<-EL[s%in%e$s&!is.na(v)][order(o)]
ge<-function(nm,lo,hi,f="first"){x<-EL[n==nm&v>=lo&v<=hi]; x<-if(f=="first") x[,.(v=v[1]),by=s] else x[,.(v=max(v)),by=s]; x$v[match(e$s,x$s)]}
e[,`:=`(lact=ge("lactate",.1,30),bun=ge("BUN",1,250),creat=ge("creatinine",.1,25),hb=ge("Hgb",3,20),wbc=ge("WBC x 1000",.1,200),plt=ge("platelets x 1000",1,1500),alb=ge("albumin",1,6),bili=ge("total bilirubin",.1,60,"max"))]
aps<-Q("select cast(patientunitstayid as bigint) s,try_cast(heartrate as double) hr,try_cast(meanbp as double) map,try_cast(respiratoryrate as double) rr from w.eaps"); a<-aps[match(e$s,aps$s)]
e[,`:=`(hr=fifelse(a$hr>0,a$hr,NA_real_),map=fifelse(a$map>0,a$map,NA_real_),rr=fifelse(a$rr>0,a$rr,NA_real_))]
ap<-Q("select * from w.eapr"); e[,sev:=as.numeric(ap$apache[match(s,ap$s)])]; e[sev<0,sev:=NA]
pt<-Q("select cast(patientunitstayid as bigint) s,apacheadmissiondx ad from w.ept"); e[,pna3:=as.integer(grepl("^Pneumonia",pt$ad[match(s,pt$s)]))]
po<-EL[n=="paO2"&v>=20&v<=700&o>=0&o<=1440,.(s,o,pao2=v)]; r2<-Q("select s,o,l,v from w.erc2 where l in ('FiO2','FIO2 (%)') and o between -240 and 1440"); r2[,v:=suppressWarnings(as.numeric(v))]
fi2<-rbind(r2[!is.na(v)&v>=21&v<=100,.(s,of=o,f=v)],EL[n=="FiO2"&v>=21&v<=100,.(s,of=o,f=v)]); pp<-merge(po,fi2,by="s",allow.cartesian=TRUE)[of<=o+15&o-of<=240][order(s,o,o-of)][,.SD[1],by=.(s,o)]
pfw<-pp[,.(pf=min(pao2/(f/100))),by=s]; e[,`:=`(pf=pfw$pf[match(s,pfw$s)],gcs=NA_real_)]
inf<-Q("select s,o,d from w.einf"); e[,vaso:=as.integer(s%in%inf[o>=-360&o<=1440]$s)]
rc<-Q("select s,vs,ve from w.erc"); e[,imv:=as.integer(s%in%rc[vs>=-360&vs<=1440]$s)]
mc<-as.data.table(readRDS("l_emicro.rds")); mc<-mc[s%in%e$s&!is.na(organism)&organism!=""&!grepl("no growth|coagulase negative|staphylococcus epidermidis|corynebacterium|bacillus|contamin|pending|other",tolower(organism))]
e[,bsi:=as.integer(s%in%mc[o>LM*60&o<=40320]$s)]
rv<-rc[s%in%e$s&vs<40320,.(vh=sum(pmax(0,pmin(ve,40320)-pmax(vs,0)))/60),by=s]; e[,vh:=rv$vh[match(s,rv$s)]]
e[,`:=`(vfd=NA_real_,los=udo/1440,hd=as.integer(hds=="Expired"),db="eICU-CRD",id=s)]
cols<-c("id","db","hosp","age","male","chf","copd","ckd","dm","liver","cancer","sepsis","ards","asp","pna3","pf","gcs","alc0","alc3","alcmin3","alc72","nalc","anc0","ster","ster0","lact","bun","creat","hb","wbc","plt","alb","bili","hr","map","rr","sev","vaso","imv","bsi","vfd","los","hd","death_h","disch_h")
S<-rbind(m[,..cols],e[,..cols]); S[,inc:=as.integer(!is.na(alc3))]; D<-S[inc==1]; D[,x:=as.integer(alc3<1)]
S[,d28:=as.integer(!is.na(death_h)&death_h<=672)]; saveRDS(S,"S.rds")
D[,`:=`(d28=as.integer(!is.na(death_h)&death_h<=672),d90=as.integer(!is.na(death_h)&death_h<=2160),d365=as.integer(!is.na(death_h)&death_h<=8760))]
D[,t28:=fifelse(d28==1,death_h,672)-LM]  # eICU-CRD: discharged alive = 28-day survivor (consistent with binary outcome)
D[db=="eICU-CRD",`:=`(d90=NA,d365=NA)]
saveRDS(D,"D.rds"); saveRDS(fl,"fl.rds"); print(unlist(fl))
print(D[,.(n=.N,d28=round(mean(d28),3),hd=round(mean(hd),3),ster=round(mean(ster),2),bsi=round(mean(bsi),3),imv=round(mean(imv),2)),by=.(db,x)])
print(D[,lapply(.SD,function(z) round(mean(is.na(z)),2)),by=db,.SDcols=c("lact","alb","bili","sev","hr","map","creat","anc0","plt")])
