# 00b extraction: serial lymphocyte/WBC counts (all times, with specimen identifiers), diagnosis sequence numbers,
# ventilation and Glasgow Coma Scale charting, corticosteroid orders, vasopressor infusions,
# procedures and blood cultures from raw MIMIC-IV v3.1 and eICU-CRD v2.0 csv.gz files
suppressMessages({library(DBI);library(data.table)})
MM<-paste0(Sys.getenv("MIMIC_DIR"),"/"); EE<-paste0(Sys.getenv("EICU_DIR"),"/")
con<-dbConnect(duckdb::duckdb(shared_home=FALSE),":memory:"); invisible(dbExecute(con,"PRAGMA threads=8")); Q<-function(s) as.data.table(dbGetQuery(con,s))
r<-function(p) sprintf("read_csv('%s',all_varchar=true,ignore_errors=true)",p)
st<-"(lower(drug) similar to '.*(hydrocortisone|methylprednisolone|prednisone|prednisolone|dexamethasone|solu-medrol|solu-cortef).*')"
saveRDS(Q(sprintf("select hadm_id::bigint hadm_id,specimen_id::bigint sp,itemid::int itemid,charttime::timestamp charttime,try_cast(valuenum as double) v from %s where itemid in ('51244','51133','51301','51256','51265') and hadm_id is not null",r(paste0(MM,"hosp/labevents.csv.gz")))),"l_mlab.rds")
saveRDS(Q(sprintf("select hadm_id::bigint hadm_id,try_cast(seq_num as int) seq,icd_code,icd_version from %s",r(paste0(MM,"hosp/diagnoses_icd.csv.gz")))),"l_mdx.rds")
saveRDS(Q(sprintf("select stay_id::bigint stay_id,itemid::int itemid,charttime::timestamp charttime,value,try_cast(valuenum as double) v from %s where itemid in ('226732','223849','229314','220739','223900','223901','223835')",r(paste0(MM,"icu/chartevents.csv.gz")))),"l_mce.rds")
saveRDS(Q(sprintf("select hadm_id::bigint hadm_id,starttime::timestamp starttime,drug,route from %s where %s",r(paste0(MM,"hosp/prescriptions.csv.gz")),st)),"l_mrx.rds")
saveRDS(Q(sprintf("select stay_id::bigint stay_id,itemid::int itemid,starttime::timestamp starttime from %s where itemid in ('221906','222315','221749','221662','221289')",r(paste0(MM,"icu/inputevents.csv.gz")))),"l_minp.rds")
saveRDS(Q(sprintf("select subject_id::bigint subject_id,charttime::timestamp charttime,chartdate::timestamp chartdate,spec_type_desc,org_name from %s where spec_type_desc='BLOOD CULTURE'",r(paste0(MM,"hosp/microbiologyevents.csv.gz")))),"l_mmicro.rds")
saveRDS(Q(sprintf("select patientunitstayid::bigint s,try_cast(labresultoffset as int) o,labname n,try_cast(labresult as double) v from %s where labname in ('-lymphs','WBC x 1000','-polys','platelets x 1000')",r(paste0(EE,"lab.csv.gz")))),"l_elab.rds")
saveRDS(Q(sprintf("select patientunitstayid::bigint s,try_cast(drugstartoffset as int) o,drugname,routeadmin from %s where lower(drugname) similar to '.*(hydrocortisone|methylprednisolone|prednisone|prednisolone|dexamethasone|solu-medrol|solu-cortef|decadron).*'",r(paste0(EE,"medication.csv.gz")))),"l_emed.rds")
saveRDS(Q(sprintf("select patientunitstayid::bigint s,try_cast(culturetakenoffset as int) o,culturesite,organism from %s where lower(culturesite) like 'blood%%'",r(paste0(EE,"microLab.csv.gz")))),"l_emicro.rds")
saveRDS(Q(sprintf("select patientunitstayid::bigint s,try_cast(treatmentoffset as int) o,lower(treatmentstring) d from %s where lower(treatmentstring) similar to '.*(ecmo|steroid|corticosteroid|hydrocortisone|methylprednisolone).*'",r(paste0(EE,"treatment.csv.gz")))),"l_etr.rds")
cat("done\n")
