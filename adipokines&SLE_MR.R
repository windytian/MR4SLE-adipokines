library(TwoSampleMR)
library(MRPRESSO)

########################SLE & adipokines univariable MR analysis##########################

#resistin
resistin <- extract_instruments(outcomes="ebi-a-GCST90012034",clump = T)
SLE<- extract_outcome_data(snps = resistin$SNP, outcomes ="finn-b-M13_SLE",proxies=T,maf_threshold =0.01,access_token=NULL)

mydata <- harmonise_data(resistin,SLE,action=2)
mydata <- mydata[which(mydata$mr_keep == TRUE),]
res <- mr(mydata)
het <- mr_heterogeneity(mydata) 
pleio <- mr_pleiotropy_test(mydata)
print(het); print(pleio)
presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdExposure = 'se.exposure',
                    SdOutcome = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                    data = mydata, NbDistribution = 1000,  SignifThreshold = 0.05)
presso$`Main MR results`

outliers <- presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
mydata_sub <- mydata[-outliers,]

res2 <- mr(mydata_sub)
het2 <- mr_heterogeneity(mydata_sub)
pleio2 <- mr_pleiotropy_test(mydata_sub)
print(het2); print(pleio2)



mr_results<-mr(mydata_sub)
mr_method_list()
generate_odds_ratios(mr_results)

res_single <- mr_singlesnp(mydata_sub)
plot2 <- mr_forest_plot(res_single)
plot2
plot4 <- mr_funnel_plot(res_single)
plot4

rm(list=ls())

#leptin
leptin <- extract_instruments(outcomes="ebi-a-GCST90007307",clump = T)
SLE<- extract_outcome_data(snps = leptin$SNP, outcomes ="finn-b-M13_SLE",proxies=T,maf_threshold =0.01,access_token=NULL)

mydata <- harmonise_data(leptin,SLE,action=2)
mydata <- mydata[which(mydata$mr_keep == TRUE),]
res <- mr(mydata)
het <- mr_heterogeneity(mydata) 
pleio <- mr_pleiotropy_test(mydata)
print(het); print(pleio)
presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdExposure = 'se.exposure',
                    SdOutcome = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                    data = mydata, NbDistribution = 1000,  SignifThreshold = 0.05)
presso$`Main MR results`

outliers <- presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
mydata_sub <- mydata[-outliers,]

res2 <- mr(mydata_sub)
het2 <- mr_heterogeneity(mydata_sub)
pleio2 <- mr_pleiotropy_test(mydata_sub)
print(het2); print(pleio2)

mr_results<-mr(mydata_sub)
mr_method_list()
generate_odds_ratios(mr_results)

res_single <- mr_singlesnp(mydata_sub)
plot2 <- mr_forest_plot(res_single)
plot2
plot4 <- mr_funnel_plot(res_single)
plot4

rm(list=ls())

#Adiponectin
adi <- extract_instruments(outcomes="ieu-a-1",clump = T)
SLE<- extract_outcome_data(snps = adi$SNP, outcomes ="finn-b-M13_SLE",proxies=T,maf_threshold =0.01,access_token=NULL)

mydata <- harmonise_data(adi,SLE,action=2)
mydata <- mydata[which(mydata$mr_keep == TRUE),]
res <- mr(mydata)
het <- mr_heterogeneity(mydata) 
pleio <- mr_pleiotropy_test(mydata)
print(het); print(pleio)
presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",SdExposure = 'se.exposure',
                    SdOutcome = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                    data = mydata, NbDistribution = 1000,  SignifThreshold = 0.05)
presso$`Main MR results`

outliers <- presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
mydata_sub <- mydata[-outliers,]

res2 <- mr(mydata_sub)
het2 <- mr_heterogeneity(mydata_sub)
pleio2 <- mr_pleiotropy_test(mydata_sub)
print(het2); print(pleio2)

mr_results<-mr(mydata_sub)
mr_method_list()
generate_odds_ratios(mr_results)

res_single <- mr_singlesnp(mydata_sub)
plot2 <- mr_forest_plot(res_single)
plot2
plot4 <- mr_funnel_plot(res_single)
plot4


#MVMR 
id_exposure <- c("ebi-a-GCST90007307","ebi-a-GCST90012034") 
id_outcome <- "finn-b-M13_SLE"
exposure_dat <- mv_extract_exposures(id_exposure)
dim(exposure_dat)

oucome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome) 
mvdat <- mv_harmonise_data(exposure_dat,  oucome_dat) 
res <- mv_multiple(mvdat) 
res






