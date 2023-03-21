# All functions needed protein atlas paper

#### 1- General functions ####
all_not_na <- function(x) any(!is.na(x))
all_na <- function(x) all(is.na(x))

#### 2- Overall associations ####

##### A- Conditional logistic regressions ####
#x is the dataframe to use for the analysis
#y is proteins to analyze vector 

###### Non-adjusted ###### 
univ_clogits <- function(x,y){
  univariate_betas <- data.frame(protein=NA,beta=NA, sd=NA, pvalue=NA, z_stat=NA, tot_obs=NA)
  prot <- y
  prot <- prot[prot %in% names(x)]
  for (i in 1:length(prot)){
    univariate_betas[i,"protein"] <- prot[i]
    #print(i)
    clog.fit <-tryCatch(clogit(as.formula(sprintf("case_status_u19~ %s", 
                                                  paste(sapply(prot[i], as.name), "+ strata(caseset)", sep=""))),x)%>%
                          broom::tidy(),
                        error=function(e){data.frame(estimate=NA, std.error=NA,statistic=NA,p.value=NA, tot_obs=NA)})
    nb_obs <- tryCatch(nobs(clogit(as.formula(sprintf("case_status_u19~ %s", paste(sapply(prot[i], as.name), "+ strata(caseset)", sep=""))),x)),error=function(e){nb_obs=NA})
    #print(clog.fit)
    univariate_betas[i,"beta"] <-clog.fit$estimate
    univariate_betas[i,"sd"] <- clog.fit$std.error
    univariate_betas[i,"pvalue"] <- clog.fit$p.value
    univariate_betas[i,"z_stat"] <- clog.fit$statistic
    univariate_betas[i,"tot_obs"] <- nb_obs
  }
  univariate_betas <- univariate_betas%>%mutate(OR=exp(beta))
  #print(y[i])
  #print(univariate_betas)
  return(univariate_betas)}

###### adjusted ###### 
univ_clogits_adj <- function(x,y,adj_par_l){
  univariate_betas <- data.frame(NULL)
  prot <- y
  prot <- prot[prot %in% names(x)]
  for (i in 1:length(prot)){
    for (adj_par in adj_par_l){
      clog.fit.inter <-tryCatch(clogit(as.formula(sprintf("case_status_u19~ %s",paste(sapply(prot[i], as.name), 
                                                                paste0("+",adj_par,"+ strata(caseset)"),
                                                                sep=""))),x)%>%broom::tidy(),
                                error=function(e){data.frame(term=NA, estimate=NA, std.error=NA,statistic=NA,p.value=NA, tot_obs=NA)})
      nb_obs <- tryCatch(nobs(clogit(as.formula(sprintf("case_status_u19~ %s", paste(sapply(prot[i], as.name), paste0("+",adj_par,"+ strata(caseset)"), sep=""))),x)),error=function(e){nb_obs=NA})
      clog.fit.inter <- clog.fit.inter %>% filter(grepl(prot[i],term)) %>% dplyr::mutate(tot_obs=nb_obs, adjustment=adj_par)
      univariate_betas <- bind_rows(univariate_betas,clog.fit.inter)
    }
  }
  univariate_betas <- univariate_betas%>%mutate(OR=exp(estimate),term=str_remove_all(term, "`") )#,adj_var=ifelse(stringr::str_detect(adj_var, "\\+"),"all",as.character(adj_var)))
  return(univariate_betas)}


##### B- proteins+ plco risk AUCs  ####
plcom2012 <- function(data){
  data <- data %>% mutate(across(c("age", "education", "bmi", "copd", "history_cancer", "family_history_binary","race_ethnicity","quit_years","intensity","years_smoked"),~as.numeric(.)))
  plcorisk <- data%>%mutate(logit_plco=(age-62)*0.0778868+(education-4)*(-0.0812744)+(bmi-27)*(-0.0274194)+
                              copd*0.3553063+history_cancer*0.4589971+family_history_binary*0.587185+
                              0*(race_ethnicity==1|race_ethnicity==5|race_ethnicity==6)+ 
                              (race_ethnicity==2)*0.3944778+ (race_ethnicity==3)*(-0.7434744)+
                              (race_ethnicity==4)*(-0.466585)+
                              (quit_years==0)*0.2597431+(((intensity/10)^-1)- 0.4021541613)*(-1.822606)+
                              (years_smoked-27)*0.0317321+
                              (quit_years-10)*(-0.0308572)+(-4.532506))%>%
    mutate(risk_plco=exp(logit_plco)/(1+exp(logit_plco)))
  return(plcorisk)}

#
## function 
aucs_model_prot <- function(x,y, name_file){
  # get table of OR by protein 
  univariate_betas <- univ_clogits(x,y)
  # using the dataset to get aucs 
  datax <- x %>% mutate(across(c("copd","family_history_binary","history_cancer"),~ifelse(is.na(.),0,.)), 
                        race_ethnicity=ifelse(is.na(race_ethnicity) & cohort_id %in% c("NSHDS", "EPIC", "HUNT", "MCCS"),"1",race_ethnicity))
  # matching factors for smoking categories 
  datax <- datax %>% mutate(smoke_cat_match=case_when(smoke_status=="2" & quit_years>10~"1",
                                                      smoke_status=="2" &quit_years<=10~"2",
                                                      smoke_status=="3"&intensity<15~"3",
                                                      smoke_status=="3"&intensity>=15~"4"))
  # get plco risk 
  datax <- plcom2012(datax) %>% mutate(risk=risk_plco)
  # get overall plco AUCs [remove matching- but adjust on matching factors]
  plco_refit <- glm(case_status_u19~logit_plco+age+as.factor(sex)+year_blood_draw+as.factor(smoke_cat_match)+cohort_id, family="binomial", datax) #refit plco
  datax <- datax %>% mutate(plco_overall=plco_refit$fitted.values) #predict refit plco
  nobs <- nobs(plco_refit) #get sample size for plco refit 
  auc_plco_full <- pROC::roc(datax$case_status_u19, datax$plco_overall) #get AUCs
  univariate_betas <- univariate_betas %>% mutate(plco_auc_fd=auc_plco_full$auc[[1]], plco_fd_nobs=nobs)
  univariate_betas <- univariate_betas %>% mutate(auc_plco_protcc_data=NA, auc_prot=NA, auc_plco_prot=NA,
                                                  auc_CI_plco_cc=NA, auc_CI_prot=NA, auc_CI_plco_prot=NA, 
                                                  nobs_plcoauc_cc=NA,nobs_auc_prot=NA, nobs_auc_plcoprot=NA) #ADDING A COLUMN FOR PROTEINS TO COMPLETE

  #get plco+1prot risk
  prot <- y
  prot <- prot[prot %in% names(datax)]
  for (i in prot){
    print(i)
    transform_data_i <- datax[complete.cases(dplyr::select(datax, c("logit_plco",i))),] # complete case protein 
    
    # PLCO AUCs in complete case for comparison // adjust on matching factors
    plco_refit_protdata <- glm(case_status_u19~logit_plco+age+as.factor(sex)+year_blood_draw+as.factor(smoke_cat_match)+cohort_id, family="binomial", transform_data_i) #refit plco
    transform_data_i <- transform_data_i %>% mutate(plco_protdata=plco_refit_protdata$fitted.values) #predict refit plco
    nobs <- nobs(plco_refit_protdata) #get sample size for plco refit 
    auc_plco_protdata<- pROC::roc(transform_data_i$case_status_u19, transform_data_i$plco_protdata, ci=TRUE) #get AUCs
    univariate_betas[univariate_betas$protein== i,"auc_plco_protcc_data"] <-auc_plco_protdata$auc[[1]]
    #CIs for AUCs
    LCI_auc_plco <- auc_plco_protdata$ci[[1]]
    UCI_auc_plco <- auc_plco_protdata$ci[[3]]
    CI_plco <- paste0(formatC(LCI_auc_plco, digits = 2, format="f"), "-",formatC(UCI_auc_plco, digits = 2, format="f") )
    univariate_betas[univariate_betas$protein== i,"auc_CI_plco_cc"] <-CI_plco
    #nobs 
    univariate_betas[univariate_betas$protein== i,"nobs_plcoauc_cc"] <-nobs
    
    #prot AUCs in complete case // adjust on matching factors
    prot_fit <- tryCatch(eval(parse(text=paste0("glm(case_status_u19~age+as.factor(sex)+year_blood_draw+
                                                as.factor(smoke_cat_match)+cohort_id+`",i,"`,family=binomial, data=transform_data_i)"))),
                         error=function(e){NA})
    transform_data_i <- transform_data_i %>% mutate(prot_predict=prot_fit$fitted.values) #predict refit protein
    nobs <- nobs(prot_fit) #get sample size for plco refit 
    auc_prot<- pROC::roc(transform_data_i$case_status_u19, transform_data_i$prot_predict,ci=TRUE) #get AUCs
    univariate_betas[univariate_betas$protein== i,"auc_prot"] <-auc_prot$auc[[1]]
    #CIs for AUCs
    LCI_auc_prot <- auc_prot$ci[[1]]
    UCI_auc_prot <- auc_prot$ci[[3]]
    CI_prot <-  paste0(formatC(LCI_auc_prot, digits = 2, format="f"), "-",formatC(UCI_auc_prot, digits = 2, format="f") )
    univariate_betas[univariate_betas$protein== i,"auc_CI_prot"] <-CI_prot
    #nobs
    univariate_betas[univariate_betas$protein== i,"nobs_auc_prot"] <-nobs
    
 
    #PLCO+prot in complete case// adjust on matching factors
    plco_prot_fit <- tryCatch(eval(parse(text=paste0("glm(case_status_u19~logit_plco+age+as.factor(sex)+
                                                     year_blood_draw+as.factor(smoke_cat_match)+cohort_id+`",i,"`,family=binomial, data=transform_data_i)"))),
                         error=function(e){NA})
    transform_data_i <- transform_data_i %>% mutate(plco_prot_predict=plco_prot_fit$fitted.values) #predict refit protein
    nobs <- nobs(plco_prot_fit) #get sample size for plco refit 
    auc_prot_plco<- pROC::roc(transform_data_i$case_status_u19, transform_data_i$plco_prot_predict, ci=TRUE) #get AUCs
    univariate_betas[univariate_betas$protein== i,"auc_plco_prot"] <-auc_prot_plco$auc[[1]]
    #CIs for AUCs
    LCI_auc_plcoprot <- auc_prot_plco$ci[[1]]
    UCI_auc_plcoprot <- auc_prot_plco$ci[[3]]
    CI_plcoprot <-  paste0(formatC(LCI_auc_plcoprot, digits = 2, format="f"), "-",formatC(UCI_auc_plcoprot, digits = 2, format="f") )
    univariate_betas[univariate_betas$protein== i,"auc_CI_plco_prot"] <-CI_plcoprot
    #nobs
    univariate_betas[univariate_betas$protein== i,"nobs_auc_plcoprot"] <-nobs
 
    #write csv 
    univariate_betas <- univariate_betas %>% mutate(delta_auc_prot=auc_plco_prot-auc_plco_protcc_data)
    write.csv(univariate_betas, paste0(path_save_clean, "tables/",name_file,".csv"), row.names = FALSE)
    
  }
  return(univariate_betas)
}

# #################################################################################################
# #################################################################################################

### 2- Resampling and selection of robust markers ####

#### A- create subsets ####

##### a.1- for common proteins : using all cohorts ##### 
all_cohort_training <- function(data, protein_list){
  #OR data frame 
  ORs_univ <- setNames(as.data.frame(matrix(nrow=1, ncol = 6)), c("protein","set","beta","sd","pvalue","OR"))
  for (i in 1:length(protein_list)) {
    #we should always have EPIC/NSHDS in full
    train_set <- data%>%filter(cohort_id %in% c("EPIC","NSHDS")) 
    #filter out the cohorts/rows were protein i was not measured=> it could affect the training and testing set (all the Na in training or testing depending on seed=> instability in protein selection)
    data_fi <- data%>%filter(!is.na(!!as.symbol(protein_list[i])))
    #calculate the nb of rows to add to the training set in order to have 70% of the data trained
    nb_rows_add <- 0.7*nrow(data_fi)-nrow(train_set)
    #we will pick casesets nb after filtering for only cases, and keep their matched control in the same set. So we sill select/2 the nb of people random and re-assign their matched control
    nb_caseset_select <- nb_rows_add/2
    rest_cohorts <- data_fi%>%filter(!cohort_id %in% c("EPIC","NSHDS"), case_status_u19=="1")
    prop <- nb_caseset_select/nrow(rest_cohorts)
    add_train <- stratified(rest_cohorts, c('cohort_id'), prop, replace = FALSE)
    add_train <- add_train%>%as.data.frame
    add_train_full <- data_fi%>%filter(caseset %in% add_train$caseset)
    train_set <- bind_rows(train_set,add_train_full)
    test_set <- data_fi%>%filter(! participant_id %in% train_set$participant_id)
    
    # Return OR table
    OR_prot <- models_sets(train_set, test_set, protein_list[i])
    #bind rows for or dataframe
    ORs_univ <- bind_rows(ORs_univ, OR_prot)
  }
  ORs_univ <- ORs_univ[-1,]
  return(ORs_univ)
}

##### a.2- for discovery only proteins : EPIC/NSHDS ##### 
splitEPICNSHDS_balanced <- function(data,y){
  #splitting the data
  cases_data <- data%>%filter(cohort_id %in% c("EPIC","NSHDS"),case_status_u19=="1")
  train_data <- stratified(cases_data, c('cohort_id'),0.7 , replace = FALSE)
  train_data <- train_data%>%as.data.frame()
  train_set <-  data%>%filter(caseset %in% train_data$caseset)
  test_set <- data%>%filter(cohort_id %in% c("EPIC","NSHDS"))%>%filter(!caseset %in% train_data$caseset)
  protein <- y
  protein <- protein[protein %in% names(data)]
  #return function for results
  ORs_univ <- models_sets(train_set, test_set, protein)
  return(ORs_univ)
}

#### B- Function for robust marker selection ####

models_sets <- function(train_set,test_set,protein){ #output ORs by training and testing sets 
  ## rescale proteins in both sets in each iteration
  train_set <- train_set %>%group_by(cohort_id) %>%  mutate(across(protein,~scale(.))) %>% ungroup()
  test_set <- test_set %>%group_by(cohort_id) %>%  mutate(across(protein,~scale(.))) %>% ungroup()
  
  #creating a dataset to save OR results
  univariate_betas <- setNames(as.data.frame(matrix(nrow=length(protein)*2, ncol = 5)), c("protein","set","beta","sd","pvalue"))
  univariate_betas[,"protein"] <- rep(protein,2)
  univariate_betas[,"set"] <- rep(c("training","testing"),each=length(protein))
  for (i in 1:length(protein)){
    #training sets
    clog.fit <-tryCatch(clogit(as.formula(sprintf("case_status_u19~ %s", paste(sapply(protein[i], as.name), "+ strata(caseset)", sep=""))),train_set)%>%broom::tidy(),error=function(e){data.frame(estimate=NA,std.error=NA,p.value=NA)})
    univariate_betas[univariate_betas$protein==protein[i] & univariate_betas$set=="training","beta"] <-clog.fit$estimate[1]
    univariate_betas[univariate_betas$protein==protein[i] & univariate_betas$set=="training","sd"] <- clog.fit$std.error[1]
    univariate_betas[univariate_betas$protein==protein[i] & univariate_betas$set=="training" ,"pvalue"] <- clog.fit$p.value[1]
    
    #testing set
    clog.fit.test <-tryCatch(clogit(as.formula(sprintf("case_status_u19~ %s", paste(sapply(protein[i], as.name), "+ strata(caseset)", sep=""))),test_set)%>%broom::tidy(),error=function(e){data.frame(estimate=NA,std.error=NA,p.value=NA)} )
    univariate_betas[univariate_betas$protein==protein[i] & univariate_betas$set=="testing","beta"] <-clog.fit.test$estimate[1]
    univariate_betas[univariate_betas$protein==protein[i] & univariate_betas$set=="testing","sd"] <- clog.fit.test$std.error[1]
    univariate_betas[univariate_betas$protein==protein[i] & univariate_betas$set=="testing" ,"pvalue"] <- clog.fit.test$p.value[1]
  }
  univariate_betas <- univariate_betas%>%mutate(OR=exp(beta))
  #print(univariate_betas)
  return(univariate_betas)
}



#### C- Run loops ####
##### c.1 - for all common proteins : in all data ##### 
sample_testing_replication <- function(data, protlist, n_resample,ent_thresh){
  #Initialisation for datasets to returm 
  ## data for only selected (yes/no selected in iteration i) and (beta in iteration i selcteed, beta replication in iteration i select)=> MERGE in betas later when iterating 
  protein_replicat_time <- setNames(as.data.frame(matrix(nrow=length(protlist), ncol = n_resample+2)), c("protein","times",seq(1:n_resample)))
  protein_replicat_time$protein <- protlist
  ## data set with proteins in column and betas in all 500 iterations in training and testing  
  protein_beta_all <- data.frame(protein=rep(protlist, each=2), set=rep(c("training","testing"),times=length(protlist)))
  print(head(protein_beta_all))
  ## dataset with proteins with nominal significance in training and testing (iteration yes/no )
  protein_nomin_sig <- setNames(as.data.frame(matrix(nrow=length(protlist), ncol = n_resample+2)), c("protein","times",seq(1:n_resample)))
  protein_nomin_sig$protein <- protlist
  
  # Bootsraps iterations
  for(i in 1:n_resample){
    set.seed(50+5*i)
    #all betas for all proteins in training and testing 
    OR_allcoh_full <- all_cohort_training(data,protlist)
    #data set with betas for replicated proteins in training and testing 
    OR_allcoh <- OR_allcoh_full%>%mutate(sig=case_when(pvalue<0.05/ent_thresh & set=="training"~1,pvalue<0.05 & set=="testing"~1, TRUE~0))%>%group_by(protein)%>%add_count(sig)%>%filter(sig==1 &n==2)
    
    ## Dataset for replicated protein
    # return dataset of replicated proteins with ENT sig : add betas testing and betas training  
    rep_prot <- OR_allcoh$protein
    protein_replicat_time <- protein_replicat_time%>%mutate("{i}" := ifelse(protein %in% rep_prot,1,0))
    #merge in train (add in beta estimates) and rename
    OR_allcoh_train <- OR_allcoh%>%filter(set=="training"& protein %in% rep_prot)%>%select(protein,beta)
    name <- paste0("beta_train",i)
    protein_replicat_time <- merge(protein_replicat_time,OR_allcoh_train, by="protein",all.x=TRUE, all.y=TRUE)%>%rename(!!name := beta)
    #merge in betas testing and rename 
    OR_allcoh_test<- OR_allcoh%>%filter(set=="testing"& protein %in% rep_prot)%>%select(protein,beta)
    name <- paste0("beta_test",i)
    protein_replicat_time <- merge(protein_replicat_time,OR_allcoh_test, by="protein",all.x=TRUE, all.y=TRUE)%>%rename(!!name := beta)
    #new column for nb of time selected
    protein_replicat_time <- protein_replicat_time%>%mutate(times=rowSums(.[as.character(c(1:i))]))
    # save it within function, if it stops nothing will be lost 
    write.csv(protein_replicat_time, paste0(path_save_clean, "bootstraps_tables/Replicated_prot_boot.csv"), row.names = FALSE)
    
    ## All betas for all proteins in training and testing
    betas_all_sets <- OR_allcoh_full%>%dplyr::select(protein, set, beta)
    protein_beta_all <- OR_allcoh_full%>%dplyr::select(protein, set, beta)%>%merge(protein_beta_all, by=c("protein","set"),all.x=TRUE, all.y=TRUE)%>%rename(!!paste0("beta_",i):=beta)
    # save 
    write.csv(protein_beta_all, paste0(path_save_clean, "bootstraps_tables/betas_allprot_bothsets.csv"), row.names = FALSE)
    
    ## proteins with nominal sig 
    OR_nomsig <- OR_allcoh_full%>%mutate(sig=case_when(pvalue<0.05 & set=="training"~1,pvalue<0.05 & set=="testing"~1, TRUE~0))%>%group_by(protein)%>%add_count(sig)%>%filter(sig==1 &n==2)
    prot_nomsig <- OR_nomsig$protein
    protein_nomin_sig <- protein_nomin_sig%>%mutate("{i}" := ifelse(protein %in% prot_nomsig,1,0))%>%mutate(times=rowSums(.[as.character(c(1:i))]))
    #save
    write.csv(protein_nomin_sig, paste0(path_save_clean, "bootstraps_tables/prot_nomsig_inboth.csv"), row.names = FALSE)
  }
  return(protein_replicat_time)
}

##### c.2- for discovery only proteins : EPIC/NSHDS ##### 
sample_testing_replication_en <- function(data, protlist, n_resample,ent_thresh){
  #Initialisation for datasets to returm 
  ## data for only selected (yes/no selected in iteration i) and (beta in iteration i selcteed, beta replication in iteration i select)=> MERGE in betas later when iterating 
  protein_replicat_time <- setNames(as.data.frame(matrix(nrow=length(protlist), ncol = n_resample+2)), c("protein","times",seq(1:n_resample)))
  protein_replicat_time$protein <- protlist
  ## data set with proteins in column and betas in all 500 iterations in training and testing  
  protein_beta_all <- data.frame(protein=rep(protlist, each=2), set=rep(c("training","testing"),times=length(protlist)))
  ## dataset with proteins with nominal significance in training and testing (iteration yes/no )
  protein_nomin_sig <- setNames(as.data.frame(matrix(nrow=length(protlist), ncol = n_resample+2)), c("protein","times",seq(1:n_resample)))
  protein_nomin_sig$protein <- protlist
  
  # Bootstrap for rep
  for(i in 1:n_resample){
    set.seed(50+5*i)
    # splitting training and testing sets and get betas for all prot 
    OR_EPICNSHDS_full <- splitEPICNSHDS_balanced(data,protlist) #splits abd returns univariable stats by set
    #data set with betas for replicated proteins in training and testing 
    OR_EPICNSHDS <- OR_EPICNSHDS_full%>%mutate(sig=case_when(pvalue<0.05/ent_thresh & set=="training"~1,pvalue<0.05 & set=="testing"~1, TRUE~0))%>%group_by(protein)%>%add_count(sig)%>%filter(sig==1 &n==2)
    print(paste0("OK split iter ",i))
    ## Dataset for replicated protein
    # return dataset of replicated proteins with ENT sig : add betas testing and betas training
    rep_prot <- OR_EPICNSHDS$protein
    protein_replicat_time <- protein_replicat_time%>%mutate("{i}" := ifelse(protein %in% rep_prot,1,0))
    #merge in train to get info on estimates and rename 
    OR_EPICNSHDS_train <- OR_EPICNSHDS%>%filter(set=="training"& protein %in% rep_prot)%>%dplyr::select(protein,beta)
    name <- paste0("beta_",i)
    protein_replicat_time <- merge(protein_replicat_time,OR_EPICNSHDS_train, by="protein",all.x=TRUE, all.y=TRUE)%>%rename(!!name := beta)
    #merge in testing and rename 
    OR_EPICNSHDS_test<- OR_EPICNSHDS%>%filter(set=="testing"& protein %in% rep_prot)%>%dplyr::select(protein,beta)
    name <- paste0("beta_test",i)
    protein_replicat_time <- merge(protein_replicat_time,OR_EPICNSHDS_test, by="protein",all.x=TRUE, all.y=TRUE)%>%rename(!!name := beta)
    #new column for nb of time selected
    protein_replicat_time <- protein_replicat_time%>%mutate(times=rowSums(.[as.character(c(1:i))]))
    # save it within function, if it stops nothing will be lost 
    write.csv(protein_replicat_time, paste0(path_save_clean, "bootstraps_tables/ReplicatedEPICNSHDS_prot_boot.csv"), row.names = FALSE)
    print(paste0("OK glm and nmerge iter ",i))
    ## All betas for all proteins in training and testing (before we were saving betas for rep protein, in replicated iterations)
    betas_all_sets <- OR_EPICNSHDS_full%>%dplyr::select(protein, set, beta)
    protein_beta_all <- OR_EPICNSHDS_full%>%dplyr::select(protein, set, beta)%>%merge(protein_beta_all, by=c("protein","set"),all.x=TRUE, all.y=TRUE)%>%rename(!!paste0("beta_",i):=beta)
    # save 
    write.csv(protein_beta_all, paste0(path_save_clean, "bootstraps_tables/betasEPICNSHDS_allprot_bothsets.csv"), row.names = FALSE)
    ## proteins with nominal sig 
    OR_nomsig <- OR_EPICNSHDS_full%>%mutate(sig=case_when(pvalue<0.05 & set=="training"~1,pvalue<0.05 & set=="testing"~1, TRUE~0))%>%group_by(protein)%>%add_count(sig)%>%filter(sig==1 &n==2)
    prot_nomsig <- OR_nomsig$protein
    protein_nomin_sig <- protein_nomin_sig%>%mutate("{i}" := ifelse(protein %in% prot_nomsig,1,0))%>%mutate(times=rowSums(.[as.character(c(1:i))]))
    #save
    write.csv(protein_nomin_sig, paste0(path_save_clean, "bootstraps_tables/protEPICNSHDS_nomsig_inboth.csv"), row.names = FALSE)
    
  }
  return(protein_replicat_time)
}

# #################################################################################################
# #################################################################################################

strat_clogits <- function(x,y,coltofilter){
  OR_strat_full <- setNames(as.data.frame(matrix(nrow=1, ncol = 9)), c("protein","beta","sd","pvalue","z_stat","Variable","strat_column","tot_obs","cases/controls"))
  for (k in 1:length(coltofilter)) { #for every column to stratify 
    strat <- full_data_scaled%>%pull(coltofilter[k])%>%as.factor()%>%levels #get strata levels as factor 
    for (i in 1:length(strat)) { #within every stratum of the column
      x_st <- x%>%filter(!!as.symbol(coltofilter[k]) == strat[i]) #filter to the needed data of strat
      #clogits of every protein within strata (x_st) and create a dataframe with protein name beta and pvalue 
      ORstrat <- univ_clogits(x_st,y)
      ORstrat <- ORstrat%>%mutate(Variable=coltofilter[k],strat_column=strat[i],`cases/controls`=NA) #add new columns for each strata 
      for(j in 1:length(y)){ #y is the list of proteins - measure the number of cases and controls on which it was measured
        N_cases <- x_st%>%filter(!is.na(!!as.symbol(y[j])))%>%filter(case_status_u19=="1")%>%nrow
        N_cases <- ifelse(is.na(N_cases),0,as.numeric(N_cases))
        N_ctrl <- x_st%>%filter(!is.na(!!as.symbol(y[j])))%>%filter(case_status_u19=="0")%>%nrow
        N_ctrl <- ifelse(is.na(N_ctrl),0,as.numeric(N_ctrl))
        ORstrat[ORstrat$protein==y[j],"cases/controls"] <- paste0(N_cases,"/",N_ctrl)
      }
      OR_strat_full <- bind_rows(OR_strat_full,ORstrat)
    }
  }
  #get the overall clogit not stratified for every protein
  OR_full <- univ_clogits(x,y)
  OR_full <- OR_full%>%mutate(strat_column="Overall")
  OR_strat_full <- bind_rows(OR_strat_full,OR_full)
  OR_strat_full <- OR_strat_full[-1,]
  return(OR_strat_full)
}

#x is the dataframe to use for the analysis
#y is proteins to analyze vector 
#strata is vector of strata 

plot_strat_OR <- function(data){
  plt <- data%>%ggplot(aes(x=protein, y=OR))+geom_pointrange(aes(ymin=exp(log(OR)-1.96*sd), ymax=exp(log(OR)+1.96*sd), color=pvalue))+geom_hline(yintercept = 1,linetype="dashed")+coord_flip()+facet_wrap(~strat_column)+theme(legend.position = "none")+theme_bw()
  print(plt)
  return(plt)
}


chisq_strata_prot <- function(data,protein_l , strata){
  data2 <- data%>%mutate(phet=NA)%>%mutate(phet=as.numeric(phet))
  for (i in 1:length(protein_l)) {
    data_prot <- data%>%filter(protein==protein_l[i])
    for(j in 1:length(strata)){
      data_strat <- data_prot%>%filter(Variable==strata[j])%>%filter(strat_column!="Other/NOS")%>%filter(!is.na(beta)|!is.na(sd))
      yi <- data_strat$beta # observed effect size estimates
      vi <- data_strat$sd 
      vi <- vi^2 # corresponding sampling variances
      het <- tryCatch(metafor::rma(yi, vi, method="FE"), error=function(e){return(NA)}) # fixed model effect
      if(!is.na(het)){
        p_het <- het$QEp
      }
      else{
        p_het <- NA
      }
      data2[data2$protein==protein_l[i]&data2$Variable==strata[j],"phet"] <- p_het
    }
    
  }
  return(data2)
}
