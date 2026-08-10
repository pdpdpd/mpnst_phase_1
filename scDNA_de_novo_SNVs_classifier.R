options(bitmapType='cairo') #to solve plotting issue on CAMP

library(tidyverse)
library(xgboost)
library(DiagrammeR)
# library(DiagrammeRsvg)

samples = c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) = c("FIT208A3", "FIT208A4", "FIT208A5", "FIT208A6", "FIT208A7", "FIT208A8")

bulk_samples = c("VER236A1", "VER236A2", "VER236A3", "VER236A4", "VER236A5", "VER236A6")

bulk_bam.dir <- paste0("mpnst/bulk/data/merged_bam/")
scDNA_bam.dir <- paste0("mpnst/10X_DNA/data/modified_BAMs/", samples, "/")
pooled_bam.dir <- paste0("mpnst/10X_DNA/data/normal_BAMs/")

snv_analysis.dir = paste0("mpnst/10X_DNA/results/snv_mnv_indel_tumouronly/SNV_analysis/")

output.dir = paste0("mpnst/10X_DNA/results/snv_mnv_indel_tumouronly/SNV_filter")
system(paste0("mkdir -p ", output.dir))
setwd(output.dir)

####################################################################################################################################
###29/1/21 Look at SNVs from bulk normal
####################################################################################################################################
#Prep SNVs locations to genotype
if (F) {
  #Load de novo SNV calls overlap information
  MPNST_regions_overlap_groups <- readRDS(paste0(snv_analysis.dir, "MPNST_bulk_scDNA_BN_scDNA_PN_SNV_regions_overlap_groups.rds"))
  
  #Extract dataframe of SNVs present in scDNA:BN and scDNA:PN but no in bulk
  for (s in samples) {
    raw_shared_SNVs <- data.frame(ID = attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN:",s,"_scDNA_PN")]]],
                                  chr = gsub("_.*", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN:",s,"_scDNA_PN")]]]),
                                  pos = gsub(".*_", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN:",s,"_scDNA_PN")]]]))
    write.table(raw_shared_SNVs[,2:3], file = paste0(s,"_scDNA_BN_scDNA_PN_SNV_loci.txt"), col.names = F, row.names = F, quote = F)
  }
  
  all_shared_snvs <- do.call(rbind, lapply(samples, function(s) {
    read.table(paste0(s,"_scDNA_BN_scDNA_PN_SNV_loci.txt"))
  }))
  all_shared_snvs <- all_shared_snvs[!duplicated(all_shared_snvs),]
  write.table(all_shared_snvs, file = paste0("All_scDNA_BN_scDNA_PN_SNV_loci.txt"), col.names = F, row.names = F, quote = F)
}

#Get allele counts in normal/tumour bulk/sc BAMs
if (F) {
  system(paste0("mkdir -p ",output.dir,"/../scDNA_BN_scDNA_PN_Allele_Freq/"))
  
  NBCORES = 8
  MIN_BASE_QUAL = 20
  # MIN_MAP_QUAL = 0 #set to 0 as cell ranger modifies MAPQ (NOTE THIS IS ONLY TRUE FOR scRNA BAMs)
  MIN_MAP_QUAL = 25
  ALLELECOUNTER = "alleleCounter"
  
  getAlleleCounts <- function (bam.file, output.file, g1000.loci, min.base.qual = 20,
                               min.map.qual = 35, allelecounter.exe = "alleleCounter") {
    cmd = paste(allelecounter.exe, "-b", bam.file, "-l", g1000.loci,
                "-o", output.file, "-m", min.base.qual, "-q", min.map.qual, "-f 0 -F 0 ") #not 10X mode
    system(cmd, wait = T)
  }
  
  #Genotype pooled normal
  print(paste0("Running allele counter on pooled normal"))
  getAlleleCounts(bam.file=paste0(pooled_bam.dir,"sc_merge_mod_normal.bam"),
                  output.file=paste0(output.dir,"/../scDNA_BN_scDNA_PN_Allele_Freq/pooled_alleleFrequencies_all.txt"),
                  g1000.loci=paste0("All_scDNA_BN_scDNA_PN_SNV_loci.txt"),
                  min.base.qual=MIN_BASE_QUAL,
                  min.map.qual=MIN_MAP_QUAL,
                  allelecounter.exe=ALLELECOUNTER)
  
  #Genotype bulk normal
  print(paste0("Running allele counter on bulk normal"))
  getAlleleCounts(bam.file=paste0(bulk_bam.dir, "VER236A7_merged_rmdup_recal.bam"),
                  output.file=paste0(output.dir,"/../scDNA_BN_scDNA_PN_Allele_Freq/normal_alleleFrequencies_all.txt"),
                  g1000.loci=paste0("All_scDNA_BN_scDNA_PN_SNV_loci.txt"),
                  min.base.qual=MIN_BASE_QUAL,
                  min.map.qual=MIN_MAP_QUAL,
                  allelecounter.exe=ALLELECOUNTER)
  
  for (s in samples[1:6]) {
    # s = samples[2]
    #Genotype scDNA tumour
    print(paste0("Running allele counter on ",s))
    getAlleleCounts(bam.file=paste0(scDNA_bam.dir[which(samples == s)],"possorted_mod_bam.bam"),
                    output.file=paste0(output.dir,"/../scDNA_BN_scDNA_PN_Allele_Freq/",s,"_scDNA_alleleFrequencies_all.txt"),
                    g1000.loci=paste0(s,"_scDNA_BN_scDNA_PN_SNV_loci.txt"),
                    min.base.qual=MIN_BASE_QUAL,
                    min.map.qual=MIN_MAP_QUAL,
                    allelecounter.exe=ALLELECOUNTER)
    
    #Genotype bulk tumour
    print(paste0("Running allele counter on ",s))
    getAlleleCounts(bam.file=paste0(bulk_bam.dir,bulk_samples[which(samples == s)],"_merged_rmdup_recal.bam"),
                    output.file=paste0(output.dir,"/../scDNA_BN_scDNA_PN_Allele_Freq/",s,"_tumour_alleleFrequencies_all.txt"),
                    g1000.loci=paste0(s,"_scDNA_BN_scDNA_PN_SNV_loci.txt"),
                    min.base.qual=MIN_BASE_QUAL,
                    min.map.qual=MIN_MAP_QUAL,
                    allelecounter.exe=ALLELECOUNTER)
  }
}

#Prep extra BN SNVs locations to genotype (NOT REURN)
if (F) {
  #Load de novo SNV calls overlap information
  MPNST_regions_overlap_groups <- readRDS(paste0(snv_analysis.dir, "MPNST_bulk_scDNA_BN_scDNA_PN_SNV_regions_overlap_groups.rds"))
  
  #Extract dataframe of SNVs present in scDNA:BN and scDNA:PN but no in bulk
  for (s in samples) {
    raw_shared_SNVs <- data.frame(ID = attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN")]]],
                                  chr = gsub("_.*", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN")]]]),
                                  pos = gsub(".*_", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN")]]]))
    write.table(raw_shared_SNVs[,2:3], file = paste0(s,"_scDNA_BN_SNV_loci.txt"), col.names = F, row.names = F, quote = F)
  }
  
  all_extra_BN_shared_snvs <- do.call(rbind, lapply(samples, function(s) {
    read.table(paste0(s,"_scDNA_BN_SNV_loci.txt"))
  }))
  all_extra_BN_shared_snvs <- all_extra_BN_shared_snvs[!duplicated(all_extra_BN_shared_snvs),]
  write.table(all_extra_BN_shared_snvs, file = paste0("All_scDNA_BN_SNV_loci.txt"), col.names = F, row.names = F, quote = F)
}

#Get extra BN SNVs allele counts in normal/tumour bulk/sc BAMs (NOT REURN)
if (F) {
  system(paste0("mkdir -p ",output.dir,"/../scDNA_BN_Allele_Freq/"))
  
  NBCORES = 8
  MIN_BASE_QUAL = 20
  # MIN_MAP_QUAL = 0 #set to 0 as cell ranger modifies MAPQ (NOTE THIS IS ONLY TRUE FOR scRNA BAMs)
  MIN_MAP_QUAL = 25
  ALLELECOUNTER = "alleleCounter"

  getAlleleCounts <- function (bam.file, output.file, g1000.loci, min.base.qual = 20,
                               min.map.qual = 35, allelecounter.exe = "alleleCounter") {
    cmd = paste(allelecounter.exe, "-b", bam.file, "-l", g1000.loci,
                "-o", output.file, "-m", min.base.qual, "-q", min.map.qual, "-f 0 -F 0 ") #not 10X mode
    system(cmd, wait = T)
  }
  
  #Genotype pooled normal
  print(paste0("Running allele counter on pooled normal"))
  getAlleleCounts(bam.file=paste0(pooled_bam.dir,"sc_merge_mod_normal.bam"),
                  output.file=paste0(output.dir,"/../scDNA_BN_Allele_Freq/pooled_alleleFrequencies_all.txt"),
                  g1000.loci=paste0("All_scDNA_BN_SNV_loci.txt"),
                  min.base.qual=MIN_BASE_QUAL,
                  min.map.qual=MIN_MAP_QUAL,
                  allelecounter.exe=ALLELECOUNTER)
  
  #Genotype bulk normal
  print(paste0("Running allele counter on bulk normal"))
  getAlleleCounts(bam.file=paste0(bulk_bam.dir, "VER236A7_merged_rmdup_recal.bam"),
                  output.file=paste0(output.dir,"/../scDNA_BN_Allele_Freq/normal_alleleFrequencies_all.txt"),
                  g1000.loci=paste0("All_scDNA_BN_SNV_loci.txt"),
                  min.base.qual=MIN_BASE_QUAL,
                  min.map.qual=MIN_MAP_QUAL,
                  allelecounter.exe=ALLELECOUNTER)
  
  for (s in samples[1:6]) {
    # s = samples[2]
    #Genotype scDNA tumour
    print(paste0("Running allele counter on ",s))
    getAlleleCounts(bam.file=paste0(scDNA_bam.dir[which(samples == s)],"possorted_mod_bam.bam"),
                    output.file=paste0(output.dir,"/../scDNA_BN_Allele_Freq/",s,"_scDNA_alleleFrequencies_all.txt"),
                    g1000.loci=paste0(s,"_scDNA_BN_SNV_loci.txt"),
                    min.base.qual=MIN_BASE_QUAL,
                    min.map.qual=MIN_MAP_QUAL,
                    allelecounter.exe=ALLELECOUNTER)
    
    #Genotype bulk tumour
    print(paste0("Running allele counter on ",s))
    getAlleleCounts(bam.file=paste0(bulk_bam.dir,bulk_samples[which(samples == s)],"_merged_rmdup_recal.bam"),
                    output.file=paste0(output.dir,"/../scDNA_BN_Allele_Freq/",s,"_tumour_alleleFrequencies_all.txt"),
                    g1000.loci=paste0(s,"_scDNA_BN_SNV_loci.txt"),
                    min.base.qual=MIN_BASE_QUAL,
                    min.map.qual=MIN_MAP_QUAL,
                    allelecounter.exe=ALLELECOUNTER)
  }
}

#Select 200SNVs, 100 for training, 100 for validation (NOTE when run didn't have R1 SNVs so just selected R2:R5 and P)
if (F) {
  MPNST_regions_overlap_groups <- readRDS(paste0(snv_analysis.dir, "MPNST_bulk_scDNA_BN_scDNA_PN_SNV_regions_overlap_groups.rds"))
  MPNST_all_SNVs <- do.call(rbind, lapply(samples, function(s) {
    return(data.frame(ID = paste0(s,"_",attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN:",s,"_scDNA_PN")]]]),
                      chr = gsub("_.*", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN:",s,"_scDNA_PN")]]]),
                      pos = gsub(".*_", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN:",s,"_scDNA_PN")]]])))
  }))
  #remove SNVs on alt/un contigs
  MPNST_all_SNVs <- MPNST_all_SNVs %>% filter(!str_detect(ID, paste(c("random", "alt", "chrUn", "chrY", "HLA"), collapse = "|")))
  # set.seed(123)
  # sample_snvs <- sample(1:nrow(MPNST_all_SNVs), 200, replace = F)
  # write.csv(MPNST_all_SNVs[sample_snvs,], file = paste0("MPNST_all_scDNA_BN_scDNA_PN_SNV_loci.csv"), row.names = F, quote = F)
  
  #Then manually looked at these 200 SNVs to classify as true or not
  MPNST_sub_SNVs_manual <- read.csv(file = paste0("MPNST_all_scDNA_BN_scDNA_PN_SNV_loci_manual.csv"), header = T)
  
  #Read in bases of 4 BAMs in SNVs and convert to ref and alt
  scDNA_BN_vcf_SNVs <- readRDS("../SNV_analysis/scDNA_BN_vcf_SNVs.rds")
  scDNA_BN_vcf_SNVs <- do.call(rbind, lapply(scDNA_BN_vcf_SNVs, function(s) {
    s %>% mutate(ref = gsub(".*_","",vcf_name) %>% gsub(pattern = "/.", replacement = "")) %>%
      mutate(alt = gsub(".*/","", vcf_name)) %>% mutate(region = gsub(".*_", "", region)) %>% mutate(name = paste0(region, "_", name))
  }))
  
  seq_samples <- c(paste0(samples, "_", c("tumour")), paste0(samples, "_", c("scDNA")), "normal", "pooled")
  #Get allele frequences for all sequencing samples
  scDNA_SNV_allele_freq <- do.call(rbind, lapply(seq_samples, function(seq) {
    ac_table <- read.table(paste0("../scDNA_BN_scDNA_PN_Allele_Freq/", seq, "_alleleFrequencies_all.txt"), header = T)
    colnames(ac_table) <- c("CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth")
    ac_table <- ac_table[,1:6] %>% mutate(seq_sample = seq)
    ac_table <- ac_table %>% rename_at(vars(starts_with("Count_")), function (r) {str_replace(r, "Count_", "")}) %>%
      pivot_longer(cols = 3:6, names_to = "bases", values_to = "Count") %>% group_by(POS) 
    return(ac_table %>% as.data.frame())
  }))
  
  #Add in ref/alt bases and region
  MPNST_sub_SNVs_manual_bases <- MPNST_sub_SNVs_manual[,1:3] %>% mutate(region = gsub("_.*","",ID)) %>%
    left_join(scDNA_BN_vcf_SNVs[,c(4,7,8)], by = c("ID"="name"))
  
  #Add in ref/alt counts
  MPNST_sub_SNVs_manual_ac <- MPNST_sub_SNVs_manual_bases %>% 
    mutate(seq_sample = paste0(region, "_tumour")) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("bulk_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("bulk_alt" = Count) %>%
    mutate(seq_sample = paste0(region, "_scDNA")) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("scDNA_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("scDNA_alt" = Count) %>% 
    mutate(seq_sample = "normal") %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("normal_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("normal_alt" = Count) %>%
    mutate(seq_sample = "pooled") %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("pooled_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("pooled_alt" = Count) %>%
    dplyr::select(-seq_sample)
  
  #Get info re: if SNVs are near centromere
  centromeres <- read_tsv("ref_files/hg38/centromeres.txt.gz", col_names = F) #many contigs each have their own
  colnames(centromeres) <- c("#bin","chrom","chromStart","chromEnd","name")
  centromeres <- centromeres %>% filter(chrom!="chrY")
  centromeres <- centromeres %>% group_by(chrom) %>% summarise(CentroStart = min(chromStart), CentroEnd = max(chromEnd))
  
  MPNST_sub_SNVs_manual_ac <- MPNST_sub_SNVs_manual_ac %>% left_join(centromeres, by = c("chr"="chrom")) %>% 
    mutate(centromere = ifelse (pos > CentroStart-1000000 & pos < CentroEnd+1000000, 1, 0)) %>% dplyr::select(-c(CentroStart,CentroEnd))
  
  #Add in if SNVs are real
  MPNST_sub_SNVs_manual_ac <- cbind(MPNST_sub_SNVs_manual_ac, real=MPNST_sub_SNVs_manual$real)
  
  #Split into train and test
  set.seed(123)
  sample_train_snvs <- sample(1:nrow(MPNST_sub_SNVs_manual_ac), 100, replace = F)
  # MPNST_train_SNVs <- MPNST_sub_SNVs_manual_ac[sample_train_snvs,-16]
  # MPNST_test_SNVs <- MPNST_sub_SNVs_manual_ac[-sample_train_snvs,-16]
  
  #Create matrix for XGBoost
  
  MPNST_train_input <- list()
  MPNST_train_input$data <- as.matrix(MPNST_sub_SNVs_manual_ac[sample_train_snvs,7:15])
  MPNST_train_input$label <- MPNST_sub_SNVs_manual_ac[sample_train_snvs,16]
  
  MPNST_test_input <- list()
  MPNST_test_input$data <- as.matrix(MPNST_sub_SNVs_manual_ac[sample_train_snvs,7:15])
  MPNST_test_input$label <- MPNST_sub_SNVs_manual_ac[sample_train_snvs,16]
  
  #Initial testing
  if (F) {
    #Create xgb.DMatrix and train data set
    MPNST_dtrain <- xgb.DMatrix(data = MPNST_train_input$data, label = MPNST_train_input$label)
    MPNST_bstDMatrix <- xgboost(data = MPNST_dtrain, max.depth = 5, eta = 1, nthread = 2, nrounds = 3, objective = "binary:logistic", verbose = 2, eval.metric = "logloss")
    
    #Predict test data set
    MPNST_pred <- predict(MPNST_bstDMatrix, MPNST_test_input$data) #Prob
    MPNST_prediction <- as.numeric(MPNST_pred > 0.5) #Convert to binary
    
    #Calculate error rate
    MPNST_err <- mean(as.numeric(MPNST_pred > 0.5) != MPNST_test_input$label)
    print(paste("test-error=", MPNST_err))
    
    #Look at output of trees
    MPNST_importance_matrix <- xgb.importance(model = MPNST_bstDMatrix)
    print(MPNST_importance_matrix)
    png(filename = paste0("MPNST_importance_matrix.png"), width = 2000, height = 2000, res = 200)
    xgb.plot.importance(importance_matrix = MPNST_importance_matrix)
    dev.off()
    
    xgb.dump(MPNST_bstDMatrix, with_stats = TRUE)
    xgb.plot.tree(model = MPNST_bstDMatrix)
    # png(filename = paste0("MPNST_scDNA_SNV_trees.png"), width = 2000, height = 2000, res = 200)
    # export_graph(xgb.plot.tree(model = MPNST_bstDMatrix), "MPNST_scDNA_SNV_trees.png", width=2000, height=2000)
    # dev.off()
  }
  
  #Maxime's code for optimsing iterations
  MPNST_xgbcv <- xgb.cv(data = MPNST_train_input$data, label = MPNST_train_input$label,
                        nrounds = 100, nfold = 5, metrics = "logloss", objective = "binary:logistic")
  min.loss.idx <- which.min(MPNST_xgbcv$evaluation_log$test_logloss_mean)
  
  #Train using train data set
  MPNST_dtrain <- xgb.DMatrix(data = MPNST_train_input$data, label = MPNST_train_input$label)
  xgb.DMatrix.save(MPNST_dtrain, "MPNST_dtrain")
  MPNST_bstDMatrix <- xgboost(data = MPNST_dtrain, max.depth = 6, eta = 1, nthread = 2, nrounds = min.loss.idx, objective = "binary:logistic", verbose = 2, eval.metric = "logloss")
  saveRDS(MPNST_bstDMatrix, "MPNST_bstDMatrix_model.rds")
  MPNST_bstDMatrix <- readRDS("MPNST_bstDMatrix_model.rds")
  # xgb.save(MPNST_bstDMatrix, "MPNST_xgboost.model") #save currentl causes bug
  # MPNST_bstDMatrix <- xgb.load("MPNST_xgboost.model")
  
  #Predict then check confusion matrix to see accuracy
  MPNST_predsGB <- predict(MPNST_bstDMatrix, MPNST_train_input$data)
  MPNST_confusionGB <- table(round(MPNST_predsGB), MPNST_train_input$label)
  MPNST_confusionGB
  
  #old
  # 0  1
  # 0 33  3
  # 1  1 63
  
  #new
  # 0  1
  # 0 55  1
  # 1  1 43
  
  #Predict on test data set
  MPNST_pred <- predict(MPNST_bstDMatrix, MPNST_test_input$data) #Prob
  MPNST_prediction <- as.numeric(MPNST_pred > 0.5) #Convert to binary
  
  #Calculate error rate
  MPNST_err <- mean(as.numeric(MPNST_pred > 0.5) != MPNST_test_input$label)
  print(paste("test-error=", MPNST_err))
  # "test-error= 0.04"
  # "test-error= 0.02" #new
  
  #Look at output of trees
  MPNST_importance_matrix <- xgb.importance(model = MPNST_bstDMatrix)
  print(MPNST_importance_matrix)
  png(filename = paste0("MPNST_importance_matrix.png"), width = 2000, height = 2000, res = 200)
  xgb.plot.importance(importance_matrix = MPNST_importance_matrix)
  dev.off()
  
  xgb.dump(MPNST_bstDMatrix, with_stats = TRUE)
  xgb.plot.tree(model = MPNST_bstDMatrix)
  
  #Test on full_manual set
  MPNST_full_manual_input <- list()
  MPNST_full_manual_input$data <- as.matrix(MPNST_sub_SNVs_manual_ac[7:15])
  MPNST_full_manual_input$label <- MPNST_sub_SNVs_manual_ac[16]
  MPNST_pred_full_manual <- predict(MPNST_bstDMatrix, MPNST_full_manual_input$data) #Prob
  MPNST_prediction_full_manual <- as.numeric(MPNST_pred_full_manual > 0.5) #Convert to binary
  MPNST_confusion_full_manual <- table(round(MPNST_pred_full_manual), unlist(MPNST_full_manual_input$label))
  print(MPNST_confusion_full_manual)
  
  #old
  # 0   1
  # 0  53  17
  # 1  19 111
  #new 
  # 0  1
  # 0 99 24
  # 1 16 61
  
  #Calculate error rate
  MPNST_err_full_manual <- mean(as.numeric(MPNST_pred_full_manual > 0.5) != MPNST_full_manual_input$label)
  print(paste("test-error=", MPNST_err_full_manual))
  # "test-error= 0.18"
  # "test-error= 0.2"
}

#Apply to full set of data
if (F) {
  MPNST_regions_overlap_groups <- readRDS(paste0(snv_analysis.dir, "MPNST_bulk_scDNA_BN_scDNA_PN_SNV_regions_overlap_groups.rds"))
  MPNST_all_SNVs <- do.call(rbind, lapply(samples, function(s) {
    return(data.frame(ID = paste0(s,"_",attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN:",s,"_scDNA_PN")]]]),
                      chr = gsub("_.*", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN:",s,"_scDNA_PN")]]]),
                      pos = gsub(".*_", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN:",s,"_scDNA_PN")]]])))
  }))
  #remove SNVs on alt/un contigs
  MPNST_all_SNVs <- MPNST_all_SNVs %>% filter(!str_detect(ID, paste(c("random", "alt", "chrUn", "chrY", "HLA"), collapse = "|")))
  
  #remove duplicates
  MPNST_all_SNVs <- distinct(MPNST_all_SNVs, chr, pos, .keep_all = T)
  
  #Read in bases of 4 BAMs in SNVs and convert to ref and alt
  scDNA_BN_vcf_SNVs <- readRDS("../SNV_analysis/scDNA_BN_vcf_SNVs.rds")
  scDNA_BN_vcf_SNVs <- do.call(rbind, lapply(scDNA_BN_vcf_SNVs, function(s) {
    s %>% mutate(ref = gsub(".*_","",vcf_name) %>% gsub(pattern = "/.", replacement = "")) %>%
      mutate(alt = gsub(".*/","", vcf_name)) %>% mutate(region = gsub(".*_", "", region)) %>% mutate(name = paste0(region, "_", name))
  }))
  
  seq_samples <- c(paste0(samples, "_", c("tumour")), paste0(samples, "_", c("scDNA")), "normal", "pooled")
  #Get allele frequences for all sequencing samples
  scDNA_SNV_allele_freq <- do.call(rbind, lapply(seq_samples, function(seq) {
    ac_table <- read.table(paste0("../scDNA_BN_scDNA_PN_Allele_Freq/", seq, "_alleleFrequencies_all.txt"), header = T)
    colnames(ac_table) <- c("CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth")
    ac_table <- ac_table[,1:6] %>% mutate(seq_sample = seq)
    ac_table <- ac_table %>% rename_at(vars(starts_with("Count_")), function (r) {str_replace(r, "Count_", "")}) %>%
      pivot_longer(cols = 3:6, names_to = "bases", values_to = "Count") %>% group_by(POS) 
    return(ac_table %>% as.data.frame())
  }))
  
  #Add in ref/alt bases and region
  MPNST_all_SNVs_bases <- MPNST_all_SNVs[,1:3] %>% mutate(region = gsub("_.*","",ID)) %>%
    left_join(scDNA_BN_vcf_SNVs[,c(4,7,8)], by = c("ID"="name"))
  MPNST_all_SNVs_bases$pos <- as.integer(levels(MPNST_all_SNVs_bases$pos))[MPNST_all_SNVs_bases$pos]
  
  #Add in ref/alt counts
  MPNST_all_SNVs_ac <- MPNST_all_SNVs_bases %>% 
    mutate(seq_sample = paste0(region, "_tumour")) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("bulk_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("bulk_alt" = Count) %>%
    mutate(seq_sample = paste0(region, "_scDNA")) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("scDNA_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("scDNA_alt" = Count) %>% 
    mutate(seq_sample = "normal") %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("normal_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("normal_alt" = Count) %>%
    mutate(seq_sample = "pooled") %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("pooled_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("pooled_alt" = Count) %>%
    dplyr::select(-seq_sample)
  
  #Get info re: if SNVs are near centromere
  centromeres <- read_tsv("ref_files/hg38/centromeres.txt.gz", col_names = F) #many contigs each have their own
  colnames(centromeres) <- c("#bin","chrom","chromStart","chromEnd","name")
  centromeres <- centromeres %>% filter(chrom!="chrY")
  centromeres <- centromeres %>% group_by(chrom) %>% summarise(CentroStart = min(chromStart), CentroEnd = max(chromEnd))
  
  MPNST_all_SNVs_ac <- MPNST_all_SNVs_ac %>% left_join(centromeres, by = c("chr"="chrom")) %>% 
    mutate(centromere = ifelse (pos > CentroStart-1000000 & pos < CentroEnd+1000000, 1, 0)) %>% dplyr::select(-c(CentroStart,CentroEnd))
  
  #Test on full_manual set
  MPNST_bstDMatrix <- readRDS("MPNST_bstDMatrix_model.rds")
  
  MPNST_all_input <- list()
  MPNST_all_input$data <- as.matrix(MPNST_all_SNVs_ac[7:15])
  MPNST_pred_all <- predict(MPNST_bstDMatrix, MPNST_all_input$data) #Prob
  MPNST_prediction_all <- as.numeric(MPNST_pred_all > 0.5) #Convert to binary
  
  MPNST_pred_all_real <- MPNST_all_SNVs_ac$ID[which(MPNST_prediction_all==1)]
  MPNST_pred_all_unreal <- MPNST_all_SNVs_ac$ID[which(MPNST_prediction_all==0)]
  
  saveRDS(MPNST_all_SNVs_ac[which(MPNST_prediction_all==1),], "MPNST_xbg_pred_real_SNVs.rds")
  #Check how many of original predicted 200 which were real are predicted as real
  print(length(MPNST_sub_SNVs_manual$ID[MPNST_sub_SNVs_manual$real==1]))
  MPNST_pred_all_real %in% MPNST_sub_SNVs_manual$ID[MPNST_sub_SNVs_manual$real==1] %>% table() #TP
  MPNST_pred_all_real %in% MPNST_sub_SNVs_manual$ID[MPNST_sub_SNVs_manual$real==0] %>% table() #FP
  MPNST_pred_all_unreal %in% MPNST_sub_SNVs_manual$ID[MPNST_sub_SNVs_manual$real==0] %>% table() #TN
  MPNST_pred_all_unreal %in% MPNST_sub_SNVs_manual$ID[MPNST_sub_SNVs_manual$real==1] %>% table() #FN
}

#Apply to SNVs only found in scDNA:BN (NOT RERUN)
if (F) {
  MPNST_regions_overlap_groups <- readRDS(paste0(snv_analysis.dir, "MPNST_bulk_scDNA_BN_scDNA_PN_SNV_regions_overlap_groups.rds"))
  MPNST_extra_BN_SNVs <- do.call(rbind, lapply(samples, function(s) {
    return(data.frame(ID = paste0(s,"_",attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN")]]]),
                      chr = gsub("_.*", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN")]]]),
                      pos = gsub(".*_", "", attr(MPNST_regions_overlap_groups[[which(samples == s)]], "elements")[MPNST_regions_overlap_groups[[which(samples == s)]][[paste0(s,"_scDNA_BN")]]])))
  }))
  #remove SNVs on alt/un contigs
  MPNST_extra_BN_SNVs <- MPNST_extra_BN_SNVs %>% filter(!str_detect(ID, paste(c("random", "alt", "chrUn", "chrY", "HLA"), collapse = "|")))
  
  #remove duplicates
  MPNST_extra_BN_SNVs <- distinct(MPNST_extra_BN_SNVs, chr, pos, .keep_all = T)
  
  #Read in bases of 4 BAMs in SNVs and convert to ref and alt
  scDNA_BN_vcf_SNVs <- readRDS("../SNV_analysis/scDNA_BN_vcf_SNVs.rds")
  scDNA_BN_vcf_SNVs <- do.call(rbind, lapply(scDNA_BN_vcf_SNVs, function(s) {
    s %>% mutate(ref = gsub(".*_","",vcf_name) %>% gsub(pattern = "/.", replacement = "")) %>%
      mutate(alt = gsub(".*/","", vcf_name)) %>% mutate(region = gsub(".*_", "", region)) %>% mutate(name = paste0(region, "_", name))
  }))
  
  seq_samples <- c(paste0(samples, "_", c("tumour")), paste0(samples, "_", c("scDNA")), "normal", "pooled")
  #Get allele frequences for all sequencing samples
  scDNA_SNV_allele_freq <- do.call(rbind, lapply(seq_samples, function(seq) {
    ac_table <- read.table(paste0("../scDNA_BN_Allele_Freq/", seq, "_alleleFrequencies_all.txt"), header = T)
    colnames(ac_table) <- c("CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth")
    ac_table <- ac_table[,1:6] %>% mutate(seq_sample = seq)
    ac_table <- ac_table %>% rename_at(vars(starts_with("Count_")), function (r) {str_replace(r, "Count_", "")}) %>%
      pivot_longer(cols = 3:6, names_to = "bases", values_to = "Count") %>% group_by(POS) 
    return(ac_table %>% as.data.frame())
  }))
  
  #Add in ref/alt bases and region
  MPNST_extra_BN_SNVs_bases <- MPNST_extra_BN_SNVs[,1:3] %>% mutate(region = gsub("_.*","",ID)) %>%
    left_join(scDNA_BN_vcf_SNVs[,c(4,7,8)], by = c("ID"="name"))
  MPNST_extra_BN_SNVs_bases$pos <- as.integer(levels(MPNST_extra_BN_SNVs_bases$pos))[MPNST_extra_BN_SNVs_bases$pos]
  
  #Add in ref/alt counts
  MPNST_extra_BN_SNVs_ac <- MPNST_extra_BN_SNVs_bases %>% 
    mutate(seq_sample = paste0(region, "_tumour")) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("bulk_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("bulk_alt" = Count) %>%
    mutate(seq_sample = paste0(region, "_scDNA")) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("scDNA_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("scDNA_alt" = Count) %>% 
    mutate(seq_sample = "normal") %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("normal_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("normal_alt" = Count) %>%
    mutate(seq_sample = "pooled") %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "ref"="bases")) %>% rename("pooled_ref" = Count) %>%
    left_join(scDNA_SNV_allele_freq, by = c("chr"="CHR", "pos"="POS", "seq_sample"="seq_sample", "alt"="bases")) %>% rename("pooled_alt" = Count) %>%
    dplyr::select(-seq_sample)
  
  #Get info re: if SNVs are near centromere
  centromeres <- read_tsv("ref_files/hg38/centromeres.txt.gz", col_names = F) #many contigs each have their own
  colnames(centromeres) <- c("#bin","chrom","chromStart","chromEnd","name")
  centromeres <- centromeres %>% filter(chrom!="chrY")
  centromeres <- centromeres %>% group_by(chrom) %>% summarise(CentroStart = min(chromStart), CentroEnd = max(chromEnd))
  
  MPNST_extra_BN_SNVs_ac <- MPNST_extra_BN_SNVs_ac %>% left_join(centromeres, by = c("chr"="chrom")) %>% 
    mutate(centromere = ifelse (pos > CentroStart-1000000 & pos < CentroEnd+1000000, 1, 0)) %>% dplyr::select(-c(CentroStart,CentroEnd))
  
  #Test on full_manual set
  MPNST_bstDMatrix <- readRDS("MPNST_bstDMatrix_model.rds")
  
  MPNST_extra_BN_input <- list()
  MPNST_extra_BN_input$data <- as.matrix(MPNST_extra_BN_SNVs_ac[7:15])
  MPNST_pred_extra_BN <- predict(MPNST_bstDMatrix, MPNST_extra_BN_input$data) #Prob
  MPNST_prediction_extra_BN <- as.numeric(MPNST_pred_extra_BN > 0.5) #Convert to binary
  
  MPNST_pred_extra_BN_real <- MPNST_extra_BN_SNVs_ac$ID[which(MPNST_prediction_extra_BN==1)]
  MPNST_pred_extra_BN_unreal <- MPNST_extra_BN_SNVs_ac$ID[which(MPNST_prediction_extra_BN==0)]
  
  saveRDS(MPNST_extra_BN_SNVs_ac[which(MPNST_prediction_extra_BN==1),], "MPNST_xbg_pred_real_extra_BN_SNVs.rds")
}

# print(paste0("Reading in allele frequencies on ",s))
# allele_freq_all <- list()
# for (i in chrom_names) {
#   cat(i)
#   allele_freq_all[[i]] = read_tsv(paste0(output.dir,"Allele_Freq/",s,"/",s,"_alleleFrequencies_chr", i, ".txt"))
# }
# allele_freq_all <- do.call(rbind, allele_freq_all)
# #Only keep counts from accepted cells (barcodes variable)
# allele_freq_all$Barcode <- gsub(paste0("(................)-1"), paste0(s, "_\\1"),allele_freq_all$Barcode)
# allele_freq_all <- allele_freq_all %>% filter(Barcode %in% unlist(barcodes[which(samples == s)]))
# 
# saveRDS(allele_freq_all,paste0(output.dir,"Allele_Freq/",s,"/allele_freq_",s,".rds"))
# rm(allele_freq_all)