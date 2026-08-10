options(bitmapType='cairo') #to solve plotting issue on CAMP

library(tidyverse)
library(VariantAnnotation)
library(VennDiagram)
library(UpSetR)

samples = c("R1", "R2", "R3", "R4", "R5", "P")
names(samples) = c("FIT208A3", "FIT208A4", "FIT208A5", "FIT208A6", "FIT208A7", "FIT208A8")

bulk_samples = c("VER236A1", "VER236A2", "VER236A3", "VER236A4", "VER236A5", "VER236A6")

bulk_snv_bulk_normal_dir = paste0("mpnst/bulk/results/snv_mnv_indel_tumouronly/", bulk_samples, "/", bulk_samples)
scDNA_snv_bulk_normal_dir = paste0("mpnst/10X_DNA/results/snv_mnv_indel_tumouronly/", samples, "_bulk_normal/")
scDNA_snv_pooled_normal_dir = paste0("mpnst/10X_DNA/results/snv_mnv_indel_tumouronly/", samples, "_scDNA_normal/")


output.dir = paste0("mpnst/10X_DNA/results/snv_mnv_indel_tumouronly/SNV_analysis")
system(paste0("mkdir -p ", output.dir))
setwd(output.dir)

####################################################################################################################################
###24/1/21 Look at SNVs from bulk normal
####################################################################################################################################
if (F) {
  samples_to_review <- c(1:6)
  
  bulk_vcf_files <- paste0(bulk_snv_bulk_normal_dir[samples_to_review], "_PASS_snvs_indels.vcf.gz")
  bulk_vcfs <- lapply(bulk_vcf_files, readVcf)
  # vcf_SNVs <- lapply(1:length(samples), function(s) {
  bulk_vcf_SNVs <- lapply(1:length(samples_to_review), function(s) {
    data.frame(chr = as.integer(seqnames(bulk_vcfs[[s]])),
               position = end(bulk_vcfs[[s]]),
               vcf_name = names(bulk_vcfs[[s]]),
               name = paste0(seqnames(bulk_vcfs[[s]]),"_",end(bulk_vcfs[[s]])),
               region = paste0("bulk_", samples[samples_to_review[s]]),
               ref_width = width(bulk_vcfs[[s]]),
               alt_width = str_length(unlist(alt(bulk_vcfs[[s]]))),
               present = T) %>% dplyr::filter(ref_width == 1 & alt_width == 1) %>% dplyr::select(-c(ref_width, alt_width)) #only keep SNVs not indels
  })
  saveRDS(bulk_vcf_SNVs, "bulk_vcf_SNVs.rds")
  
  scDNA_vcf_files <- paste0(scDNA_snv_bulk_normal_dir[samples_to_review], samples[samples_to_review], "_PASS_snvs_indels.vcf.gz")
  scDNA_vcfs <- lapply(scDNA_vcf_files, readVcf)
  scDNA_vcf_SNVs <- lapply(1:length(samples_to_review), function(s) {
    data.frame(chr = as.integer(seqnames(scDNA_vcfs[[s]])),
               position = end(scDNA_vcfs[[s]]),
               vcf_name = names(scDNA_vcfs[[s]]),
               name = paste0(seqnames(scDNA_vcfs[[s]]),"_",end(scDNA_vcfs[[s]])),
               region = paste0("scDNA_", samples[samples_to_review[s]]),
               ref_width = width(scDNA_vcfs[[s]]),
               alt_width = str_length(unlist(alt(scDNA_vcfs[[s]]))),
               present = T) %>% dplyr::filter(ref_width == 1 & alt_width == 1) %>% dplyr::select(-c(ref_width, alt_width)) #only keep SNVs not indels
  })
  saveRDS(scDNA_vcf_SNVs, "scDNA_BN_vcf_SNVs.rds")
  
  ###NOTE (didn't have R1 so excluded)
  #also look at scDNA pooled normal calls (didn't have R1 so excluded)
  scDNA_PN_vcf_files <- paste0(scDNA_snv_pooled_normal_dir[samples_to_review], samples[samples_to_review], "_PASS_snvs_indels.vcf.gz")
  scDNA_PN_vcfs <- lapply(scDNA_PN_vcf_files, readVcf)
  scDNA_PN_vcf_SNVs <- lapply(1:length(samples_to_review), function(s) {
    data.frame(chr = as.integer(seqnames(scDNA_PN_vcfs[[s]])),
               position = end(scDNA_PN_vcfs[[s]]),
               vcf_name = names(scDNA_PN_vcfs[[s]]),
               name = paste0(seqnames(scDNA_PN_vcfs[[s]]),"_",end(scDNA_PN_vcfs[[s]])),
               region = paste0("scDNA_", samples[samples_to_review[s]]),
               ref_width = width(scDNA_PN_vcfs[[s]]),
               alt_width = str_length(unlist(alt(scDNA_PN_vcfs[[s]]))),
               present = T) %>% dplyr::filter(ref_width == 1 & alt_width == 1) %>% dplyr::select(-c(ref_width, alt_width)) #only keep SNVs not indels
  })
  saveRDS(scDNA_PN_vcf_SNVs, "scDNA_PN_vcf_SNVs.rds")
  
  #Plot for each region bulk vs scDNA
  lapply(1:length(samples_to_review), function(s) {
    png(filename = paste0(samples_to_review[s], "_bulk_vs_scDNABN_SNV_venn.png"), width = 2000, height = 2000, res = 200)
    grid.draw(venn.diagram(x = list(as.character(bulk_vcf_SNVs[[s]]$name), as.character(scDNA_vcf_SNVs[[s]]$name)),
                           category.names = c(paste0("bulk_", samples[samples_to_review[s]]), paste0("scDNA_", samples[samples_to_review[s]])),
                           filename = NULL))
    dev.off()
  })
  
  #Plots for all regions
  # png(filename = paste0("All_bulk_vs_scDNA_SNV_venn.png"), width = 2000, height = 2000, res = 200)
  # grid.draw(venn.diagram(x = c(lapply(bulk_vcf_SNVs, function(n) {n$name}), lapply(scDNA_vcf_SNVs, function(n) {n$name})),
  #                        category.names = c(paste0("Bulk_", samples[samples_to_review]), paste0("scDNA_", samples[samples_to_review])),
  #                        filename = NULL))
  # dev.off()
  
  bulk_all_regions_upset_input <- lapply(bulk_vcf_SNVs, function(n) {as.character(n$name)})
  names(bulk_all_regions_upset_input) <- paste0(samples, "_Bulk")
  png(filename = paste0("All_bulk_SNV_upset.png"), width = 2000, height = 2000, res = 200)
  upset(fromList(bulk_all_regions_upset_input), 
        nintersects = 40, 
        nsets = length(bulk_all_regions_upset_input), 
        sets = rev(sort(names(bulk_all_regions_upset_input))), #sets primary first
        order.by = "freq", 
        decreasing = T, 
        mb.ratio = c(0.7, 0.3),
        number.angles = 0,
        text.scale = 1.5,
        point.size = 2.2, 
        line.size = 0.7,
        keep.order = TRUE)
  dev.off()
  
  scDNA_BN_all_regions_upset_input <- lapply(scDNA_vcf_SNVs, function(n) {as.character(n$name)})
  names(scDNA_BN_all_regions_upset_input) <- paste0(samples, "_scDNA_BN")
  png(filename = paste0("All_scDNA_BN_SNV_upset.png"), width = 2000, height = 2000, res = 200)
  upset(fromList(scDNA_BN_all_regions_upset_input), 
        nintersects = 40, 
        nsets = length(scDNA_BN_all_regions_upset_input), 
        sets = rev(sort(names(scDNA_BN_all_regions_upset_input))), #sets primary first
        order.by = "freq", 
        decreasing = T, 
        mb.ratio = c(0.7, 0.3),
        number.angles = 0,
        text.scale = 1.5,
        point.size = 2.2, 
        line.size = 0.7,
        keep.order = T)
  dev.off()
  
  # png(filename = paste0("All_bulk_scDNA_BN_SNV_upset.png"), width = 2000, height = 2000, res = 200)
  # upset(fromList(c(bulk_all_regions_upset_input, scDNA_BN_all_regions_upset_input)), 
  #       nintersects = 40, 
  #       nsets = length(c(bulk_all_regions_upset_input, scDNA_BN_all_regions_upset_input)), 
  #       sets = rev(sort(names(c(bulk_all_regions_upset_input, scDNA_BN_all_regions_upset_input)))),
  #       order.by = "freq", 
  #       decreasing = T, 
  #       mb.ratio = c(0.7, 0.3),
  #       number.angles = 0,
  #       text.scale = 1.5,
  #       point.size = 2.2, 
  #       line.size = 0.7,
  #       keep.order = T)
  # dev.off()
  
  scDNA_PN_all_regions_upset_input <- lapply(scDNA_PN_vcf_SNVs, function(n) {as.character(n$name)})
  names(scDNA_PN_all_regions_upset_input) <- paste0(samples, "_scDNA_PN")
  png(filename = paste0("All_scDNA_PN_SNV_upset.png"), width = 2000, height = 2000, res = 200)
  upset(fromList(scDNA_PN_all_regions_upset_input), 
        nintersects = 40, 
        nsets = length(scDNA_PN_all_regions_upset_input), 
        sets = rev(sort(names(scDNA_PN_all_regions_upset_input))), #sets primary first
        order.by = "freq", 
        decreasing = T, 
        mb.ratio = c(0.7, 0.3),
        number.angles = 0,
        text.scale = 1.5,
        point.size = 2.2, 
        line.size = 0.7,
        keep.order = T)
  dev.off()
  
  #Plot for each region bulk vs scDNA_BN vs scDNA_PN
  lapply(1:length(samples_to_review), function(s) {
    # lapply(1:length(samples_to_review), function(s) {
    png(filename = paste0(samples[s],"_bulk_scDNA_BN_scDNA_PN_SNV_upset.png"), width = 2000, height = 2000, res = 200)
    print(upset(fromList(c(bulk_all_regions_upset_input[s], scDNA_BN_all_regions_upset_input[s], scDNA_PN_all_regions_upset_input[s])), 
                nintersects = 40, 
                nsets = length(c(bulk_all_regions_upset_input[s], scDNA_BN_all_regions_upset_input[s], scDNA_PN_all_regions_upset_input[s])), 
                sets = rev(sort(names(c(bulk_all_regions_upset_input[s], scDNA_BN_all_regions_upset_input[s], scDNA_PN_all_regions_upset_input[s])))), #sets primary first
                order.by = "freq", 
                decreasing = T, 
                mb.ratio = c(0.7, 0.3),
                number.angles = 0,
                text.scale = 1.5,
                point.size = 2.2, 
                line.size = 0.7,
                keep.order = T))
    dev.off()
  })
  
  png(filename = paste0("All_bulk_scDNA_BN_scDNA_PN_SNV_upset.png"), width = 4000, height = 2000, res = 200)
  upset(fromList(c(bulk_all_regions_upset_input, scDNA_BN_all_regions_upset_input, scDNA_PN_all_regions_upset_input)), 
        nintersects = 80, 
        nsets = length(c(bulk_all_regions_upset_input, scDNA_BN_all_regions_upset_input, scDNA_PN_all_regions_upset_input)), 
        sets = rev(sort(names(c(bulk_all_regions_upset_input, scDNA_BN_all_regions_upset_input, scDNA_PN_all_regions_upset_input)))),
        order.by = "freq", 
        decreasing = T, 
        mb.ratio = c(0.7, 0.3),
        number.angles = 0,
        text.scale = 1.5,
        point.size = 2.2, 
        line.size = 0.7, 
        keep.order = TRUE)
  dev.off()
  
  #Extract overlap groups
  overlapGroups <- function (listInput, sort = TRUE) {
    # listInput could look like this:
    # $one
    # [1] "a" "b" "c" "e" "g" "h" "k" "l" "m"
    # $two
    # [1] "a" "b" "d" "e" "j"
    # $three
    # [1] "a" "e" "f" "g" "h" "i" "j" "l" "m"
    listInputmat    <- fromList(listInput) == 1
    #     one   two three
    # a  TRUE  TRUE  TRUE
    # b  TRUE  TRUE FALSE
    #...
    # condensing matrix to unique combinations elements
    listInputunique <- unique(listInputmat)
    grouplist <- list()
    print(paste0("Going through ",nrow(listInputunique)," combinations"))
    # going through all unique combinations and collect elements for each in a list
    for (i in 1:nrow(listInputunique)) {
      if (i%%10 == 0) {print(i)}
      currentRow <- listInputunique[i,]
      myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
      attr(myelements, "groups") <- currentRow
      grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
      myelements
      # attr(,"groups")
      #   one   two three 
      # FALSE FALSE  TRUE 
      #  f  i 
      # 12 13 
    }
    if (sort) {
      grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
    }
    attr(grouplist, "elements") <- unique(unlist(listInput))
    return(grouplist)
    # save element list to facilitate access using an index in case rownames are not named
  }
  
  MPNST_overlap_groups <- overlapGroups(c(bulk_all_regions_upset_input, scDNA_BN_all_regions_upset_input, scDNA_PN_all_regions_upset_input))
  saveRDS(MPNST_overlap_groups, "MPNST_bulk_scDNA_BN_scDNA_PN_SNV_overlap_groups.rds")
  #To get SNV names:
  attr(MPNST_overlap_groups, "elements")[MPNST_overlap_groups[["R2_scDNA_BN:R2_scDNA_PN"]]][c(201,1001)]
  attr(MPNST_overlap_groups, "elements")[MPNST_overlap_groups[["R2_scDNA_BN"]]][c(583,2188)]
  attr(MPNST_overlap_groups, "elements")[MPNST_overlap_groups[["R2_scDNA_PN"]]][c(1032,3042)]
  attr(MPNST_overlap_groups, "elements")[MPNST_overlap_groups[["R2_Bulk:R2_scDNA_BN"]]][c(31,203)]
  attr(MPNST_overlap_groups, "elements")[MPNST_overlap_groups[["R2_Bulk:R2_scDNA_BN:R2_scDNA_PN"]]][c(173,318)]
  
  MPNST_regions_overlap_groups <- lapply(1:length(samples_to_review), function(s) {
    overlapGroups(c(bulk_all_regions_upset_input[s], scDNA_BN_all_regions_upset_input[s], scDNA_PN_all_regions_upset_input[s]))
  })
  saveRDS(MPNST_regions_overlap_groups, "MPNST_bulk_scDNA_BN_scDNA_PN_SNV_regions_overlap_groups.rds")
}


####################################################################################################################################
###12/1/21 Look at SNVs from bulk normal (1st attempt using venn diagrams - not used)
####################################################################################################################################
if (F) {
  samples_to_review <- c(1:6)
  
  bulk_vcf_files <- paste0("mpnst/bulk/results/snv_mnv_indel_tumouronly/", 
                           bulk_samples[samples_to_review], "/", bulk_samples[samples_to_review], "_PASS_snvs_indels.vcf.gz")
  bulk_vcfs <- lapply(bulk_vcf_files, readVcf)
  # vcf_SNVs <- lapply(1:length(samples), function(s) {
  bulk_vcf_SNVs <- lapply(1:length(samples_to_review), function(s) {
    data.frame(chr = as.integer(seqnames(bulk_vcfs[[s]])),
               position = end(bulk_vcfs[[s]]),
               name = paste0(seqnames(bulk_vcfs[[s]]),"_",end(bulk_vcfs[[s]])),
               region = paste0("bulk_", samples[samples_to_review[s]]),
               width = width(bulk_vcfs[[s]]),
               present = T) %>% dplyr::filter(width == 1) %>% dplyr::select(-width) #only keep SNVs not indels
  })
  
  #Tried to look at unfiltered SNVs to see if many of the scDNA ones were those
  # bulk_raw_vcf_files <- paste0("mpnst/bulk/results/snv_mnv_indel_tumouronly/", 
  #                          bulk_samples[samples_to_review], "/", bulk_samples[samples_to_review], "_unfiltered_snvs_indels.vcf.gz")
  # bulk_raw_vcfs <- lapply(bulk_raw_vcf_files, readVcf)
  # bulk_raw_vcf_SNVs <- lapply(1:length(samples_to_review), function(s) {
  #   data.frame(chr = as.integer(seqnames(bulk_raw_vcfs[[s]])),
  #              position = end(bulk_raw_vcfs[[s]]),
  #              name = paste0(seqnames(bulk_raw_vcfs[[s]]),"_",end(bulk_raw_vcfs[[s]])),
  #              region = paste0("bulk_", samples[samples_to_review[s]]),
  #              width = width(bulk_raw_vcfs[[s]]),
  #              present = T) %>% dplyr::filter(width == 1) %>% dplyr::select(-width) #only keep SNVs not indels
  # })
  
  scDNA_vcf_files <- paste0(scDNA_snv_bulk_normal_dir[samples_to_review], samples[samples_to_review], "_PASS_snvs_indels.vcf.gz")
  scDNA_vcfs <- lapply(scDNA_vcf_files, readVcf)
  scDNA_vcf_SNVs <- lapply(1:length(samples_to_review), function(s) {
    data.frame(chr = as.integer(seqnames(scDNA_vcfs[[s]])),
               position = end(scDNA_vcfs[[s]]),
               name = paste0(seqnames(scDNA_vcfs[[s]]),"_",end(scDNA_vcfs[[s]])),
               region = paste0("scDNA_", samples[samples_to_review[s]]),
               width = width(scDNA_vcfs[[s]]),
               present = T) %>% dplyr::filter(width == 1) %>% dplyr::select(-width) #only keep SNVs not indels
  })
  
  lapply(1:length(samples_to_review), function(s) {
    png(filename = paste0(samples_to_review[s], "_bulk_vs_scDNA_SNV_venn.png"), width = 2000, height = 2000, res = 200)
    grid.draw(venn.diagram(x = list(as.character(bulk_vcf_SNVs[[s]]$name), as.character(scDNA_vcf_SNVs[[s]]$name)),
                           category.names = c(paste0("bulk_", samples[samples_to_review[s]]), paste0("scDNA_", samples[samples_to_review[s]])),
                           filename = NULL))
    dev.off()
  })
  
  # lapply(1:length(samples_to_review), function(s) {
  #   png(filename = paste0(samples_to_review[s], "_raw_vs_bulk_vs_scDNA_SNV_venn.png"), width = 2000, height = 2000, res = 200)
  #   grid.draw(venn.diagram(x = list(as.character(bulk_raw_vcf_SNVs[[s]]$name), as.character(bulk_vcf_SNVs[[s]]$name), as.character(scDNA_vcf_SNVs[[s]]$name)),
  #                          category.names = c(paste0("raw_", samples[samples_to_review[s]]), paste0("bulk_", samples[samples_to_review[s]]), paste0("scDNA_", samples[samples_to_review[s]])),
  #                          filename = NULL))
  #   dev.off()
  # })
  
  #only works for up to 5 samples
  # png(filename = paste0("All_bulk_vs_scDNA_SNV_venn.png"), width = 2000, height = 2000, res = 200)
  # grid.draw(venn.diagram(x = c(lapply(bulk_vcf_SNVs, function(n) {n$name}), lapply(scDNA_vcf_SNVs, function(n) {n$name})),
  #                        category.names = c(paste0("Bulk_", samples[samples_to_review]), paste0("scDNA_", samples[samples_to_review])),
  #                        filename = NULL))
  # dev.off()
}