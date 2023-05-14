library(biomaRt)
library(VariantAnnotation)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(RCurl)


load_data <- function(control_vcf_file_path, test_vcf_file_path){
  #' load test and control SNPs and remove like SNPs from the dataset
  ctrl_vcf <- readVcf(open(VcfFile(file = control_vcf_file_path,
                                   index = "mouse-snps-all.annots.vcf.gz.tbi")))
  
  test_vcf <- readVcf(open(VcfFile(file = test_vcf_file_path,
                                   index = "mouse-snps-all.annots.vcf.gz.tbi")))
  # Filter Control SNPs from the test SNP dataset
  # This also removes duplicate SNPs
  vcf <- test_vcf[!ctrl_vcf@rowRanges@ranges %in% test_vcf@rowRanges@ranges]
  return(vcf)
}


get_annos <- function(csv_file_path){
  #' Finds the appropriate column of gene IDs, 
  #' extracts the gene_IDs and converts them to Ensemble Mouse IDs
  
  raw_anno <- read.csv(csv_file_path, stringsAsFactors=F)
  
  entrezIDs <- mapIds(org.Mm.eg.db, raw_anno$Symbol,"ENTREZID","SYMBOL")
  
  return(entrezIDs)
}

SNP_filter <- function(vcf, entrez_ID_list, var_type){
  vcf <- vcf[vcf@assays@data@listData$FI==1]
  vcf <- vcf[lapply(vcf@info@listData$CSQ, function(x) strsplit(x,split="\\|")[[1]][4]) %in% names(annos)]
  list_1 <- lapply(vcf@info@listData$CSQ, function(x) strsplit(x, split="\\|")[[1]][3]) == "HIGH"
  list_2 <- lapply(vcf@info@listData$CSQ, function(x) strsplit(x, split="\\|")[[1]][25]) >= 4
  list_2[is.na(list_2)] <- F
  list_3 <- unlist(lapply(vcf@info@listData$CSQ, function(x) str_detect(x, "deleterious\\(")))
  
  vcf <- vcf[list_1|list_2|list_3]
  
  #load known all mouse genes
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  
  #match chromosome names to txdb
  vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf)))
  seqlevels(vcf)[20] <- "chrM"
  # Locate Variants
  loc <- locateVariants(rowRanges(vcf), txdb, AllVariants(promoter = PromoterVariants(400000, 20000), intergenic = IntergenicVariants(400000, 20000)))

  #So this is a bit shit, but txdb mm10 has out of bounds in 
  loc <- loc[!duplicated(loc$QUERYID),]
  loc$CSQ <- vcf@info@listData$CSQ
  
  loc$GENENAME <- sapply(vcf@info@listData$CSQ, function(x) strsplit(x, split="\\|")[[1]][4])
  
  
  
  # filter by mouse ensembl gene IDs
  #loc[mcols(loc)$GENEID %in% na.omit(entrez_ID_list)]
  return(loc)
}

# init ensembl
ensembl = useEnsembl(biomart="ensembl",dataset="mmusculus_gene_ensembl")

BM_annotate <- function(vcf){
  
  #print(searchAttributes(ensembl))
  #values <- getBM(attributes=c('description','goslim_goa_description',, 'wikigene_description', "transcript_tsl",
  #                             'name_1006','phenotype_description'),filters = 'ensembl_gene_id', values = names(vcf), mart = ensembl)
  
  # convert gene_name to entreziDs for biomart search
  entrezIDs <- na.omit(unlist(mapIds(org.Mm.eg.db, unlist(vcf$GENENAME),"ENTREZID","SYMBOL")))
  

  # get a bunch of annotations from biomart
  tmp <- getBM(attributes=c('external_gene_name'), values = na.omit(entrezIDs), filters = 'entrezgene_id', mart = ensembl, verbose=T)
  
  phenos <- getBM(attributes=c('external_gene_name','phenotype_description','strain_name'), values = tmp$external_gene_name, filters = 'external_gene_name', mart = ensembl, verbose=T)
  
  GO <- getBM(attributes=c('external_gene_name', 'name_1006','goslim_goa_description'), values = tmp$external_gene_name, filters = 'external_gene_name', mart = ensembl, verbose=T)
  
  reactome <- getBM(attributes=c('external_gene_name','reactome', 'biogrid'), values = tmp$external_gene_name, filters = 'external_gene_name', mart = ensembl, verbose=T)
  
  mrna <- getBM(attributes=c('external_gene_name', 'rnacentral', 'rfam', 'transcript_appris', 'transcript_tsl', 'arrayexpress'), values = tmp$external_gene_name, filters = 'external_gene_name', mart = ensembl, verbose=T)
  
  # merge into a crazy dataframe 
  df_list <- list(phenos,GO,reactome,mrna)
  databases <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list) 
  
  
  # group by gene name and reduce to all unique values for each column
  grouped_SNPs <- databases %>% group_by(external_gene_name)
  annoed_SNPs <- grouped_SNPs %>% summarise_all(function(x) list(unique(x)))
  
  
  annoed_SNPs <- merge(vcf,annoed_SNPs, by.x="GENENAME", by.y="external_gene_name", all=T)                      
  
  return(annoed_SNPs)
   
}
 


SNP_data <- load_data("results_C57Bl6_v7.vcf.gz","results_C3H_HeH_v6.vcf.gz")


annos <- c(get_annos("Nodal_Signalling_GO.csv"), 
           get_annos("Sonic_HedgeHog_GO.csv"), 
           get_annos("FGF_Signalling_GO.csv"),
           get_annos("Pusapati_cilia_or_pathy.csv"))

nodal <- get_annos("Nodal_Signalling_GO.csv")
ssh <-get_annos("Sonic_HedgeHog_GO.csv")
# Quality Control
cil <- get_annos("Pusapati_cilia_or_pathy.csv")

# CSQs <- lapply(vcf@info@listData$CSQ, function(x) {
#   strsplit(x,split="\\|")[[1]][4]
# })x,split="\\|"


# SNPs in annotated genes


Exon_SNPs <- SNP_filter(SNP_data, annos, AllVariants())

BM_SNPs <- BM_annotate(Exon_SNPs)


fwrite(BM_SNPs, "GRCm38_results.csv")

# need to free up memory
rm(SNP_data, Exon_SNPs)

## ENSEMBL MOUSE STRAINS



m_ensembl = useEnsembl(biomart="snp",dataset="mmusculus_snp")


# these are just results from MGI SNP query that I did manually
f1 <- read.table(file = 'MGI_SNP_Query_20230119_071431.txt', sep = '\t', header = TRUE)
f2 <- read.table(file = 'MGI_SNP_Query_20230119_071923.txt', sep = '\t', header = TRUE)
f3 <- read.table(file = 'MGI_SNP_Query_20230119_071724.txt', sep = '\t', header = TRUE)


#combine all the things
file <- merge(f1, f2, by=intersect(names(f1), names(f2)), all=T)
file <- merge(file, f3, by=intersect(names(file), names(f3)), all=T)
rm(f1, f2, f3)

# SNPs 
# filtered to be like in C3H and 129 but unlike C57BL6
file_1 <- file %>% dplyr::filter((C3H.HeJ != C57BL.6J & C3H.HeJ != "" & C57BL.6J != "" & X129S6.SvEvTac != C57BL.6J & X129S6.SvEvTac != "")) %>% 
  filter(Gene.Symbol %in% names(annos))

# filtered to be like in C3H and FVB but unlike C57BL6
file_2 <- file %>% dplyr::filter((C3H.HeJ != C57BL.6J & C3H.HeJ != "" & C57BL.6J != "" & FVB.NJ != C57BL.6J & FVB.NJ != "")) %>% 
  filter(Gene.Symbol %in% names(annos))

# just unlike between C57Bl6/N and C3H/HeH

file_3 <- file %>% dplyr::filter((C3H.HeJ != C57BL.6J & C3H.HeJ != "" & C57BL.6J != "")) %>%  filter(Gene.Symbol %in% names(annos))

rm(file)

fetch_SNPs_from_BM <- function(mgi_query_snps, annos){
  
  pathway_ensemblIDs <- mapIds(org.Mm.eg.db, names(annos),"ENSEMBL","SYMBOL")
  
  
  functional_info_1 <- getBM(attributes=c('refsnp_id', 'sift_prediction', 'reg_consequence_types'),
                           filters='snp_filter', values=unique(mgi_query_snps$SNP.ID..GRCm39.), mart=m_ensembl, verbose=T)
  
  functional_info_2 <- getBM(attributes=c('refsnp_id', 'motif_consequence_types','ensembl_transcript_stable_id','ensembl_gene_name'),
                             filters='snp_filter', values=unique(mgi_query_snps$SNP.ID..GRCm39.), mart=m_ensembl, verbose=T)
  
  functional_info <- merge(functional_info_1, functional_info_2)
  
  trans_stab_info <- getBM(attributes=c('ensembl_transcript_id','transcript_tsl','transcript_is_canonical','transcript_appris'),
                           filters='ensembl_transcript_id', values=unique(functional_info$ensembl_transcript_stable_id), mart=ensembl, verbose=T)
  
  functional_info <- merge(functional_info, trans_stab_info, by.x='ensembl_transcript_stable_id', by.y = 'ensembl_transcript_id', all.x = T)
  
  filtered_functional_info <- functional_info %>% 
    filter(ensembl_gene_name %in% pathway_ensemblIDs) %>%
    filter((str_detect(sift_prediction, 'deleterious') | (reg_consequence_types != "" | motif_consequence_types!= "")) & (str_detect(transcript_tsl, 'tsl1')| transcript_is_canonical == 1)) 
  
  
  genome_pos <- getBM(attributes=c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'ensembl_transcript_chrom_strand','ensembl_gene_name'),
                             filters='snp_filter', values=unique(filtered_functional_info$refsnp_id), mart=m_ensembl, verbose=T)

  
  snp_information <- getBM(attributes=c('refsnp_id','title','ensembl_peptide_allele', 'motif_name','reg_feature_stable_id', 'ensembl_gene_name', 'sift_score','consequence_type_tv'),
                           filters='snp_filter', values=unique(filtered_functional_info$refsnp_id), mart=m_ensembl, verbose=T)
  
  
  phenos <- getBM(attributes=c('external_gene_name','phenotype_description', 'ensembl_gene_id'), values = unique(snp_information$ensembl_gene_name), 
                  filters = 'ensembl_gene_id', mart = ensembl, verbose=T)
  
  GO <- getBM(attributes=c('external_gene_name','name_1006', 'goslim_goa_description','ensembl_gene_id'), values = unique(snp_information$ensembl_gene_name),
              filters = 'ensembl_gene_id', mart = ensembl, verbose=T)
  
  reactome <- getBM(attributes=c('external_gene_name','reactome','biogrid','ensembl_gene_id'), values = unique(snp_information$ensembl_gene_name), 
                    filters = 'ensembl_gene_id', mart = ensembl, verbose=T)
  
  ensembl_list <- lapply(list(phenos,GO,reactome), as.data.table)
  
  
  
  ens_databases <- rbindlist(ensembl_list, use.names = TRUE, fill = TRUE) %>% group_by(ensembl_gene_id) %>% summarise_all(function(x) paste(unique(na.omit(x)), collapse = "|"))
  
  
  m_ensembl_list <- lapply(list(filtered_functional_info, genome_pos, snp_information), as.data.table)
  
  m_ens_databases <- rbindlist(m_ensembl_list, use.names = TRUE, fill = TRUE) %>% group_by(refsnp_id, ensembl_gene_name) %>%
    filter(!is.na(ensembl_gene_name)) %>% summarise_all(function(x) paste(unique(na.omit(x)), collapse = "|"))
  
  annoed_SNPs <- merge(ens_databases, m_ens_databases, by.x='ensembl_gene_id', by.y = 'ensembl_gene_name', all.y=T)
   
  grouped_SNPs <- annoed_SNPs %>% group_by(refsnp_id)
  
  annoed_SNPs <- grouped_SNPs %>% summarise_all(function(x) list(unique(x)))
  
  annoed_SNPs <- merge(annoed_SNPs, as.data.table(mgi_query_snps), by.x = 'refsnp_id', by.y='SNP.ID..GRCm39.', all.x=T)
  
  return(annoed_SNPs)
}

#'goslim_goa_description', 'reactome'




snps_filtered_from_just_C3H_HEH <- fetch_SNPs_from_BM(file_3, annos)





MGI_Results <- unnest(snps_filtered_from_just_C3H_HEH)


fwrite(MGI_Results, "MGI_query_results_GRCm39.csv")


consequences_of_int <- c("NMD_transcript_variant", "3_prime_UTR_variant","non_coding_transcript_exon_variant", "5_prime_UTR_variant", "splice_region_variant")

consequences_from_C3H <- snps_filtered_from_just_C3H_HEH %>% filter(str_detect(consequence_type_tv,
                                                                               paste(consequences_of_int, collapse = "|")))

sifted_just_C3H_HEH <- snps_filtered_from_just_C3H_HEH %>% filter(str_detect(sift_prediction, 'deleterious'))

fwrite(sifted_just_C3H_HEH, "MGI_query_grcm39_results_sifted.csv")


UTR_variants_from_C3H <- snps_filtered_from_just_C3H_HEH %>% filter(str_detect(consequence_type_tv,"3_prime_UTR_variant|5_prime_UTR_variant"))

splicing_variants_from_C3H <- snps_filtered_from_just_C3H_HEH %>% filter(str_detect(consequence_type_tv,"donor|acceptor"))



fwrite(splicing_variants_from_C3H, "GRCm39_splice_variants.csv")



filename <- "gerp_conservation_scores.mus_musculus.GRCm39.bw"

#download bw file if it doesn't exist
if (!file.exists(filename)) {
  url <- "http://ftp.ensembl.org/pub/current_compara/conservation_scores/90_mammals.gerp_conservation_score/gerp_conservation_scores.mus_musculus.GRCm39.bw"
  
  download.file(url, destfile = filename, mode = "wb", 
                method = "curl", extra = "-L", 
                quiet = FALSE, 
                timeout = 60, 
                progress = function(downloaded, total) {
                  message(paste0("Downloaded ", round(downloaded/total * 100, 1), "% of ", file_out))
                  return(TRUE)
                })
}



tester <- snps_filtered_from_just_C3H_HEH %>% mutate(chr_name=unlist(chr_name),
                                                       chrom_start=unlist(chrom_start), 
                                                       chrom_end=unlist(chrom_end))


tester$ensembl_transcript_chrom_strand[1 %in% tester$ensembl_transcript_chrom_strand] <- "+"

tester$ensembl_transcript_chrom_strand[-1 %in% tester$ensembl_transcript_chrom_strand] <- "-"



snp_gr <-makeGRangesFromDataFrame(tester, seqnames.field = "chr_name", start.field = "chrom_start", end.field = "chrom_end", strand.field = "ensembl_transcript_chrom_strand")


bw_data <- import(filename)

#bw_data <- renameSeqlevels(bw_data, paste0("chr", seqlevels(bw_data)))

#seqlevels(bw_data)[27] <- "chrM"


bw_data_subset <- subsetByOverlaps(bw_data, snp_gr, ignore.strand=T)

snp_gr_subset <- subsetByOverlaps(snp_gr, bw_data, ignore.strand=T)

sorted_snp_gr <- sort(snp_gr_subset)

sorted_bw <- sort(bw_data_subset)


mcols(sorted_snp_gr)$score <- mcols(sorted_bw)$score

GERP_df <-data.frame(sorted_snp_gr)

colnames(GERP_df)[6] <- "GERP_score"



merged_df <- merge(GERP_df, snps_filtered_from_just_C3H_HEH, by.x = c("seqnames", "start", "end"), by.y = c("chr_name","chrom_start","chrom_end"))


fwrite(merged_df, "MGI_query_grcm39_results_with_GERP.csv")


merged_df_filt <- merged_df %>% filter(GERP_score < -5)

fwrite(merged_df_filt, "MGI_query_grcm39_results_with_GERP_filt.csv")




# Get data from MouseMine


#, token = , "YOUR-API-KEY"




MouseMine_GXD_query<-function(db){
  service2 <- import("intermine.webservice")
  
  service_obj2 <- service2$Service("https://www.mousemine.org/mousemine/service", token = "71s8Q7g49aj8HbOfFbS3")
  
  query <- service_obj2$new_query("GXDExpression")
  
  query$add_view("assayType", "feature.symbol","feature.primaryIdentifier", "stage", "age",
                 "structure.name", "strength", "pattern", "genotype.symbol", "sex",
                 "assayId", "probe", "image", "publication.mgiJnum")
  
  
  query$add_sort_order("GXDExpression.assayId", "ASC")
  query$add_constraint("feature.organism.taxonId", "=", "10090", code = "B")
  query$add_constraint("feature", "IN", "link_1", code = "A")
  query$add_constraint("feature.symbol", "ONE OF", db$Gene.Symbol, code = "C")
  
  
  results <- query$results()
  
  res_expression <- reticulate::iterate(results)
  return(res_expression)
}



MouseMine_query <- function(db){
  library(reticulate)
  # Use the use_python() function to specify the path to your Python executable
  use_python("/usr/bin/python3")
  
  # Import the intermine.webservice module
  service <- import("intermine.webservice")
  
  # Create a Service object
  service_obj <- service$Service("https://www.mousemine.org/mousemine/service", token = "71s8Q7g49aj8HbOfFbS3")
  
  # Create a new OntologyAnnotation query object
  query <- service_obj$new_query("SequenceFeature")

  query$add_constraint("ontologyAnnotations.ontologyTerm", "GOTerm")
  # Add constraints to the query object

  query$add_view(
    "primaryIdentifier", "symbol", "name", "mgiType",
    "ontologyAnnotations.ontologyTerm.namespace",
    "ontologyAnnotations.qualifier",
    "ontologyAnnotations.ontologyTerm.identifier",
    "ontologyAnnotations.ontologyTerm.name",
    "ontologyAnnotations.evidence.code.code",
    "ontologyAnnotations.evidence.withText",
    "ontologyAnnotations.evidence.publications.mgiJnum",
    "ontologyAnnotations.evidence.publications.pubMedId",
    "ontologyAnnotations.evidence.publications.citation",
    "ontologyAnnotations.evidence.comments.type",
    "ontologyAnnotations.evidence.comments.description"
  )

  query$add_constraint("ontologyAnnotations.dataSets.name", "=", "GOTerm to Mouse Feature Annotations from MGI", code = "B")
  query$add_constraint("organism.taxonId", "=", "10090", code = "D")
  query$add_constraint("SequenceFeature", "IN", "link_1", code = "C")
  query$add_constraint("symbol", "ONE OF", db$Gene.Symbol, code = "E")

  query$outerjoin("ontologyAnnotations.evidence.comments")

  # Execute the query and print the results
  results <- query$results()
  res_ontology <- reticulate::iterate(results)
  return(res_ontology)
}


onto <- MouseMine_query(sifted_just_C3H_HEH)


gxd <- MouseMine_GXD_query(sifted_just_C3H_HEH)




gxd_terms <- lapply(gxd, function(x){
  gene_symb <- x$feature$symbol
  struct <- ifelse(is.null(x$structure$name), "", x$structure$name)
  stage <- ifelse(is.null(x$stage), "", x$stage)
  age <- ifelse(is.null(x$age), "", x$age)
  assay_type <- ifelse(is.null(x$assayType), "", x$assayType)
  assay_id <- ifelse(is.null(x$assayId), "", x$assayId)
  pattern <- ifelse(is.null(x$pattern), "", x$pattern)
  exp_strength <- ifelse(is.null(x$strength), "", x$strength)
  genotype <- ifelse(is.null(x$genotype$symbol), "", x$genotype$symbol)
  probe <- ifelse(is.null(x$probe), "", x$probe)
  emaps <- ifelse(is.null(x$emaps), "", x$emaps)
  publication <- ifelse(is.null(x$publication$mgiJnum), "", x$publication$mgiJnum)
  spec_no. <- ifelse(is.null(x$specimenNum), "", x$specimenNum)
  data.frame(gene_symb, struct, stage, age, spec_no., assay_id, assay_type, probe, pattern, exp_strength, genotype, emaps, publication)
  
})

gxd_df <- do.call(rbind, gxd_terms) %>% filter(str_detect(age, "embryonic")) 

gxd_sum <-gxd_df %>% group_by(gene_symb) %>% summarise_all(function(x) paste(unique(na.omit(x)), collapse = "|"))


full_df_sifted <- merge(sifted_just_C3H_HEH, gxd_sum, by.x="Gene.Symbol", by.y="gene_symb", all.x=T)


gene_terms <- lapply(onto, function(x) {
  gene <- x$subject$symbol
  term <- x$ontologyTerm$name
   id <-  x$ontologyTerm$identifier
  data.frame(gene = gene, description = term)
})

onto_df <- do.call(rbind, gene_terms)








session <- browserSession("UCSC")


genome(session) <- "mm39"







tester$chr_name <- unlist(tester$chr_name)

tester$chrom_start <- unlist(tester$chrom_start)

#tester$chrom_end <- unlist(tester$chrom_end)+1

tester$end <- unlist(tester$end)+1
tester$strand <- unlist(tester$strand)

snp_gr_v2 <-makeGRangesFromDataFrame(tester, seqnames.field = "seqnames", start.field = "start", end.field = "end", strand.field = "strand")



#snp_gr <-makeGRangesFromDataFrame(tester, seqnames.field = "chr_name", start.field = "chrom_start", end.field = "chrom_end", strand.field = "strand")

snp_gr_v2 <- renameSeqlevels(snp_gr_v2, paste0("chr", seqlevels(snp_gr_v2)))


for_query <- GRangesForUCSCGenome("mm39", chrom = snp_gr_v2@seqnames, ranges = snp_gr_v2@ranges, strand=snp_gr_v2@strand)




# so lapply the table query 

session <- browserSession("UCSC")
genome(session) <- "mm39"


# test rtracklayer with conserversation 35-way scores

con_scores <- sapply(seq_along(for_query), function(i) {
  Sys.sleep(10)
  getTable(ucscTableQuery(session, track="cons35way",range=for_query[i], table="multiz35way"))})




con_df <- data.frame(con_scores)
con_df_t <-t(as.matrix(con_df))

con_df_t_df <- data.frame(matrix(unlist(con_df_t), ncol = 7, byrow = F))

colnames(con_df_t_df) <- c("bin", "chrom", "chromStart", "chromEnd", "name", "multi35way_score", "strand")


good_con_scores <- select(con_df_t_df, "chrom", "chromStart", "chromEnd", "multi35way_score")

good_con_scores$chromEnd <- as.numeric(good_con_scores$chromEnd) - 1

df_merged_v2<- merge(merged_df, good_con_scores,
                     by.x = c("seqnames", "start", "end"), 
                     by.y=c("chrom", "chromStart", "chromEnd"))




# download mammals.99 GERP conservation file: 




# MGP projects dumb merged variant file




test_vcf <- readVcf(open(VcfFile(file = "mgp_REL2021_snps_C3H_BL6.vcf.gz",
                                 index = "mgp_REL2021_snps_C3H_BL6.vcf.gz.tbi")))




