library(data.table)
library(readr)
library(readxl)
library(dplyr)
library(progress)
library(DescTools)
library(ggplot2)
library(reshape2)
library(grid)

if(!require(DescTools)){install.packages("DescTools")}

# biomart has been fucked somehow lately
#library(biomaRt)
if (!requireNamespace("progress", quietly = TRUE)) {
  install.packages("progress")
}


library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)


mapfun <- function(mousegenes){
  gns <- mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
  naind <- is.na(mapped$Homo_sapiens)
  hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
  out$Human_symbol[!naind] <- hsymb
  return(out$Human_symbol)
}





#Initialise house and mouse marts for BiomaRt- nope fuck BiomaRt lately
#human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",  mirror = "useast")
#mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "asia")






# Read datasets

Mouse_Models <- read_csv("HDisease_MModel.csv")



test <- Mouse_Models %>% filter()
Human_Genes_and_Diseases <- read_csv("HGene_HDisease.csv")
diseases_for_HP_0000006_AD <- read_excel("diseases_for_HP_0000006_AD.xlsx")
diseases_for_HP_0000007_AR <- read_excel("diseases_for_HP_0000007_AR.xlsx")
diseases_for_HP_0001419_XLR <- read_excel("diseases_for_HP_0001419_XLR.xlsx")
diseases_for_HP_0001423_XLD <- read_excel("diseases_for_HP_0001423_XLD.xlsx")



# Prepare data - i.e. only get rows where the synonmym is an OMIM ID
Mouse_Models <- subset(Mouse_Models, grepl("^OMIM", OntologyAnnotation.ontologyTerm.synonyms.name))
Human_Genes_OMIM <- subset(Human_Genes_and_Diseases, grepl("^OMIM", Gene.ontologyAnnotations.ontologyTerm.synonyms.name))




Filter_by_Human_Inheritance <- function(Hum_Genes, Mouse_Data, Hum_Inheritance_Data, same_gene_in_man){
  # Get shared OMIM numbers between Hum_Genes and Inheritance Data and filter
  intersection_human_disease_names <-  dplyr::intersect(Hum_Genes$Gene.ontologyAnnotations.ontologyTerm.synonyms.name, Hum_Inheritance_Data$DISEASE_ID)
  
  Human_Inh_genes <- Hum_Genes[Hum_Genes$Gene.ontologyAnnotations.ontologyTerm.synonyms.name %in% intersection_human_disease_names,]

  
  # Now get shared DOIDs between Human Genes and Mouse Models
  
  Hum_Mouse_intersection <-dplyr::intersect(Human_Inh_genes$Gene.ontologyAnnotations.ontologyTerm.identifier, Mouse_Data$OntologyAnnotation.ontologyTerm.identifier)
  
  # Filter for human genes where there is a mouse  model of the diseases
  
  Human_Inh_genes_with_MM <- Human_Inh_genes[Human_Inh_genes$Gene.ontologyAnnotations.ontologyTerm.identifier %in% Hum_Mouse_intersection,]

  
  # Filter Mouse_Models using the intersected disease names
  filtered_Mouse_Models <- Mouse_Data[Mouse_Data$OntologyAnnotation.ontologyTerm.identifier %in% Hum_Mouse_intersection, ] %>% 
    distinct(OntologyAnnotation.ontologyTerm.identifier, OntologyAnnotation.subject.symbol, OntologyAnnotation.subject.background.name, .keep_all = T)
  
  #stuff it  - let's just do the thing
  
  if (same_gene_in_man){
    #pb <- progress_bar$new(total = length(filtered_Mouse_Models$OntologyAnnotation.subject.alleles.feature.symbol))
    
    filtered_Mouse_Models$HGNC_Symbol <- mapfun(filtered_Mouse_Models$OntologyAnnotation.subject.alleles.feature.symbol)
    #pb$terminate()
    
    # reshape HGNC_symbol
    for (i in seq_along(filtered_Mouse_Models$HGNC_Symbol)) {
      filtered_Mouse_Models$HGNC_Symbol[[i]] <- gsub("[<>]", "", filtered_Mouse_Models$HGNC_Symbol[[i]][1])
    }
    
    print(filtered_Mouse_Models$HGNC_Symbol)
    
    Hum_sub <- Human_Inh_genes_with_MM %>% dplyr::select(Gene.symbol, Gene.ontologyAnnotations.ontologyTerm.identifier)
    
    
    
    final_dataset <- inner_join(filtered_Mouse_Models, Hum_sub, 
                              by = c("HGNC_Symbol" = "Gene.symbol",
                                     "OntologyAnnotation.ontologyTerm.identifier" = "Gene.ontologyAnnotations.ontologyTerm.identifier"))
    
    print(str(final_dataset))
    
    View(final_dataset)
    return(final_dataset)
    
  }
  else{
    return(filtered_Mouse_Models)
  }
  
  
  
  
  # Split dataset into different conditions
  
  
  # so that that we have our data split, convert the MGI symbol to Human HGNC and make sure that that it's in the human gene list (for that inheritance model)
  
  
}



filtered_Mouse_Models_AD <- Filter_by_Human_Inheritance(Human_Genes_OMIM, Mouse_Models, diseases_for_HP_0000006_AD, TRUE)


filtered_Mouse_Models_AR <- Filter_by_Human_Inheritance(Human_Genes_OMIM, Mouse_Models, diseases_for_HP_0000007_AR, TRUE) 

filtered_Mouse_Models_XLR <- Filter_by_Human_Inheritance(Human_Genes_OMIM, Mouse_Models, diseases_for_HP_0001419_XLR, TRUE)


filtered_Mouse_Models_XLD <- Filter_by_Human_Inheritance(Human_Genes_OMIM, Mouse_Models, diseases_for_HP_0001423_XLD, TRUE)
  


Create_Confusion_Matrix <- function(AD_data, AR_data){
  
  AR_subsets <- split(AR_data, AR_data$OntologyAnnotation.subject.zygosity)
  AR_summary <- lapply(AR_subsets, function(subset) nrow(subset))
  
  AD_subsets <- split(AD_data, AD_data$OntologyAnnotation.subject.zygosity)
  AD_summary <- lapply(AD_subsets, function(subset) nrow(subset))
  
  AR_df <- data.frame(Zygosity = names(AR_summary), AR_Count = unlist(AR_summary))
  AD_df <- data.frame(Zygosity = names(AD_summary), AD_Count = unlist(AD_summary))
  
  # Merge the two data frames based on Zygosity
  combined_df <- merge(AR_df, AD_df, by = "Zygosity", all = TRUE)
  
  colnames(combined_df) <- c("Mouse_Zygosity", "Human_AR", "Human_AD")
  str(combined_df)
  
  melted_df <- melt(combined_df, id.vars = "Mouse_Zygosity",
                    measure.vars = c("Human_AR", "Human_AD"),
                    variable.name = "Human_Category",
                    value.name = "Count")
  
  dotplot <- ggplot(melted_df, aes(x = Mouse_Zygosity, y = Human_Category, size = Count, color=Human_Category)) +
    geom_point() +
    labs(title = "Dot Plot of Mouse Zygosity across Human AR and AD",
         x = "Mouse Zygosity", y = "Human Inheritance") +
    scale_size_continuous(range = c(1, 10)) +  # Adjust the range as needed
    theme_minimal()
  # Print the dot plot
  print(dotplot)
  
  subset_df <- combined_df[combined_df$Mouse_Zygosity %in% c("ht", "hm"), ]
  
  print(subset_df)
  
  
  # Create the confusion matrix
  # Create the confusion matrix manually
  
  # Print the confusion matrix
  subset_matrix <- as.matrix(subset_df[, c("Human_AR", "Human_AD")])
  
  print(subset_matrix)
  
  # Create the confusion matrix
  confusion_matrix <- subset_matrix
  rownames(confusion_matrix) <- c("Recessive", "Dominant")
  colnames(confusion_matrix) <- c("Recessive", "Dominant")
  # Print the confusion matrix
  print(confusion_matrix)
  
  confusion_matrix_normalized <-  sweep(confusion_matrix,2,colSums(confusion_matrix),`/`) * 100
  
  print(confusion_matrix_normalized)
  
  melted_confusion_matrix <- melt(confusion_matrix_normalized)
  melted_confusion_matrix$percentage <- melted_confusion_matrix$value
  melted_confusion_matrix$formatted_label <- lapply(seq_len(nrow(melted_confusion_matrix)), function(i) {
    paste(melt(confusion_matrix)$value[i],
          " (", sprintf("%.2f%%", melted_confusion_matrix$percentage[i]), ")", sep = "")
  })
  
  confusion_matrix_plot <- ggplot(melted_confusion_matrix, aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value), color = "white") +
    geom_text(aes(label = formatted_label),
              vjust = 0.5, hjust = 0.5, size = 5)+
    scale_fill_gradient(low="white", high="#009194")+
    labs(title = "",
         x = "Inheritance Pattern in Mouse",
         y = "Inheritance Pattern in Human", 
         fill = "Frequency (%)")+
    theme(axis.title.x = element_text(size = 14),  # Adjust the font size for the x-axis label
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),  # Adjust the font size for x-axis tick labels
          axis.text.y = element_text(size = 12), 
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 12),
          legend.position = "right")  # Adjust the font size for the y-axis label
    
  
  # Print the confusion matrix plot
  print(confusion_matrix_plot)
  
  
  
  #final_bit_is stat time!
  #get rowsums for expected_values
  
  rowSumA <- sum(confusion_matrix[1, ])
  rowSumB <- sum(confusion_matrix[2, ])
  
  expected_matrix <- matrix(c(0, rowSumA, rowSumB, 0), nrow = 2, byrow = TRUE)
  
  #result <- chisq.test(confusion_matrix, p = expected_matrix, correct = FALSE)  # Set correct = FALSE if you don't want Yates' continuity correction
  
  GTest(x=confusion_matrix,
          p=expected__matrix,
          correct="none")  
  
  # Print the result
  
}

Create_Confusion_Matrix(filtered_Mouse_Models_XLD, filtered_Mouse_Models_XLR)


Create_Confusion_Matrix(filtered_Mouse_Models_AD, filtered_Mouse_Models_AR)


### create a confusion matrix specific to HPE ###


Human_HPE_genes <- Human_Genes_OMIM %>% filter(str_detect(Gene.ontologyAnnotations.ontologyTerm.name, "holoprosencephaly"))




Missed_Human_genes <- Human_Genes_OMIM %>% filter(Gene.symbol %in% c("CENPF", "DISP1", "DLL1", "FGFR1", "KMT2D", "PPP1R12A", "SMC1A", "SMC3", "STIL"))



corrected_HPE_genes <- rbind(Human_HPE_genes, Missed_Human_genes)



test <- Mouse_Models %>% filter(str_detect(OntologyAnnotation.ontologyTerm.name, "holoprosencephaly"))

filtered_mouse_HPE_AD <- filtered_Mouse_Models_AD %>% filter(HGNC_Symbol %in% corrected_HPE_genes$Gene.symbol)

filtered_mouse_HPE_AR <- filtered_Mouse_Models_AR %>% filter(HGNC_Symbol %in% corrected_HPE_genes$Gene.symbol)






# Determine background differences in HPE:
library(stringr)


mouse_genes_pheno <- read_csv("_Genotype_Phenotype.csv")


library(tidyr)


HPE_mouse <- mouse_genes_pheno %>%
  filter(stringr::str_detect(OntologyAnnotation.ontologyTerm.name, "holoprosencephaly"))



HPE_causing_allele_mouse <- mouse_genes_pheno %>% 
  #filter(OntologyAnnotation.subject.alleles.symbol %in% unique(HPE_mouse$OntologyAnnotation.subject.alleles.symbol), !(OntologyAnnotation.subject.zygosity %in% c("cx", "cn"))) %>%
  filter(OntologyAnnotation.subject.alleles.symbol %in% unique(HPE_mouse$OntologyAnnotation.subject.alleles.symbol)) %>%
  group_by(OntologyAnnotation.subject.symbol, OntologyAnnotation.subject.background.name) %>%
  summarize(Zygosity = unique(OntologyAnnotation.subject.zygosity), 
            ContainsHoloprosencephaly = any(str_detect(OntologyAnnotation.ontologyTerm.name, "holoprosencephaly")),
            CombinedPhenotypes = paste(OntologyAnnotation.ontologyTerm.name, collapse = ", "),
            CombinedAlleles = paste(unique(OntologyAnnotation.subject.alleles.symbol), collapse = ", "),
            CombinedPubmed = paste(unique(OntologyAnnotation.evidence.publications.pubMedId), collapse = ", ")) %>% 
            filter(!str_detect(OntologyAnnotation.subject.symbol, "Del\\(|In\\(|Dp\\(|Gt\\(")) # this is needed to be a second line
  
HPE_strains <- str_remove_all(HPE_causing_allele_mouse$OntologyAnnotation.subject.background.name,"involves: ") %>%
  str_remove_all("\\(") %>%
  str_remove_all("\\)") %>%
  str_remove_all("either:") %>%
  str_remove_all(" ") %>%
  str_remove("-.*$") %>%
  str_remove("or.*$") %>%
  str_replace( "\\.", "\\*") %>%
  str_split( "\\*") 
  




HPE_causing_allele_mouse$HPE_strains <- HPE_strains


unique_alleles <- unique(HPE_causing_allele_mouse$CombinedAlleles)


test_v2 <- HPE_causing_allele_mouse %>% filter(CombinedAlleles == "Smad2<+>")

for (i in 1:length(unique_alleles)) {
  current_row_name <- unique_alleles[i]
  print(current_row_name)
  rows_to_average <- HPE_causing_allele_mouse %>% filter(CombinedAlleles == current_row_name)
  print(rows_to_average)
  
  
  if (nrow(rows_to_average) > 1) {
    aggregated_row <- apply(rows_to_average, 2, min, na.rm = TRUE)
    aggregated_matrix[i, ] <- aggregated_row}
  else {
    aggregated_matrix[i, ] <- rows_to_average}
  
}

indeed_HPE <- HPE_causing_allele_mouse %>% filter(ContainsHoloprosencephaly == TRUE)



fwrite(HPE_causing_allele_mouse, "HPE_alleles_mouse_cx_cn.csv")


fwrite(indeed_HPE, "HPE_mouse_cx_cn.csv")




unique_strains <- sort(unique(unlist(HPE_causing_allele_mouse$HPE_strains)))

heatmap_matrix <- matrix(NA, nrow = nrow(HPE_causing_allele_mouse), ncol = length(unique_strains))


# Loop through the dataframe
for (i in 1:nrow(HPE_causing_allele_mouse)) {
  # Check if ContainsHoloprosencephaly is TRUE or FALSE
  
  # Get the relevant strains for this row
  row_strains <- HPE_causing_allele_mouse$HPE_strains[[i]]
  
  
  # Loop through unique_strains and populate the matrix
  for (j in 1:length(unique_strains)) {
    if ((unique_strains[j] %in% row_strains) & (HPE_causing_allele_mouse$ContainsHoloprosencephaly[i])){
     
      heatmap_matrix[i, j] <- (ifelse(HPE_causing_allele_mouse$Zygosity[i] == "ht", 3, ifelse(HPE_causing_allele_mouse$Zygosity[i] == "hm", 1, NA)) / length(row_strains))
      
    }
    else if (unique_strains[j] %in% row_strains){
      heatmap_matrix[i, j] <- 0
      
    }
  }
}


# Find rows where all values are either 0 or NA

# Remove the identified rows from the heatmap_matrix
rownames(heatmap_matrix) <- HPE_causing_allele_mouse$OntologyAnnotation.subject.symbol

colnames(heatmap_matrix) <- unique_strains


merge_columns_in_heatmap <- function(heatmap_matrix, desired_col, undesired_columns){
  
  for (col_name in undesired_columns){
    heatmap_matrix[is.na(heatmap_matrix[,desired_col]), desired_col] <- heatmap_matrix[is.na(heatmap_matrix[, desired_col]), col_name]
  }
  

  heatmap_matrix <- heatmap_matrix[, !colnames(heatmap_matrix) %in% undesired_columns]
  return(heatmap_matrix)
}

# "B6J", "C57BL/6Jx129P2/OlaHsdF1", "C57BL/6JCrlj", "C57BL/Gr", "B6N",, 
heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "C57BL", c("B6", "B6JTyr;B6N", "B10Rl", "C57BL/6","C57BL/6J","C57BL/6JAnu",
                                    "C57BL/6N","C57BL/6NCrl","C57BL/6NTac","C57BL/LiA","C57BL10/Rl", "C57BL/6Cr","C57BL/6NHsd","C57BL/6Ola",
                                    "C57BL/10Rl","C57BL/RlGso", "C57BL/10Rl", "B6JJclC3JJcl", "B6Cg", "C57BL/6Ola", "B6SJL"
                                    ))
# "FVB/NRj"

heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "FVB", c("FVB/N", "FVB/NJ"))


# 
heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "101", c("101/Rl", "101/H"))


heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "CD", c("D1"))

heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "DBA", c("DBA/2", "DBA/2J", "D2"))

heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "BALB/c", c("BALB/cJ", "CAnNCrl", "CAnNCrlC3", "BALB/cAnNCrl", "BALB/cByJ", "BALB/cOlaHsd"))

heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "SJL", c("SJL/J"))

#heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "m.domesticusposchiavinus", c("M"))

heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "C3H", c("C3H/He", "C3H/HeH", "C3H/HeJ","C3H/HeN", "C3H/Rl", "C3HeB/FeJ", "C3Hf/HeA", "C3", "C3Sn", "C3H/HeSn"))

# "CBA/Gr", 
heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "CBA", c("CBA/CaH", "CBA/H","CBA/J", "CBA/JNCrlj"))

heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "SEC", c("SEC/RlGso"))

#"Orl:Swiss"
heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "SW", c("SwissWebster", "SWR"))

heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "BlackSwiss", c("NIHBlackSwiss"))

heatmap_matrix <- merge_columns_in_heatmap(heatmap_matrix, "129", c("129/Ola", "129S", "129/Sv","129Cg","129P2", "129P2/OlaHsd", "129S",
                                                                    "129S/Sv", "129S/SvEv","129S1", "129S1/Sv", "129S1/SvImJ","129S2","129S2/SvPas",
                                                                    "129S4","129S4/SvJae","129S4/SvJaeS","129S6","129S6/SvEvTac",
                                                                    "129S7","129S7/SvEvBrd","129X1/Sv", "129X1/SvJ"))



#time to fix some matrix 

#heatmap_matrix[is.na(heatmap_matrix[, 'C57BL/6']), 'C57BL/6'] <- heatmap_matrix[is.na(heatmap_matrix[, 'C57BL/6']), 'B6']

#heatmap_matrix[is.na(heatmap_matrix[, 'C57BL/6J']), 'C57BL/6J'] <- heatmap_matrix[is.na(heatmap_matrix[, 'C57BL/6J']), 'B6J']


#heatmap_matrix[is.na(heatmap_matrix[, 'C57BL/6N']), 'C57BL/6N'] <- heatmap_matrix[is.na(heatmap_matrix[, 'C57BL/6N']), 'B6N']


heatmap_matrix <- heatmap_matrix[, !colnames(heatmap_matrix) %in% c("STOCKShh<tm1Amc>/J", "Trappc10<b2b2416Clo>", "2Gtva>Ptch1<tm1Zim>/Cnrm", 
                                                                    "m.domesticusposchiavinus", "chimera129X1/SvJ", "Cg", "B6SJL", "STOCKGpr37l1<tm1", "STOCKShh<tm1Amc>/J")]


row_names <- rownames(heatmap_matrix)

unique_row_names <- unique(rownames(heatmap_matrix))

# Initialize a new matrix to store the aggregated values
aggregated_matrix <- matrix(NA, nrow = length(unique_row_names), ncol = ncol(heatmap_matrix))
rownames(aggregated_matrix) <- unique_row_names
colnames(aggregated_matrix) <- colnames(heatmap_matrix)

# Loop through unique row names and calculate the average without NAs
for (i in 1:length(unique_row_names)) {
  current_row_name <- unique_row_names[i]
  print(current_row_name)
  rows_to_average <- heatmap_matrix[row_names == current_row_name, ,drop = FALSE]
  print("rows_to_average")
  print(rows_to_average)
  
  if (nrow(rows_to_average) > 1) {
    aggregated_row <- apply(rows_to_average, 2, min, na.rm = TRUE)
    aggregated_matrix[i, ] <- aggregated_row}
  else {
    aggregated_matrix[i, ] <- rows_to_average}
  
}



heatmap_matrix <- aggregated_matrix



#part 2 


rownames(heatmap_matrix)




#heatmap_matrix <- heatmap_matrix[, colSums(heatmap_matrix, na.rm=T)> 0]



heatmap_matrix[is.infinite(heatmap_matrix)] <- NA

remove_second_substring <- function(input) {
  result <- sapply(input, function(label) {
    parts <- unlist(strsplit(label, "<"))
    substring_ <- parts[1]
    print(substring_)
    rest <- paste(parts[-1], collapse = "<")
    print(rest)
    print(grep(substring_, rest))
    if (length(grep(substring_, rest)) > 0) {
      rest <- sub(substring_, "", rest, fixed = TRUE)
      print(rest)
      paste(substring_, rest, sep = "<")
    }
    else {paste(substring_, rest, sep = "<")}
    
  })
  return(result)
}

# Apply the function to the input vector
rownames(heatmap_matrix) <-remove_second_substring(rownames(heatmap_matrix))


substrings_to_remove <- c("Del\\(", "In\\(")

# Identify rows to remove based on row names
rows_to_remove <- rownames(heatmap_matrix) %in% rownames(heatmap_matrix)[grepl(paste(substrings_to_remove, collapse = "|"), rownames(heatmap_matrix))]


test <- heatmap_matrix[!rows_to_remove, ]


test <- test[rowSums(test, na.rm=T)> 0,]


#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#install.packages("grDevices")



my_color_palette <- colorRampPalette(c("white", "#009194"))
# Create a heatmap using ComplexHeatmap


png("test.png",width=25,height=65,units="cm",res=1200)



ht<- Heatmap(test, 
        col=my_color_palette(100),
        name = "HPE Contribution Score",  # Add your own column/row names
        column_title = "Mouse Strain",  # Customize column title
        row_title = "Mouse Models",  # Customize row title
        cluster_rows = F,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 6))

draw(ht,  heatmap_legend_side="top", padding = unit(c(2, 2, 2, 6), "mm"))

dev.off()


#turn_off and on
rownames(test) <- NULL


colPercentages <- function(matrix) {
  num_cols <- ncol(matrix)
  percentages <- numeric(num_cols)
  
  for (i in 1:num_cols) {
    column <- matrix[, i]
    valid_values <- column[!is.na(column)]
    
    str(valid_values)
    str(valid_values[valid_values> 0])
    percentages[i] <- (length(valid_values[valid_values> 0]) / length(valid_values)) * 100
  }
  
  return(percentages)
}



colCounts <- function(matrix) {
  num_cols <- ncol(matrix)
  counts <- numeric(num_cols)
  
  for (i in 1:num_cols) {
    column <- matrix[, i]
    valid_values <- column[!is.na(column)]
    
    str(valid_values)
    str(valid_values[valid_values> 0])
    counts[i] <- length(valid_values[valid_values> 0])
  }
  
  return(counts)
}



just_counts <- colCounts(heatmap_matrix)

# Calculate the percentages for your matrix
percentage_excluding_NA <- colPercentages(heatmap_matrix)

formatted_percentages <- gsub("NaN%", "N/A", paste(round(percentage_excluding_NA, 2), "%", sep=""))


column_annotation <- HeatmapAnnotation(
  custom = anno_text(formatted_percentages, rot=0, gp = gpar(fontsize = 9), just = "center" )
)


column_annotation <- HeatmapAnnotation(
  custom = anno_text(just_counts, rot=0, gp = gpar(fontsize = 9), just = "center" )
)
ht<- Heatmap(heatmap_matrix, 
             col=my_color_palette(100),
             top_annotation = column_annotation,
             name = "HPE Score",  # Add your own column/row names
             column_title = "Mouse Strain",  # Customize column title
             row_title = "",  # Customize row title
             cluster_rows = F,
             cluster_columns = F,
             column_title_side = "bottom",
             row_title_side = "right",
             heatmap_legend_param = list(title_position = "lefttop-rot"),
             row_names_gp = gpar(fontsize = 10), 
             column_names_gp = gpar(fontsize = 12))

draw(ht,  heatmap_legend_side="right", padding = unit(c(3, 1, 3, 5), "mm"))


draw(ht,  heatmap_legend_side="left", padding = unit(c(3, 2, 3, 60), "mm"))
# Add color bar scale




legend <- HeatmapLegend(
  title = "HPE",
  angle = 45  # Set the rotation angle to 45 degrees
)




