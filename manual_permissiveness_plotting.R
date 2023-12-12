BiocManager::install("ComplexHeatmap", force=T)


library(ComplexHeatmap)


library(readr)

library(dplyr)
library(tidyr)
library(ggplot2)

manual_permissiveness_values <- read_csv("manual_permissiveness_better_major.csv")

head(manual_permissiveness_values)


df <- manual_permissiveness_values

# Define the order you want for 'Human_Freq'
priority_order <- c("Major", "Minor", "Non-Reproducible", "Mouse_Only")

# Create an ordering vector based on the priority order
order_vector <- match(df$Human_Freq, priority_order)

# Order the dataframe by the ordering vector
df_ordered <- df[order(order_vector, df$Gene), ]


test <- as.matrix(df_ordered[, -c(1:4)]) 





#rownames(test) <- sapply(seq_along(df_ordered$Gene), function(i) {
#  gene <- df_ordered$Gene[i]
#  allele <- df_ordered$Allele[i]
#  zygosity <- df_ordered$Zyosity[i]
  #
  #if (zygosity == "hm") {
#    allele_name <- paste(allele, allele, sep = "/")
#  } else {
#    allele_name <- paste(allele, "+", sep = "/")
#  }
  
  # Create an expression with the gene name in italic and the allele name in superscript
#  expression_string <- expression(paste0(gene, "^", allele_name)) 
#  eval(parse(text = expression_string))
##})
##


rownames(test) <-df_ordered$Gene

test[is.na(test)] <- ""

col = c(HPE_Resistant = "#338AB6", HPE_Resistant_Mixed = "#009194", 
        HPE_Permissive = "#B17359", HPE_Permissive_Mixed = "#B96D71",
        HPE_Truncation = "#940300", HPE_Truncation_Mixed="pink",
        Completely_Penetrant_Embryonic_Lethality="#9F59B1", 
        Incompletely_Penetrant_Embryonic_Lethality="#E7D6EB"
        )

df_ordered$Zygosity <- factor(df_ordered$Zyosity, levels = c("ht", "hm"))


df_ordered$Human_Freq <- factor(df_ordered$Human_Freq, levels =c("Major", "Minor", "Non-Reproducible", "Mouse_Only"))

zyg_freq_comb <- interaction(df_ordered$Zygosity, df_ordered$Human_Freq, drop=TRUE)

OP <-oncoPrint(test, 
          alter_fun = list(
            HPE_Resistant = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                 gp = gpar(fill = col["snv"], fill = "#338AB6")),
            HPE_Permissive = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                            gp = gpar(fill = col["snv"], fill = "#B17359")),
            HPE_Truncation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                            gp = gpar(fill = col["snv"], fill = "#940300")),
            
  
            HPE_Resistant_Mixed = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                   gp = gpar(fill = col["indel"], fill = "#009194")),
            
            HPE_Permissive_Mixed = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                                 gp = gpar(fill = col["indel"], fill = "#B96D71")),
            
            
            HPE_Truncation_Mixed = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                                  gp = gpar(fill = col["indel"], fill = "pink")),
            Completely_Penetrant_Embryonic_Lethality=function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                                    gp = gpar(fill = col["indel"], fill = "#9F59B1")),
            
            Incompletely_Penetrant_Embryonic_Lethality=function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                                    gp = gpar(fill = col["indel"], fill = "#E7D6EB"))
            
            ), col = col, show_column_names = TRUE, show_pct = FALSE, row_names_side = "left", column_title_side ="bottom", row_title = "Gene",
          heatmap_legend_param = list(title = "Effect to Holoprosencephaly (HPE)", labels=c("HPE Resistant", "HPE Permissive", "Anterior Truncation", "HPE Resistant (Mixed Background)", 
                                                                                            "HPE Permissive (Mixed Background)", "Completely Penetrant Embryonic Lethality", "Incompletely Penetrant Embryonic Lethality"), 
                                      legend_label_gp = gpar(fontsize = 8)),
          column_names_rot = 45, column_names_gp = gpar(fontsize = 10), column_title = "Genetic Background", row_order = rownames(test), gap = unit(6, "mm"),
          row_names_gp = gpar(fontsize = 10, fontface = "bold.italic"),  column_order = c("C57BL6", "129", "C3H", "FVB", "CD-1","NMRI"), row_split = zyg_freq_comb)

#draw(OP, heatmap_legend_side = "right")


draw(OP, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "cm"))





group_end_positions <- cumsum(table(zyg_freq_comb))


label_positions <- group_end_positions - 0.5

# Normalize the label positions to the range of 0 to 1
total_rows <- nrow(test)
label_positions <- label_positions / total_rows

# Scale and shift the normalized positions to fit within the range of 0.2 to 0.8
label_positions <- (label_positions * 0.7) -0.1

# Invert the positions so that the first label (top of the heatmap) has a higher value than subsequent labels
label_positions <- 1 - label_positions



label_positions <- c(0.89, 0.778, 0.665, 0.578, 0.48)


cate_list <- c("Major Gene - Heterozygous", "Major Gene - Homozygous", "Minor Gene - Homozygous", "Non Reproducible in Humans - Homozygous", "Mouse Models Only - Homozygous")

ht_list <- draw(OP, heatmap_legend_side = "bottom")

# Add text annotations using the appropriate decorate function
for (i in seq_along(label_positions)) {
  # Add the label using grid.text
  grid.text(label = cate_list[i],
            x = unit(30, "lines"),  # Adjust the x position as needed
            y = unit(label_positions[i], "native"),
            just = "center",
            gp = gpar(fontsize = 11, fontface = "bold"),
            vp = viewport(name = "heatmap_body", width = unit(1, "npc"), height = unit(1, "npc")))
}




current_vp <- current.viewport()

for (i in seq_along(label_positions)) {
  # Calculate the y-coordinate in the native space of the oncoprint
  y_coord <- unit(label_positions[i], "native") + unit(2, "lines")  # Adjust this offset as needed
  
  # Add the label using grid.text
  grid.text(label = names(table(zyg_freq_comb))[i],
            x = unit(-3, "lines"),  # Adjust the x position as needed
            y = y_coord,
            just = "right",
            gp = gpar(fontsize = 10, fontface = "bold"),
            vp = current_vp)
}

# Push the viewport to add labels
#pushViewport(oncoprint_vp)

# Add labels using grid.text
#for (i in seq_along(label_positions)) {
#  grid.text(label = names(table(zyg_freq_comb))[i],
#            x = unit(2, "lines"),  # Adjust the x position as needed
#            y = unit(label_positions[i], "native"),
#            just = "right",
#            gp = gpar(fontsize = 10, fontface = "bold"))
#}

# Pop the viewport when done
#popViewport()

decorate_heatmap_body(OP, {
  for (i in seq_along(label_positions)) {
    grid.text(label = names(table(zyg_freq_comb))[i],
              x = unit(-3, "lines"),  # Adjust the x position as needed
              y = unit(label_positions[i], "native"),
              just = "right",
              gp = gpar(fontsize = 10, fontface = "bold"))
  }
}, slice = TRUE)
