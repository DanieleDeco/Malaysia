# Malaysia


## Figure 5

```
# Load required libraries for phyloseq analysis, plotting, and arrangement
library(phyloseq)      # Core microbiome data structures and analysis
library(speedyseq)     # Fast tax_glom implementation for large datasets
library(ggplot2)       # Publication-quality stacked bar plots
library(ggpubr)        # ggarrange for multi-panel figures
library(tidyr)         # pivot_wider for CSV export
library(scales)        # percent formatting for y-axis
# Define color palette (must exist before plotting)
my_colors <- c("#1f77b4", "Firmicutes" = "#ff7f0e", ...)
```

Pannel A
```
##############################
Phylum-Level Analysis Workflow
################################
# STEP 1: Agglomerate ASVs/OTUs to Phylum level
# tax_glom merges all taxa sharing the same Phylum into single entries
# NArm=TRUE removes taxa with NA at Phylum level (unclassified)
glom_m <- speedyseq::tax_glom(ps, taxrank = 'Phylum', NArm = TRUE)

# STEP 2: Transform counts to relative abundance within each sample
# Each sample's Phylum abundances now sum to 1.0 (100%)
carbom_relative_m <- transform_sample_counts(glom_m, function(x) x / sum(x))

# STEP 3: Melt phyloseq object to long-format data frame for ggplot
# Creates columns: Sample, Description, Phylum, Abundance (rel. abund. 0-1)
data_m <- psmelt(carbom_relative_m)

# STEP 4: Identify and collapse rare Phyla (<0.1 total relative abundance)
# aggregate sums each Phylum's abundance across ALL samples
abund_phylum <- aggregate(Abundance ~ Phylum, data_m, sum)
# Flag Phyla whose total contribution across dataset <= 10%
rare_phylum <- abund_phylum$Phylum[abund_phylum$Abundance <= 0.1]
# Replace all rare Phylum occurrences with "Others"
data_m$Phylum[data_m$Phylum %in% rare_phylum] <- "Others"

# STEP 5: Calculate mean relative abundance by Description group
# aggregate computes mean abundance per Description-Phylum combination
data_avg <- aggregate(Abundance ~ Description + Phylum, data_m, mean)
# Normalize within each Description so bars sum to 100%
data_avg$rel_abund <- ave(data_avg$Abundance, data_avg$Description, 
                         FUN = function(x) x/sum(x))

# STEP 6: Reorder Phyla by total relative abundance (largest first)
# tapply sums each Phylum across all Descriptions; sort decreasing
data_avg$Phylum <- factor(data_avg$Phylum, 
                         levels = names(sort(tapply(data_avg$rel_abund, 
                                                    data_avg$Phylum, sum), 
                                             decreasing = TRUE)))


```
```
#########
Plotting
########
# STEP 7: Create base stacked bar plot
q <- ggplot(data_avg, aes(x = Description, y = rel_abund, fill = Phylum))

# STEP 8: Finalize Phylum plot with publication styling
r_16s_phylum <- q +
  geom_bar(stat = "identity", position = "stack") +  # Stack segments to 100%
  scale_fill_manual(values = my_colors) +           # Custom colors per Phylum
  labs(fill = NULL) +                              # Remove legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-labels
        legend.position = "bottom", legend.text = element_text(size = 5),
        legend.key.size = unit(0.5, "lines"), legend.spacing.x = unit(0.1, "cm")) +
  guides(fill = guide_legend(reverse = TRUE))       # Reverse legend order

# STEP 9: Apply custom ggplot2 styling (remove gridlines, etc.)
zz <- ggplot2.customize(r_16s_phylum, backgroundColor = "white",
                       removePanelGrid = TRUE, removePanelBorder = TRUE,
                       axisLine = c(0.5, "solid", "black"),
                       legendFontSize = 1, legendKeyWidth = 0.2, legendKeyHeight = 0.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# STEP 10: Arrange and label single panel
phylum_plot <- ggarrange(zz, labels = c("A"), legend = "right",
                        widths = c(5,5), heights = c(2,2), ncol = 1, nrow = 1)
```
Pannel B
##############################
Class-Level Analysis Workflow
################################


