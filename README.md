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
ps=phyloseq object
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
########
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
```
##############################
Class-Level Analysis Workflow
##############################

# STEP 1: Agglomerate ASVs/OTUs to Class level
# tax_glom merges all taxa sharing the same Class into single entries
# NArm=TRUE removes taxa with NA at Class level (unclassified)
glom_class <- speedyseq::tax_glom(ps, taxrank = 'Class', NArm = TRUE)

# STEP 2: Transform counts to relative abundance within each sample
# Each sample's Class abundances now sum to 1.0 (100%)
carbom_relative_class <- transform_sample_counts(glom_class, function(x) x / sum(x))

# STEP 3: Melt phyloseq object to long-format data frame for ggplot
# Creates columns: Sample, Description, Class, Abundance (rel. abund. 0-1)
data_class <- psmelt(carbom_relative_class)

# STEP 4: Identify and collapse rare Classes (<0.1 total relative abundance)
# aggregate sums each Class's abundance across ALL samples
abund_class <- aggregate(Abundance ~ Class, data_class, sum)
# Flag Classes whose total contribution across dataset <= 10%
rare_class <- abund_class$Class[abund_class$Abundance <= 0.1]
# Replace all rare Class occurrences with "Others"
data_class$Class[data_class$Class %in% rare_class] <- "Others"

# STEP 5: Calculate mean relative abundance by Description group
# aggregate computes mean abundance per Description-Class combination
data_avg_class <- aggregate(Abundance ~ Description + Class, data_class, mean)
# Normalize within each Description so bars sum to 100%
data_avg_class$rel_abund <- ave(data_avg_class$Abundance, data_avg_class$Description, 
                               FUN = function(x) x/sum(x))

# STEP 6: Reorder Classes by total relative abundance (largest first)
# tapply sums each Class across all Descriptions; sort decreasing
data_avg_class$Class <- factor(data_avg_class$Class, 
                              levels = names(sort(tapply(data_avg_class$rel_abund, 
                                                         data_avg_class$Class, sum), 
                                                  decreasing = TRUE)))

```
```
########
Plotting
########
# STEP 7: Create base stacked bar plot for CLASS
q_class <- ggplot(data_avg_class, aes(x = Description, y = rel_abund, fill = Class))

# STEP 8: Finalize CLASS plot with publication styling
r_16s_class <- q_class +
  geom_bar(stat = "identity", position = "stack") +  # Stack segments to 100%
  scale_fill_manual(values = my_colors) +           # Custom colors per Class  
  labs(fill = NULL) +                              # Remove legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-labels
        legend.position = "bottom", legend.text = element_text(size = 5),
        legend.key.size = unit(0.5, "lines"), legend.spacing.x = unit(0.1, "cm")) +
  guides(fill = guide_legend(reverse = TRUE))       # Reverse legend order

# STEP 9: Apply custom ggplot2 styling (remove gridlines, etc.)
zz_class <- ggplot2.customize(r_16s_class, backgroundColor = "white",
                             removePanelGrid = TRUE, removePanelBorder = TRUE,
                             axisLine = c(0.5, "solid", "black"),
                             legendFontSize = 1, legendKeyWidth = 0.2, legendKeyHeight = 0.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# STEP 10: Arrange and label single CLASS panel
class_plot <- ggarrange(zz_class, labels = c("B"), legend = "right",
                       widths = c(5,5), heights = c(2,2), ncol = 1, nrow = 1)
```
Pannel C
```
##############################
Order-Level Analysis Workflow
##############################

# STEP 1: Agglomerate ASVs/OTUs to Order level
# tax_glom merges all taxa sharing the same Order into single entries
# NArm=TRUE removes taxa with NA at Order level (unclassified)
glom_order <- speedyseq::tax_glom(ps, taxrank = 'Order', NArm = TRUE)

# STEP 2: Transform counts to relative abundance within each sample
# Each sample's Order abundances now sum to 1.0 (100%)
carbom_relative_order <- transform_sample_counts(glom_order, function(x) x / sum(x))

# STEP 3: Melt phyloseq object to long-format data frame for ggplot
# Creates columns: Sample, Description, Order, Abundance (rel. abund. 0-1)
data_order <- psmelt(carbom_relative_order)

# STEP 4: Identify and collapse rare Orders (<0.2 total relative abundance)
# aggregate sums each Order's abundance across ALL samples
abund_order <- aggregate(Abundance ~ Order, data_order, sum)
# Flag Orders whose total contribution across dataset <= 20%
rare_order <- abund_order$Order[abund_order$Abundance <= 0.2]
# Replace all rare Order occurrences with "Others"
data_order$Order[data_order$Order %in% rare_order] <- "Others"

# STEP 5: Calculate mean relative abundance by Description group
# aggregate computes mean abundance per Description-Order combination
data_avg_order <- aggregate(Abundance ~ Description + Order, data_order, mean)
# Normalize within each Description so bars sum to 100%
data_avg_order$rel_abund <- ave(data_avg_order$Abundance, data_avg_order$Description, 
                               FUN = function(x) x/sum(x))

# STEP 6: Reorder Orders by total relative abundance (largest first)
# tapply sums each Order across all Descriptions; sort decreasing
data_avg_order$Order <- factor(data_avg_order$Order, 
                              levels = names(sort(tapply(data_avg_order$rel_abund, 
                                                         data_avg_order$Order, sum), 
                                                  decreasing = TRUE)))
```
```
########
Plotting
########
# STEP 8: Finalize ORDER plot with publication styling
r_16s_order <- q_order +
  geom_bar(stat = "identity", position = "stack") +  # Stack segments to 100%
  scale_fill_manual(values = my_colors) +           # Custom colors per Order  
  labs(fill = NULL) +                              # Remove legend title
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-labels
        legend.position = "bottom", legend.text = element_text(size = 5),
        legend.key.size = unit(0.5, "lines"), legend.spacing.x = unit(0.1, "cm")) +
  guides(fill = guide_legend(reverse = TRUE))       # Reverse legend order

# STEP 9: Apply custom ggplot2 styling (remove gridlines, etc.)
zz_order <- ggplot2.customize(r_16s_order, backgroundColor = "white",
                             removePanelGrid = TRUE, removePanelBorder = TRUE,
                             axisLine = c(0.5, "solid", "black"),
                             legendFontSize = 1, legendKeyWidth = 0.2, legendKeyHeight = 0.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# STEP 10: Arrange and label single ORDER panel
order_plot <- ggarrange(zz_order, labels = c("C"), legend = "right",
                       widths = c(5,5), heights = c(2,2), ncol = 1, nrow = 1)
```
```
#####################################
Combine the 3 pannels in a single plot
#####################################
combined_plot <- ggarrange(zz_phylum, zz_class, zz_order, 
                          labels = c("A", "B", "C"),
                          legend = "none",
                          widths = c(1),
                          heights = c(2,2,2),
                          nrow = 3, ncol = 1,  # Vertical stack: 3 rows, 1 column
                          common.legend = FALSE)

```
