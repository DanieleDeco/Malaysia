# Malaysia


## Figure 5 barplot

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
## Figure S PCA
```
# Load required libraries
library(phyloseq)
library(vegan)
library(ggbiplot)
library(easyGgplot2)
```
```
# STEP 1: Normalize OTU abundances per sample to zero mean and unit variance
# decostand "standardize" with MARGIN=2 scales each sample (column) independently
otu_abund_norm <- decostand(otu_table(ps), method="standardize", MARGIN=2, na.rm=TRUE)  # Fixed: na.rm=TRUE (quotes unnecessary)

# STEP 2: Convert to data frame
data_abund_meta <- as.data.frame(otu_abund_norm)

# STEP 3: Remove empty rows (OTUs with all zeros across samples)
data_abund_meta <- data_abund_meta[rowSums(data_abund_meta) > 0, ]

# STEP 4: Transpose: OTUs as rows -> samples as rows (PCA on samples)
data_PCA_otu <- t(data_abund_meta)

# STEP 5: Perform PCA (center=TRUE centers vars, scale.=FALSE keeps standardized scales)
PCA_data_norm_otu <- prcomp(data_PCA_otu, center = TRUE, scale. = FALSE)  # OK: already standardized [web:36]

# STEP 6: Extract metadata from phyloseq
metadata <- sample_data(ps)

# STEP 8: Define grouping and shape variables
groups_name <- metadata$River_1  # Colors by River
habitat <- metadata$Season        # Shapes by Season

# STEP 9: Create PCA biplot
plot_PCA_otu <- ggbiplot::ggbiplot(PCA_data_norm_otu, choices = c(1,2), ellipse = TRUE, circle = FALSE, 
                                   scale = 0, obs.scale = 1, var.axes = FALSE, var.scale = 1, groups = groups_name) +
  geom_point(aes(colour = groups_name, shape = habitat), size = 4) +
  scale_color_manual(values = c("darkorange", "darkblue", "darkgreen")) + 
  scale_fill_manual(values = c(NA, NA, NA)) 

# STEP 10: Customize plot appearance
PCA_OTUs <- ggplot2.customize(plot_PCA_otu,
                              backgroundColor = "white",
                              removePanelGrid = TRUE, removePanelBorder = TRUE, showLegend = TRUE, 
                              axisLine = c(0.5, "solid", "black")) 
```
##Figure S Diversity
```
# Load required libraries
library(phyloseq)     # plot_richness()
library(ggplot2)      # Base plotting 
library(ggpubr)       # ggarrange() 
```
```
#############Diversity Box plot#############
# STEP 1: Subset phyloseq object by Season for separate analysis
# InterMonsoon:
InterMonsoon <- subset_samples(ps, Season %in% c("InterMonsoon"))

# Monsoon: 
Monsoon <- subset_samples(ps, Season %in% c("Monsoon"))

# STEP 2: Create alpha diversity boxplot for InterMonsoon season
# plot_richness computes Observed, Chao1, Shannon indices across Rivers
InterMonsoon <- plot_richness(InterMonsoon,
                   x = "River_1",                    # X-axis: River factor
                   measures = c("Observed",          # Richness (raw OTUs)
                                "Chao1",             # Estimated richness  
                                "Shannon"),          # Diversity (evenness weighted)
                   color = "River_1") +              # Color bars/points by River
  geom_boxplot() +                             # Overlay boxplots on violin/jitter
  theme_bw()                                   # Clean white background theme



# STEP 3: Create alpha diversity boxplot for Monsoon season 
Monsoon <- plot_richness(Monsoon,
                   x = "River_1",
                   measures = c("Observed",
                                "Chao1", 
                                "Shannon"),
                   color = "River_1") +
  geom_boxplot() +
  theme_bw()


# STEP 4: Combine plots vertically with labels
# A = InterMonsoon (top), B = Monsoon (bottom)
ggarrange(InterMonsoon, Monsoon,
          ncol = 1,                         # 1 column
          nrow = 2,                         # 2 rows (vertical)
          labels = c("A", "B"))             # Panel labels
```
