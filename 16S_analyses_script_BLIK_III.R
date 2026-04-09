library(phyloseq)
library(decontam)
library(dplyr)
library(vegan)
library(microbiomeMarker)
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)
library(rstatix)
library(car)
library(apeglm)
library(DESeq2)
library(scales)
library(microbiomeMarker)
library(Biostrings)
library(tidyr)

larvae_phyloseq_prefilt = readRDS("phyloseq_larvae_prefilt.rds")

colors = c("FO" = "#008080", "EO" = "#EFDDB5", "IO" = "#8B4513") #Color palette

## 1b. Store ASV sequences in refseq before renaming ----
dna <- DNAStringSet(taxa_names(larvae_phyloseq_prefilt))
names(dna) <- taxa_names(larvae_phyloseq_prefilt)
physeq.larvae.filt <- merge_phyloseq(larvae_phyloseq_prefilt, dna)

## 1c.Then rename ASVs
taxa_names(larvae_phyloseq_prefilt) <- paste0("ASV", seq(ntaxa(larvae_phyloseq_prefilt)))

# 2.Removal of contaminants from the negative controls ----
sample_data(larvae_phyloseq_prefilt)$is.neg = sample_data(larvae_phyloseq_prefilt)$Site %in% c("extraction_negative", "pbs_negative", "pcr_negative")
contamdf.prev.01 = isContaminant(larvae_phyloseq_prefilt, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev.01$contaminant)
phy.larvae.filt = prune_taxa(!contamdf.prev05$contaminant, larvae_phyloseq_prefilt) #remove contaminants from the negatives

ntaxa(phy.larvae.filt) #Check number of taxa before filtering

phy.larvae.filt = phy.larvae.filt %>%
  subset_taxa(
    (is.na(Kingdom) | !grepl("Eukaryota|^NA$", Kingdom, ignore.case = TRUE)) &
      !grepl("Chloroplast|Mitochondria", Order, ignore.case = TRUE) &
      !grepl("Chloroplast|Mitochondria", Family, ignore.case = TRUE)
  ) #Remove ASVs assigned to Chloroplasts, Eukaryota, Mitochondria and NAs

ntaxa(phy.larvae.filt) #Check number of taxa after filtering

phy.larvae.filt = subset_samples(phy.larvae.filt, is.neg !=TRUE) #Remove negative controls
sample_data(phy.larvae.filt)$is.neg = NULL
phy.larvae.filt = prune_taxa(taxa_sums(phy.larvae.filt) > 1, phy.larvae.filt) #Remove singletons
saveRDS(phy.larvae.filt, file = "phyloseq_larvae_filt.rds") #Save filtered phyloseq object

# 3.Normalization ----
## 3a. Normalization with Rarefaction ----
rarecurve(tab, step=50, lwd=2, ylab="ASVs", label=FALSE, xlim=c(0, 10000))
asv.blik3 = phyloseq::otu_table(physeq.blik3.filt) %>% as.matrix %>% t()

depths = sample_sums(phy.larvae.filt)
summary(depths) # Check the minimum number of reads for rarefaction
min(depths)

phy.rare = prune_samples(sample_sums(phy.larvae.filt) >= 5047, phy.larvae.filt)

phy.rare = rarefy_even_depth(phy.rare, sample.size = 5047,#5047
                             rngseed = 42, replace = FALSE) #Rarefy

## 3b. Normalization with CSS ----
otu_tab = otu_table(physeq.larvae.filt) #Transpose ASV table so that taxa are rows
otu_tab_t = t(otu_tab)

#Rewrite ASV table with correct orientation
otu_table(physeq.larvae.filt) <- otu_table(otu_tab_t, taxa_are_rows = TRUE) #Rewrite ASV table with correct orientation
physeq.larvae.filt.css = microbiomeMarker::normalize(physeq.larvae.filt, method = "CSS") #Normalized with CSS for beta-div analyses

# 4. Alpha diversity ----
## 4a. Calculation of alpha diversity through Hill numbers ----
hillq0 = estimate_richness(phy.rare, measures = "Observed")

shannon.hill = exp(vegan::diversity(otu_table(phy.rare), index = "shannon"))
shannon.whole.df = as.data.frame(shannon.hill)

simpson.hill = 1/(1-(vegan::diversity(otu_table(phy.rare), index = "simpson")))
simpson.whole.df = as.data.frame(simpson.hill)

hill.div.df = data.frame(hillq0, shannon.whole.df, simpson.whole.df,
                         sample_data(phy.rare))

hill.lg.long = hill.div.df %>%
  pivot_longer(cols = c(Observed, shannon.hill, simpson.hill),
               names_to = "metric",
               values_to = "value")

## 4b. LM & ANOVA for alpha diversity ----
lm.q0 = lm(Observed ~ Site * Sampling_time, data = hill.div.df)
anova.q0 = anova(lm.q0)

lm.q1 <- lm(shannon.hill ~ Site * Sampling_time, data = hill.div.df)
anova.q1 <- anova(lm.q1)

lm.q2 = lm(simpson.hill ~ Site * Sampling_time, data = hill.div.df)
anova.q2 = anova(lm.q2)

#Check at the results
anova.q0
anova.q1
anova.q2

## 4c. Plotting the results ----

lm_labels <- data.frame(
  metric = factor(c("Observed", "shannon.hill", "simpson.hill"), 
                  levels = c("Observed", "shannon.hill", "simpson.hill")),
  label = c("Site:ST P = 0.567", "Site:ST P = 0.470", "Site:ST P = 0.465"),
  site_p = c("Site P = 0.4390", "Site P = 0.5012", "Site P = 0.6361"),
  time_p = c("Sampling time P = 0.166", "Sampling time P = 0.283", "Sampling time P = 0.236"))

lm_labels$full_label = paste(lm_labels$label, lm_labels$site_p, lm_labels$time_p, sep = "\n")

hill.lg.long$Site <- factor(hill.lg.long$Site, levels = c("EO", "IO", "FO"))

# 2. Generate the Plot
hill.scatter = ggplot(hill.lg.long, aes(x = Sampling_time, y = value, color = Site)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "crossbar", 
               position = position_dodge(width = 0.75),
               width = 0.4, linewidth = 0.3, show.legend = FALSE) +
  geom_text(data = lm_labels, aes(x = 1.5, y = Inf, label = full_label), 
            inherit.aes = FALSE, vjust = 1.5, size = 3.5, fontface = "italic") +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_color_manual(values = colors) +
  labs(x = "Sampling Time",
       y = "Diversity Value") +
  theme_pubr() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

hill.scatter

ggsave(filename = "FigS1_lm_Hill.pdf", plot = hill.scatter, height = 5, width = 7.5)

# 5. Beta Diversity ----
## 5a. Calculations of dissimilarity and PERMANOVA ----
mat = as(otu_table(phy.larvae.filt.css), "matrix")
if (taxa_are_rows(phy.larvae.filt.css)) mat = t(mat)

# 2. Bray-Curtis distance
bc.lg.ott = vegdist(mat, method = "bray")
bc.mat = as.matrix(bc.lg.ott)

# 3. Metadata
meta = data.frame(sample_data(phy.larvae.filt.css))
meta$Site = factor(meta$Site)
meta$Sampling_time = factor(meta$Sampling_time)

# Align distance matrix to metadata
bc.mat = bc.mat[rownames(meta), rownames(meta)]
all(rownames(bc.mat) == rownames(meta))  # must be TRUE

# 4. PERMANOVA
adonis2(bc.mat ~ Site, data = meta, permutations = 999)
adonis2(bc.mat ~ Sampling_time, data = meta, permutations = 999)
adonis2(bc.mat ~ Site * Sampling_time, data = meta, permutations = 999)

# 5. PCoA via phyloseq
ord.lg = ordinate(phy.larvae.filt.css, method = "PCoA", distance = "bray")
scores = as.data.frame(ord.lg$vectors[, 1:2])
colnames(scores) = c("Axis1", "Axis2")
scores$SampleID = rownames(scores)

# 6. Merge scores with metadata
meta$SampleID = rownames(meta)
pcoa.df = merge(scores, meta, by = "SampleID")
pcoa.df$Site = factor(pcoa.df$Site)
pcoa.df$Sampling_time = factor(pcoa.df$Sampling_time)

## 5b. Plot results ----

permanova.label = "PERMANOVA\nSite: R² = 0.087, p = 0.130\nSampling time: R² = 0.160, p < 0.001\nSite × Sampling time: R² = 0.326, p < 0.001"

# 7. Plot
pcoa.lg = ggplot(pcoa.df, aes(x = Axis1, y = Axis2, color = Site)) +
  geom_point(aes(shape = factor(Sampling_time)), size = 4, alpha = 0.8) +
  scale_color_manual(values = c("#EFDDB5","#008080", "#8B4513")) +
  scale_fill_manual(values = c("#EFDDB5","#008080", "#8B4513")) +
  #stat_ellipse(aes(group = Site), linetype = 2) +
  ggside::geom_xsidedensity(aes(fill = Site, group = Site), alpha = 0.3, color = NA) +
  ggside::geom_ysidedensity(aes(fill = Site, group = Site), alpha = 0.3, color = NA) +
  labs(
    x = paste0("PCoA1 (", round(ord.lg$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("PCoA2 (", round(ord.lg$values$Relative_eig[2] * 100, 2), "%)"),
    title = "PCoA"
  ) +
  annotate("text",
           x = Inf, y = Inf,
           label = permanova.label,
           hjust = 1.05, vjust = 1.5,
           size = 9 / .pt,
           family = "sans",
           lineheight = 1.3) +
  ggpubr::theme_classic2() +
  theme(panel.grid = element_blank(), strip.background = element_blank())

pcoa.lg

ggsave(filename = "Fig1A_PCoA_BC_020426_hcl.pdf", plot = pcoa.lg, width = 8, height = 5)

# 6. Beta dispersion ----
## 6a. Beta-disp calculation ----
mat = as(otu_table(phy.larvae.filt.css), "matrix")

if(taxa_are_rows(phy.larvae.filt.css)) {
  mat = t(mat) 
}
bc.lg.ott = vegdist(mat, method = "bray")  # rows = samples, columns = taxa

# 4. Sample data as data.frame
lg.ott.sd = as.matrix.data.frame(sample_data(phy.larvae.filt.css))

lg.ott.sd = as.data.frame(lg.ott.sd)

lg.ott.sd$Site = factor(lg.ott.sd$Site)

bd.ott = betadisper(bc.lg.ott, lg.ott.sd$Site)
df.distances = data.frame(dist_to_centroid = bd.ott$distances)
df.ott.dist = cbind(df.distances, lg.ott.sd)

leveneTest(dist_to_centroid ~ Site, data = df.ott.dist)

df.ott.dist %>%
  group_by(Site) %>%
  shapiro_test(dist_to_centroid)

kruskal = df.ott.dist %>%
  group_by(Sampling_time) %>%
  kruskal_test(dist_to_centroid ~ Site) %>%
  add_significance()

## 6b. Plot results

custom.order = c("EO", "IO", "FO") #Custom order of sites
df.ott.dist$Site <- factor(df.ott.dist$Site, levels = custom.order)

beta.bp = ggboxplot(df.ott.dist,
                    x = "Site",
                    y = "dist_to_centroid",
                    fill = "Site",
                    alpha = 0.7, width = 0.5,
                    color = "grey30") +
  geom_jitter(width=0.2, size=1.5, alpha=0.5) +
  # Match your color values to EO, IO, FO order
  scale_fill_manual(values = c("FO" = "#008080", "EO" = "#EFDDB5", "IO" = "#8B4513")) + 
  stat_pvalue_manual(dunn, hide.ns = TRUE) +
  facet_wrap(~ Sampling_time, scales = "free_x") +
  labs(
    title = "Distance to Centroid") +
  ggpubr::stat_compare_means(
    method = "kruskal",
    label = "p.format",   # custom label
    hide.ns = TRUE
  ) +
  theme_pubr() +
  theme(
    plot.title = element_text(size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    legend.position = "none" # Optional: hide legend since X-axis has labels
  )

beta.bp

ggsave("Fig_S2_dist_to_centroid.pdf", plot = beta.bp, width = 6, height = 4)

# 7. Differential abundance analyses with DESeq2 ----
## 7a. Sub-setting and so ----
#Use the un-normalized filtered dataset!!!!

otl.t1 = phyloseq::subset_samples(phy.larvae.filt, Sampling_time == "T1")
otl.t2 = phyloseq::subset_samples(phy.larvae.filt, Sampling_time == "T2")

otl.t1.filt = prune_taxa(taxa_sums(otl.t1) > 5, otl.t1)
otl.t2.filt = prune_taxa(taxa_sums(otl.t2) > 5, otl.t2)

all(colnames(t(otu_table(otl.t1))) == rownames(sample_data(otl.t1)))
all(colnames(t(otu_table(otl.t2))) == rownames(sample_data(otl.t2)))

## 7b. Running DESeq2 with T1 samples ----
lg.io.fo.t1 = phyloseq::subset_samples(otl.t1, Site != "EO")
lg.io.fo.t1

ds.lg.io.fo.t1 = phyloseq_to_deseq2(lg.io.fo.t1, ~ Site)
ds.lg.io.fo.t1 = estimateSizeFactors(ds.lg.io.fo.t1, type = "poscounts")
ds.lg.io.fo.t1 = DESeq(ds.lg.io.fo.t1, test="Wald", fitType="local")
resultsNames(ds.lg.io.fo.t1 )

# Shrink log2 fold changes
res.IO.FO.t1 <- lfcShrink(ds.lg.io.fo.t1, coef="Site_IO_vs_FO", type="normal")

# Filter significant results
res.IO.FO.t1.sig <- res.IO.FO.t1[!is.na(res.IO.FO.t1$padj) & 
                                   res.IO.FO.t1$padj <= 0.05 & 
                                   abs(res.IO.FO.t1$log2FoldChange) >= 2.5, ]


res.IO.FO.t1.sig.tax = cbind(as(res.IO.FO.t1.sig, "data.frame"),
                             as(tax_table(otl.t1)[rownames(res.IO.FO.t1.sig), ], "matrix"))

write.csv(res.IO.FO.t1.sig.tax, "res_deseq2_lg_io_fo_t1.csv")

lg.eo.fo.t1 = phyloseq::subset_samples(otl.t1, Site != "IO")
lg.eo.fo.t1

ds.lg.eo.fo.t1 = phyloseq_to_deseq2(lg.eo.fo.t1, ~ Site)
ds.lg.eo.fo.t1 = estimateSizeFactors(ds.lg.eo.fo.t1, type = "poscounts")
ds.lg.eo.fo.t1 = DESeq(ds.lg.eo.fo.t1, test="Wald", fitType="local")
resultsNames(ds.lg.eo.fo.t1 )

# Shrink log2 fold changes
res.EO.FO.t1 <- lfcShrink(ds.lg.eo.fo.t1, coef="Site_FO_vs_EO", type="normal")

# Filter significant results
res.EO.FO.t1.sig <- res.EO.FO.t1[!is.na(res.EO.FO.t1$padj) & 
                                   res.EO.FO.t1$padj <= 0.05 & 
                                   abs(res.EO.FO.t1$log2FoldChange) >= 2.5, ]


res.EO.FO.t1.sig.tax = cbind(as(res.EO.FO.t1.sig, "data.frame"),
                             as(tax_table(otl.t1)[rownames(res.EO.FO.t1.sig), ], "matrix"))

write.csv(res.EO.FO.t1.sig.tax, "res_deseq2_lg_fo_eo_t1.csv")

lg.eo.io.t1 = phyloseq::subset_samples(otl.t1, Site != "FO")
lg.eo.io.t1

ds.lg.eo.io.t1 = phyloseq_to_deseq2(lg.eo.io.t1, ~ Site)
ds.lg.eo.fo.t1 = estimateSizeFactors(ds.lg.eo.io.t1, type = "poscounts")
ds.lg.eo.io.t1 = DESeq(ds.lg.eo.io.t1, test="Wald", fitType="local")
resultsNames(ds.lg.eo.io.t1 )

# Shrink log2 fold changes
res.EO.IO.t1 <- lfcShrink(ds.lg.eo.io.t1, coef="Site_IO_vs_EO", type="normal")

# Filter significant results
res.EO.IO.t1.sig <- res.EO.IO.t1[!is.na(res.EO.IO.t1$padj) & 
                                   res.EO.IO.t1$padj <= 0.05 & 
                                   abs(res.EO.IO.t1$log2FoldChange) >= 2.5, ]


res.EO.IO.t1.sig.tax = cbind(as(res.EO.IO.t1.sig, "data.frame"),
                             as(tax_table(otl.t1)[rownames(res.EO.IO.t1.sig), ], "matrix"))

write.csv(res.EO.IO.t1.sig.tax, "res_deseq2_lg_io_eo_t1.csv")

res.IO.FO.t1.sig.tax = read.csv("res_deseq2_lg_io_fo_t1.csv")
res.EO.FO.t1.sig.tax = read.csv("res_deseq2_lg_eo_fo_t1.csv")
res.IO.EO.t1.sig.tax = read.csv("res_deseq2_lg_io_eo_t1.csv")
deseq2.lg.t1 = rbind(res.IO.FO.t1.sig.tax, res.EO.FO.t1.sig.tax, res.IO.EO.t1.sig.tax)

## 7c. Plot results of T1 ----

invert_trans = trans_new("invert", 
                         transform = function(x) -x, 
                         inverse = function(x) -x)

deseq2.lg.t1$Taxa = forcats::fct_reorder(deseq2.lg.t1$Taxa, deseq2.lg.t1$log2FoldChange, .desc = TRUE)

custom_palette = c("EOvsFO" = "#004E69", "IOvsFO" = "#562400", "IOvsEO" = "#E0F2F1")

bb.deseq2.t1 = ggplot(deseq2.lg.t1, aes(x = log2FoldChange, y = Taxa)) + 
  geom_point(aes(fill = comparison, size = padj), shape = 21, alpha = 0.75, color = "transparent") + 
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  theme(
    panel.grid.major = element_line(linetype = 1, color = "grey"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, family = "Arial", size = 11),
    axis.text.y = element_text(face = "italic", family = "Arial", size = 10),
    axis.title = element_text(family = "Arial", size = 11),
    legend.title = element_text(family = "Arial", size = 11),
    legend.text = element_text(family = "Arial", size = 11),
    legend.key = element_rect(fill = "transparent"),
    plot.title = element_text(family = "Arial", size = 11),
    panel.background = element_blank()
  ) +
  ylab("Taxa") +
  xlab("log2FoldChange") +
  ggpubr::theme_pubr() +
  scale_fill_manual(values = custom_palette, drop = FALSE) +
  scale_size_continuous(trans = "reverse", range = c(3,10)) +   # map padj
  guides(
    fill = guide_legend(override.aes = list(size = 5)),  
    size = guide_legend(override.aes = list(shape = 21, fill = "grey", color = "black"))
  ) +
  labs(
    title = "T1",
    size = "Significance (-padj)",
    fill = "Comparison"
  ) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_vline(xintercept = c(-2.5, 2.5), size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 10, size = 0.5, linetype = "solid", color = "black") +
  scale_x_continuous(
    limits = c(-10, 10),
    breaks = c(-10, -7.5,-5,-2.5, 0, 2.5, 5, 7.5, 10)
  ) + coord_fixed(ratio = 1)

bb.deseq2.t1

ggsave(filename = "DESeq2_T1.pdf", plot = bb.deseq2.t1)

## 7d. Running DESeq2 with T2 samples ----

lg.io.fo.t2 = phyloseq::subset_samples(otl.t2, Site != "EO")
lg.io.fo.t2

ds.lg.io.fo.t2 = phyloseq_to_deseq2(lg.io.fo.t2, ~ Site)
ds.lg.io.fo.t2 = estimateSizeFactors(ds.lg.io.fo.t2, type = "poscounts")
ds.lg.io.fo.t2 = DESeq(ds.lg.io.fo.t2, test="Wald", fitType="local")
resultsNames(ds.lg.io.fo.t2 )

# Shrink log2 fold changes
res.IO.FO.t2 <- lfcShrink(ds.lg.io.fo.t2, coef="Site_IO_vs_FO", type="normal")

# Filter significant results
res.IO.FO.t2.sig <- res.IO.FO.t2[!is.na(res.IO.FO.t2$padj) & 
                                   res.IO.FO.t2$padj <= 0.05 & 
                                   abs(res.IO.FO.t2$log2FoldChange) >= 2.5, ]


res.IO.FO.t2.sig.tax = cbind(as(res.IO.FO.t2.sig, "data.frame"),
                             as(tax_table(otl.t2)[rownames(res.IO.FO.t2.sig), ], "matrix"))

write.csv(res.IO.FO.t2.sig.tax, "res_deseq2_lg_io_fo_t2.csv")

lg.eo.fo.t2 = phyloseq::subset_samples(otl.t2, Site != "IO")
lg.eo.fo.t2

ds.lg.eo.fo.t2 = phyloseq_to_deseq2(lg.eo.fo.t2, ~ Site)
ds.lg.eo.fo.t2 = estimateSizeFactors(ds.lg.eo.fo.t2, type = "poscounts")
ds.lg.eo.fo.t2 = DESeq(ds.lg.eo.fo.t2, test="Wald", fitType="local")
resultsNames(ds.lg.eo.fo.t2 )

# Shrink log2 fold changes
res.EO.FO.t2 <- lfcShrink(ds.lg.eo.fo.t2, coef="Site_FO_vs_EO", type="normal")

# Filter significant results
res.EO.FO.t2.sig <- res.EO.FO.t2[!is.na(res.EO.FO.t2$padj) & 
                                   res.EO.FO.t2$padj <= 0.05 & 
                                   abs(res.EO.FO.t2$log2FoldChange) >= 2.5, ]


res.EO.FO.t2.sig.tax = cbind(as(res.EO.FO.t2.sig, "data.frame"),
                             as(tax_table(otl.t2)[rownames(res.EO.FO.t2.sig), ], "matrix"))

write.csv(res.EO.FO.t2.sig.tax, "res_deseq2_lg_fo_eo_t2.csv")

lg.eo.io.t2 = phyloseq::subset_samples(otl.t2, Site != "FO")
lg.eo.io.t2

ds.lg.eo.io.t2 = phyloseq_to_deseq2(lg.eo.io.t2, ~ Site)
ds.lg.eo.io.t2 = DESeq(ds.lg.eo.io.t2, test="Wald", fitType="local")
resultsNames(ds.lg.eo.io.t2 )

# Shrink log2 fold changes
res.EO.IO.t2 <- lfcShrink(ds.lg.eo.io.t2, coef="Site_IO_vs_EO", type="normal")

# Filter significant results
res.EO.IO.t2.sig <- res.EO.IO.t2[!is.na(res.EO.IO.t2$padj) & 
                                   res.EO.IO.t2$padj <= 0.05 & 
                                   abs(res.EO.IO.t2$log2FoldChange) >= 2.5, ]


res.EO.IO.t2.sig.tax = cbind(as(res.EO.IO.t2.sig, "data.frame"),
                             as(tax_table(otl.t2)[rownames(res.EO.IO.t2.sig), ], "matrix"))

write.csv(res.EO.IO.t2.sig.tax, "res_deseq2_lg_io_eo_t2.csv")

res.IO.FO.t2.sig.tax = read.csv("res_deseq2_lg_io_fo_t2.csv")
res.EO.FO.t2.sig.tax = read.csv("res_deseq2_lg_fo_eo_t2.csv")
res.IO.EO.t2.sig.tax = read.csv("res_deseq2_lg_io_eo_t2.csv")

deseq2.lg.t2 = rbind(res.IO.FO.t2.sig.tax, res.EO.FO.t2.sig.tax, res.IO.EO.t2.sig.tax)

## 7e. Plot results of T2 ----
invert_trans = trans_new("invert", 
                         transform = function(x) -x, 
                         inverse = function(x) -x)

deseq2.lg.t2$Taxa = forcats::fct_reorder(deseq2.lg.t2$Taxa, deseq2.lg.t2$log2FoldChange, .desc = TRUE)

bb.deseq2.t2 = ggplot(deseq2.lg.t2, aes(x = log2FoldChange, y = Taxa)) + 
  geom_point(aes(fill = comparison, size = padj), shape = 21, alpha = 0.75, color = "transparent") +
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  theme(
    panel.grid.major = element_line(linetype = 1, color = "grey"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, family = "Arial", size = 11, face = "italic"),
    axis.text.y = element_text(family = "Arial", size = 11),   # italic is now handled by scale_y_discrete
    axis.title = element_text(family = "Arial", size = 11),
    legend.title = element_text(family = "Arial", size = 11),
    legend.text = element_text(family = "Arial", size = 11),
    legend.key = element_rect(fill = "transparent"),
    plot.title = element_text(family = "Arial", size = 11),
    panel.background = element_blank()
  ) +
  ylab("Taxa") +
  xlab("log2FoldChange") +
  ggpubr::theme_pubr() +
  scale_fill_manual(values = custom_palette, drop = FALSE) +
  scale_size_continuous(trans = "reverse", range = c(3,10)) +
  guides(
    fill = guide_legend(override.aes = list(size = 5)),  
    size = guide_legend(override.aes = list(shape = 21, fill = "grey", color = "black"))
  ) +
  labs(
    title = "T2",
    size = "Significance (-padj)",
    fill = "Comparison"
  ) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_vline(xintercept = c(-2.5, 2.5), size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 10, size = 0.5, linetype = "solid", color = "black") +
  scale_x_continuous(
    limits = c(-10, 10),
    breaks = c(-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10)
  ) +
  coord_fixed(ratio = 1)

bb.deseq2.t2

ggsave(filename = "DESeq2_T2.pdf", plot = bb.deseq2.t2)

## 7f. Merging DESeq2 plots -----

both.deseq2 = ggpubr::ggarrange(
  bb.deseq2.t1, bb.deseq2.t2,
  ncol = 2, labels = c("A", "B"),
  align = "hv", common.legend = TRUE
)
both.deseq2

ggsave(filename = "Fig_1B_DESeq2.pdf", plot = both.deseq2, width = 12, height = 10)