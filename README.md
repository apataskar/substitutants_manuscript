# substitutants_manuscript
Code for cancer datasets substitutant analysis

Large Scale analysis of CPTAC exemplified by LSCC analysis

1. Datasets for LSCC are downloaded in MZML format from PDC commmons webpage https://proteomic.datacommons.cancer.gov/pdc/study/PDC000234

2. Datasets are formatted per TMT channel (eg, 01CPTAC) as instructed on philosopher webpage https://github.com/Nesvilab/philosopher/wiki/Pipeline-mode-for-TMT-analysis

3. Philospher parameters are set according to Supplementary Tables in the manuscript.

3. Philospher was run;

bin/philosopher pipeline --config params/philosopher.yml 01CPTAC.....
(philosopher.yml parameters are submitted in supplementary tables, Database was created by substituting W-> to all other amino acids

4. ratio_peptide_None.tsv (output from philosopher) is used for further downstream analysis.

Extracted substitutants for all W>X combinations using UNIX command.
grep -Fwf W_F.txt ratio_peptide_None.tsv > W_F2.txtratio_peptide_None.tsv

5. Make a Rscript to get counts;

library ("dplyr")
wf = read.table ("W_F2.txtratio_peptide_None.tsv", head=T)
wy = read.table ("W_Y2.txtratio_peptide_None.tsv", head=T)
wa = read.table ("W_A2.txtratio_peptide_None.tsv", head=T)
wg = read.table ("W_G2.txtratio_peptide_None.tsv", head=T)
wh = read.table ("W_H2.txtratio_peptide_None.tsv", head=T)
wc = read.table ("W_C2.txtratio_peptide_None.tsv", head=T)
wd = read.table ("W_D2.txtratio_peptide_None.tsv", head=T)
we = read.table ("W_E2.txtratio_peptide_None.tsv", head=T)
wi = read.table ("W_I2.txtratio_peptide_None.tsv", head=T)
wl = read.table ("W_L2.txtratio_peptide_None.tsv", head=T)
wm = read.table ("W_M2.txtratio_peptide_None.tsv", head=T)
wn = read.table ("W_N2.txtratio_peptide_None.tsv", head=T)
wp = read.table ("W_P2.txtratio_peptide_None.tsv", head=T)
wq = read.table ("W_Q2.txtratio_peptide_None.tsv", head=T)
ws = read.table ("W_S2.txtratio_peptide_None.tsv", head=T)
wt = read.table ("W_T2.txtratio_peptide_None.tsv", head=T)
wv = read.table ("W_V2.txtratio_peptide_None.tsv", head=T)


wf = wf [,-c(1:6)]
wy = wy [,-c(1:6)]
wa = wa [,-c(1:6)]
wg = wg [,-c(1:6)]
wh = wh [,-c(1:6)]
wc = wc [,-c(1:6)]
wd = wd [,-c(1:6)]
we = we [,-c(1:6)]
wi = wi [,-c(1:6)]
wl = wl [,-c(1:6)]
wm = wm [,-c(1:6)]
wn = wn [,-c(1:6)]
wp = wp [,-c(1:6)]
wq = wq [,-c(1:6)]
ws = ws [,-c(1:6)]
wt = wt [,-c(1:6)]
wv = wv [,-c(1:6)]


wf_c= wf%>% mutate_if (is.numeric, ~1 * (. > 0))
wy_c= wy%>% mutate_if (is.numeric, ~1 * (. > 0))
wa_c= wa%>% mutate_if (is.numeric, ~1 * (. > 0))
wg_c= wg%>% mutate_if (is.numeric, ~1 * (. > 0))
wh_c= wh%>% mutate_if (is.numeric, ~1 * (. > 0))
wc_c= wc%>% mutate_if (is.numeric, ~1 * (. > 0))
wd_c= wd%>% mutate_if (is.numeric, ~1 * (. > 0))
we_c= we%>% mutate_if (is.numeric, ~1 * (. > 0))
wi_c= wi%>% mutate_if (is.numeric, ~1 * (. > 0))
wl_c= wl%>% mutate_if (is.numeric, ~1 * (. > 0))
wm_c= wm%>% mutate_if (is.numeric, ~1 * (. > 0))
wn_c= wn%>% mutate_if (is.numeric, ~1 * (. > 0))
wp_c= wp%>% mutate_if (is.numeric, ~1 * (. > 0))
wq_c= wq%>% mutate_if (is.numeric, ~1 * (. > 0))
ws_c= ws%>% mutate_if (is.numeric, ~1 * (. > 0))
wt_c= wt%>% mutate_if (is.numeric, ~1 * (. > 0))
wv_c= wv%>% mutate_if (is.numeric, ~1 * (. > 0))

wf_c[is.na (wf_c)] <- 0
wy_c[is.na (wy_c)] <- 0
wa_c[is.na (wa_c)] <- 0
wg_c[is.na (wg_c)] <- 0
wh_c[is.na (wh_c)] <- 0
wc_c[is.na (wc_c)] <- 0
wd_c[is.na (wd_c)] <- 0
we_c[is.na (we_c)] <- 0
wi_c[is.na (wi_c)] <- 0
wl_c[is.na (wl_c)] <- 0
wm_c[is.na (wm_c)] <- 0
wn_c[is.na (wn_c)] <- 0
wp_c[is.na (wp_c)] <- 0
wq_c[is.na (wq_c)] <- 0
ws_c[is.na (ws_c)] <- 0
wt_c[is.na (wt_c)] <- 0
wv_c[is.na (wv_c)] <- 0

pdf ("densities_total_number_of_samples_None.pdf")
plot (density (rowSums(wf_c)), col="darkred", lwd=2, ylim=c(0,0.09))
lines (density (rowSums(wy_c)))
lines (density (rowSums(wa_c)))
lines (density (rowSums(wg_c)))
lines (density (rowSums(wh_c)))
lines (density (rowSums(wc_c)))
lines (density (rowSums(we_c)))
lines (density (rowSums(wi_c)))
lines (density (rowSums(wl_c)))
lines (density (rowSums(wm_c)))
lines (density (rowSums(wn_c)))
lines (density (rowSums(wp_c)))
lines (density (rowSums(wq_c)))
lines (density (rowSums(ws_c)))
lines (density (rowSums(wt_c)))
lines (density (rowSums(wv_c)))
lines (density (rowSums(random_c)), col="darkgreen",lwd=2)
dev.off()

#Note in the manuscript
# For figure 3f the filter for wf_c is added, while for figure 3a, this filter is not added;
#wf_c= wy_c [ rowSums (wy_c) < 20,]

# As described in the methods section “analysis of large-scale proteomics data of human cancer”, for intra-tumour type analysis (Fig.3f), a filter for maximum number of samples was applied to retain peptides with higher specificity in expression. However, in the tumor-specific analysis (Fig.3a) this filter was not applied for W>F Substitutants both because of their exclusive significant and specific distribution (Extended Fig. 3, p.val 1.13E-09) and in order to optimize inclusion of signal for gene expression correlation analysis (Fig.3c-e). 
wy_c = wy_c [ rowSums (wy_c) < 20,]
wa_c = wa_c [ rowSums (wa_c) < 20,]
wg_c = wg_c [ rowSums (wg_c) < 20,]
wh_c = wh_c [ rowSums (wh_c) < 20,]
wc_c = wc_c [ rowSums (wc_c) < 20,]
wd_c = wd_c [ rowSums (wd_c) < 20,]
we_c = we_c [ rowSums (we_c) < 20,]
wi_c = wi_c [ rowSums (wi_c) < 20,]
wl_c = wl_c [ rowSums (wl_c) < 20,]
wm_c = wm_c [ rowSums (wm_c) < 20,]
wn_c = wn_c [ rowSums (wn_c) < 20,]
wp_c = wp_c [ rowSums (wp_c) < 20,]
wq_c = wq_c [ rowSums (wq_c) < 20,]
ws_c = ws_c [ rowSums (ws_c) < 20,]
wt_c = wt_c [ rowSums (wt_c) < 20,]
wv_c = wv_c [ rowSums (wv_c) < 20,]

counts_wf = as.data.frame ( colSums (wf_c))
counts_wy = as.data.frame (colSums (wy_c))
counts_wa = as.data.frame (colSums (wa_c))
counts_wg = as.data.frame (colSums (wg_c))
counts_wh = as.data.frame (colSums (wh_c))
counts_wc = as.data.frame (colSums (wc_c))
counts_wd = as.data.frame (colSums (wd_c))
counts_we = as.data.frame (colSums (we_c))
counts_wi = as.data.frame (colSums (wi_c))
counts_wl = as.data.frame (colSums (wl_c))
counts_wm = as.data.frame (colSums (wm_c))
counts_wn = as.data.frame (colSums (wn_c))
counts_wp = as.data.frame (colSums (wp_c))
counts_wq = as.data.frame (colSums (wq_c))
counts_ws = as.data.frame (colSums (ws_c))
counts_wt = as.data.frame (colSums (wt_c))
counts_wv = as.data.frame (colSums (wv_c))

counts_wf$id = rownames  (counts_wf)
counts_wy$id = rownames  (counts_wy)
counts_wa$id = rownames  (counts_wa)
counts_wg$id = rownames  (counts_wg)
counts_wh$id = rownames  (counts_wh)
counts_wc$id = rownames  (counts_wc)
counts_wd$id = rownames  (counts_wd)
counts_we$id = rownames  (counts_we)
counts_wi$id = rownames  (counts_wi)
counts_wl$id = rownames  (counts_wl)
counts_wm$id = rownames  (counts_wm)
counts_wn$id = rownames  (counts_wn)
counts_wp$id = rownames  (counts_wp)
counts_wq$id = rownames  (counts_wq)
counts_ws$id = rownames  (counts_ws)
counts_wt$id = rownames  (counts_wt)
counts_wv$id = rownames  (counts_wv)

colnames (counts_wf)=c("wf","id")
colnames (counts_wy)=c("wy","id")
colnames (counts_wa)=c("wa","id")
colnames (counts_wg)=c("wg","id")
colnames (counts_wh)=c("wh","id")
colnames (counts_wc)=c("wc","id")
colnames (counts_wd)=c("wd","id")
colnames (counts_we)=c("we","id")
colnames (counts_wi)=c("wi","id")
colnames (counts_wl)=c("wl","id")
colnames (counts_wm)=c("wm","id")
colnames (counts_wn)=c("wn","id")
colnames (counts_wp)=c("wp","id")
colnames (counts_wq)=c("wq","id")
colnames (counts_ws)=c("ws","id")
colnames (counts_wt)=c("wt","id")
colnames (counts_wv)=c("wv","id")

counts= merge (counts_wf, counts_wy, by="id")
counts= merge (counts, counts_wa, by="id")
counts= merge (counts, counts_wg, by="id")
counts= merge (counts, counts_wh, by="id")
counts= merge (counts, counts_wc, by="id")
counts= merge (counts, counts_wd, by="id")
counts= merge (counts, counts_we, by="id")
counts= merge (counts, counts_wi, by="id")
counts= merge (counts, counts_wl, by="id")
counts= merge (counts, counts_wm, by="id")
counts= merge (counts, counts_wn, by="id")
counts= merge (counts, counts_wp, by="id")
counts= merge (counts, counts_wq, by="id")
counts= merge (counts, counts_ws, by="id")
counts= merge (counts, counts_wt, by="id")
counts= merge (counts, counts_wv, by="id")

#tumor.txt contains list of CPTAC ids with tumor entries from PDC commons (link above). Same for normal.txt with CPTAC ids for adjacent normal tissue entries.

write.table (counts, "COUNTS_FINAL_None.txt",sep="\t", quote=F)
pdf ("Barplot_counts.pdf")
barplot (sort (colSums (counts [,-1])))
dev.off()

tumor = read.table ("tumor.txt", head=F)
colnames (tumor)="id"

normal = read.table ("normal.txt", head=F)
colnames (normal)="id"

tumor =merge (tumor, counts, by="id")
normal =merge (normal, counts, by="id")

exp = read.table ("T_EXP.txt", head=T)
tumor =merge (tumor, exp, by="id")
normal =merge (normal, exp, by="id")

write.table (tumor,"TUMOR_none.txt",sep="\t", quote=F, row.names=F)
write.table (normal,"Normal_none.txt",sep="\t", quote=F, row.names=F)

pdf ("Boxplot_tumor_normal_NONE.pdf")
par (mfrow=c(1,2))
boxplot (normal$wf, tumor$wf, col=c("darkred","darkgreen"), ylim=c(0,50), main="WF")
boxplot (normal$wy, tumor$wy, col=c("darkred","darkgreen"), ylim=c(0,50), main="WY")
dev.off()





