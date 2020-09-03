# METHYLATION CODE FOUNDIN

##### Cornelis
##### July 2020

##### using code of => https://github.com/perishky/meffil => big thank you

##### set working dir
```
cd /PATH/TO/METH/FOUNDIN_Methylation/
```

### create variant matrix needed for validation

```
module load plink
plink --bim /PATH/TO/WGS/wgs_hg38_ppmi.july2018_ORIGINAL.bim \
--fam /PATH/TO/WGS/wgs_hg38_ppmi.july2018.fam \
--bed /PATH/TO/WGS/wgs_hg38_ppmi.july2018.bed \
--extract /PATH/TO/METH/FOUNDIN_Methylation/snp-names.txt --recodeA \
--out /PATH/TO/METH/FOUNDIN_Methylation/genotypes --noweb
```

### Optional create map file for later....
```
module load plink
plink --bim /PATH/TO/WGS/wgs_hg38_ppmi.july2018_ORIGINAL.bim \
--fam /PATH/TO/WGS/wgs_hg38_ppmi.july2018.fam \
--bed /PATH/TO/WGS/wgs_hg38_ppmi.july2018.bed \
--extract /PATH/TO/METH/FOUNDIN_Methylation/snp-names.txt --recode \
--out /PATH/TO/METH/FOUNDIN_Methylation/genotypes_ped --noweb
# manually updated names to match METH data
# general FOUNDIN format... => METH_PPMI3459_1662_da65_v1
```

# start R with METH data

```
module load R
R
#install.packages("devtools") # if the devtools package is not installed
library(devtools)
#install_github("perishky/meffil")
library(meffil)

options(mc.cores=10)
# Read in samplesheet
samplesheet <- meffil.read.samplesheet(base=".",pattern="METH_FOUNDIN_sample_sheet_v4.csv")
# old one pre-sample switches/bad samples (see below)
# samplesheet <- meffil.read.samplesheet(base=".",pattern="METH_FOUNDIN_sample_sheet.csv")

# load in data
qc.objects <- meffil.qc(samplesheet, verbose=TRUE)

# save QC matrix
save(qc.objects,file="qc.FOUNDIN_first_pass.Robj")
# note if you need to load again use => load("qc.FOUNDIN_first_pass.Robj")

# set QC parameters...
qc.parameters <- meffil.qc.parameters(
	beadnum.samples.threshold             = 0.1,
	detectionp.samples.threshold          = 0.1,
	detectionp.cpgs.threshold             = 0.1, 
	beadnum.cpgs.threshold                = 0.1,
	sex.outlier.sd                        = 5,
	snp.concordance.threshold             = 0.95,
	sample.genotype.concordance.threshold = 0.8
)

# load in WGS genotypes
plink.files <- "genotypes.raw"
genotypes <- meffil.extract.genotypes(plink.files)

# create summary
qc.summary <- meffil.qc.summary(
	qc.objects,
	parameters = qc.parameters,
	genotypes=genotypes
)

# saving the report
meffil.qc.report(qc.summary, output.file="qc/FINAL_report_FOUNDIN_July2020.html")

# OPTIONAL saving genotypes in matrix...
snp.betas <- meffil.snp.betas(qc.objects)
genotypes <- meffil:::calculate.beta.genotypes(snp.betas)
temp <- t(genotypes)
write.table(temp,file="genotypes_METH.txt", row.names = T, quote = F, sep = "\t")

# Remove outlier samples if necessary
outlier <- qc.summary$bad.samples
table(outlier$issue)
index <- outlier$issue %in% c("Control probe (dye.bias)", 
                              "Methylated vs Unmethylated",
                              "X-Y ratio outlier",
                              "Low bead numbers",
                              "Detection p-value",
                              "Sex mismatch",
                              "Genotype mismatch",
                              "Control probe (bisulfite1)",
                              "Control probe (bisulfite2)")

outlier <- outlier[index,]

length(qc.objects)
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
length(qc.objects)
save(qc.objects,file="qc.objects.clean.Robj")

# rerun QC

qc.summary <- meffil.qc.summary(qc.objects,parameters=qc.parameters,genotypes=genotypes)

meffil.qc.report(qc.summary, output.file="qc/FINAL_report_FOUNDIN_July2020_after_QC.html")

# Perform functional normalisation
# note that this is here done combining day 0 and day 65 
# for certain analysis it might be better to split these and analyse separately

y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename="pc.fit.pdf",height=6,width=6)

# set PC to 10 for now

pc <- 10

# Perform quantile normalization
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)
save(norm.objects,file=paste("norm.obj.pc",pc,".Robj",sep=""))

# Generate normalized probe values
norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)

# Generate normalization report
pcs <- meffil.methylation.pcs(norm.beta,probe.range=20000)
norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs)
meffil.normalization.report(norm.summary, output.file="normalization/FINAL_report.html")

# saving data in normalized matrix
meth <- meffil.normalize.samples(norm.objects)
write.table(meth, file="normalization/FINAL_normalized_FOUNDIN.txt", col.names= TRUE, row.names = TRUE, quote = F, sep = "\t")


# ----->  samples
# |
# |
# |
# v
# probes

####### Get annotation

meffil.list.featuresets()
## [1] "450k"   "common" "epic"
y<-meffil.get.features("epic")
head(y)
write.table(y, file="normalization/EPIC_annotation.txt", col.names= TRUE, row.names = F, quote = F, sep = "\t")
```
####### Liftover annotation

```
cut -f 4,5 EPIC_annotation.txt > temp
cut -f 3 EPIC_annotation.txt > temp2
cut -f 5 EPIC_annotation.txt > temp4
paste temp temp4 temp2 > prep_for_liftover.bed

module load crossmap #(from here => https://pubmed.ncbi.nlm.nih.gov/24351709/)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
crossmap bed hg19ToHg38.over.chain.gz prep_for_liftover.bed > prep_for_liftover_hg38.bed
cut -f 6,7,8,9 prep_for_liftover_hg38.bed > EPIC_annotation_hg38.txt 
```

