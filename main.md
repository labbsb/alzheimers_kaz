# Calling genotypes again: 
```bash
echo 'export PATH=/home/aygera/tools/array-analysis-cli-linux-x64-v2.1.0/array-analysis-cli/:$PATH' >> ~/.bashrc
source ~/.bashrc
array-analysis-cli genotype call \
    --num-threads 32 \
    --idat-sample-sheet ./sample_sheet2/sample_sheet2.csv \
    --bpm-manifest ./manifest_bpm/InfiniumImmunoArray-24v2-0_A.bpm    \
    --cluster-file ./cluster/InfiniumImmunoArray-24v2-0_A_ClusterFile.egt \    
    --output-folder ./gtc/

array-analysis-cli genotype gtc-to-vcf \
    --bpm-manifest ./manifest_bpm/InfiniumImmunoArray-24v2-0_A.bpm \
    --genome-fasta-file /home/aygera/biostar/gwas/redo_october/making_vcf/GRCh37_genome/GRCh37_genome.fa \
    --gtc-sample-sheet ./sample_sheet2/sample_sheet2.csv \
    --csv-manifest ./manifest_csv/InfiniumImmunoArray-24v2-0_A.csv \
    --output-folder ./vcf2/

bgzip vcf2/*vcf*
for f in vcf2/*.vcf.gz; do tabix -p vcf -f $f;done
bcftools merge /home/aygera/biostar/gwas/alzheimer/vcf2/*.vcf.gz -o alz.vcf
bcftools reheader -s names.tsv -o alz1.vcf alz.vcf
cd all_samples_plink
plink -vcf ../alz1.vcf --pheno phenotypes.tsv --make-bed --out alz2

awk -F'\t' '!seen[$1]++' InfiniumImmunoArray-24v2-0_A_b138_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > IIA-dictionary.txt
plink --bfile alz2 --update-name IIA-dictionary.txt --make-bed --out alz3
awk '$2 !~ /^rs/' alz3.bim | sort -k2,2 > custom_non_rs_SNP.txt
plink --bfile alz3 --exclude custom_non_rs_SNP.txt --make-bed --out alz4

# Remove individuals with low call rate < 0.9
awk '$1 > 0.9 {print $7, $7}' OFS="\t" callRate.tsv | tr -d '\r' > keep_call_rate.tsv
plink --bfile alz4 --keep keep_call_rate.tsv --make-bed --out alz4.5

plink --bfile alz4.5 --geno 0.05 --make-bed --out alz5
plink --bfile alz5 --mind 0.02 --make-bed --out alz6
plink --bfile alz6 --genome --min 0.2 --out pihat_min0.2
plink --bfile alz6 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt

echo "D180    D180
AK015    AK015
C132    C132
C124    C124
C077    C077
C067    C067
C170    C170
C136    C136
AK045   AK045
C165    C165
A016    A016
A003    A003
AK039   AK039
D136    D136
" > relatives_to_remove.tsv

plink --bfile alz6 --remove relatives_to_remove.tsv --allow-no-sex --make-bed --out alz7
plink --bfile alz7 --maf 0.001 --make-bed --out alz8
```

all ethnicities
```bash
plink2 --bfile alz8 --pca 10 --out alz_pca
plink --bfile alz8 --covar alz_pca.eigenvec --allow-no-sex  --model --out model_mafs_filtered
cat model_mafs_filtered.model | awk '$10 != "NA" && $10 < 1e-3' | sort -gk 9,9
```

only kazakhs:
```bash
cd only_kaz
plink2 --bfile alz8 --keep only_kaz.tsv --out alz_kaz
plink2 --bfile alz_kaz --pca 10 --out alz_pca
plink --bfile alz_kaz --covar alz_pca.eigenvec --allow-no-sex  --model --out model
cat model_mafs_filtered.model | awk '$10 != "NA" && $10 < 1e-3' | sort -gk 9,9
```

Getting QC call rate graphs:
```bash
cat alz_kaz.fam | awk '{print $1}' > samples_after.tsv
cat callRate.tsv | head -n 1 > afterQC.tsv
grep -wf samples_after.tsv callRate.tsv >> afterQC.tsv
```

The graphs for QC before after are in Downlaods folder. they need some work (Specified int he file)

Gene names were identified in https://biit.cs.ut.ee/gprofiler/snpense
