# ancestry1KG

Inferring ancestry using the 1000 Genomes high coverage reference and admixture (https://dalexander.github.io/admixture/)
```
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
tar -xvf admixture_linux-1.3.0.tar.gz
export PATH="$HOME/admixture/dist/admixture_linux-1.3.0:$PATH"
```

## Step 1: Merge your plink files with the 1KG reference
Note that the file needs to be QCed and on hg38

```bash
plink --bfile 1KG_high_coverage_20130606_g1k_3202.merged --bmerge $prefix --make-bed --out 1kg_$prefix.merged
plink --bfile 1kg_$prefix.merged --missing --out 1kg_$prefix.merged
echo "SNPs with a missing rate > 5%: $(awk '$5 > 0.05' 1kg_$prefix.merged.lmiss | wc -l)"
echo "excluding SNPs with a missing rate > 5%..."
plink --bfile 1kg_$prefix.merged --geno 0.05 --make-bed --out 1kg_$prefix.geno.05.merged
```

## Step 2: Pruning
The admixture manual recommends to prune the input files. <br>
The first command targets for removal each SNP that has an R2 value of greater than 01 with any other SNP within a 50-SNP sliding window (advanced by 10 SNPs each time).
```bash
plink --bfile 1kg_$prefix.geno.05.merged --indep-pairwise 50 10 0.1 
plink --bfile 1kg_$prefix.geno.05.merged --extract plink.prune.in --make-bed --out 1kg_$prefix.geno.05.merged.pruned
``` 
## Step 3: Generate .pop file
In order to run supervised ancestry inference a .pop file indicating the known ancestry is needed.

```bash
# count ancestries in the reference
echo "Ancestral populations in the reference file: $(awk '{print $1}' 1kg_$prefix.geno.05.merged.pruned.fam | grep -v 'con' | grep -v 'cas' | sort | uniq -c)"
awk '{
    # Exclude FIDs containing "con" or "cas" and output "-"
    if ($1 ~ /con/ || $1 ~ /cas/) {
        print "-";
    } else {
        # Otherwise, output the ancestry name based on the FID
        # In this case, output the FID or whatever ancestry you expect from FID
        # Extract ancestry from the FID (Assuming FID format: "AFR_ACB", "EUR_GBR", etc.)
        print $1;
    }
}' 1kg_$prefix.geno.05.merged.pruned.fam > 1kg_$prefix.geno.05.merged.pruned.pop
```

## Step 4: Run ADMIXTURE
```bash
admixture --supervised --seed 666 -C 10  -j16 1kg_scz_pash1.geno.05.merged.pruned.bed 26
```
