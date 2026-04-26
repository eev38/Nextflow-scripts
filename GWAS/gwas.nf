params.SNPs = "~/AT_ITA_imputed_EC.vcf.gz"
params.Pheno = "~/pheno.txt"




// convert vcf to plink bed file - if your file is not imputed can clean missing genotypes here

process VC2Plink {

input:
path(SNPs)

output:
path("AT_ITA_clean.{bed,bim,fam,log,nosex}"), emit: cleanbed 

"""
plink --vcf AT_ITA_imputed_EC.vcf.gz --make-bed --out AT_ITA_clean

"""    
}


// Generate phenotype file for plink

process GenPheno {

input:
path(Pheno)
path(cleanbed)


output:
path("input4GEMMA.{bed,bim,fam,log,nosex}"), emit: gemmaFiles 


"""
plink \
--bfile "AT_ITA_clean" \
--no-parents --allow-no-sex \
--pheno pheno.txt \
--make-bed \
--out input4GEMMA

"""
}


// Create relationship matrix for gemma 

process RelMat {

conda '~/.conda/envs/GWAS/bin/'

input:
path(GemmaFiles)


output:
path("RelMatr.{cXX.txt,log.txt}"), emit: relFile

"""
gemma -bfile "input4GEMMA" -gk 1 -outdir . -o RelMatr
"""

}


// Run GWAS p_wald - MAF = 0.05

process Gwas {

conda '~/.conda/envs/GWAS/bin/'   

input:
path(GemmaFiles)
path(RelFile)

output:
path("GWASResults.lmm.{assoc.txt, log.txt}"),  emit: gResults 

"""
gemma -bfile "input4GEMMA" -k RelMatr.cXX.txt -lmm 1 -maf 0.05 -outdir . -o GWASResults.lmm

"""

}


//Create Manhattan plot, qq plot and extract TOP hits based on FDR q-value < 0.05

process ManPLot {


input:
path(gResults)

output:
file("ManhattanPlot.pdf")
file("qqPlot.pdf")
file("TopResults.csv")

"""
module load anaconda3

set +eu

source activate qqm

set -eu
  

Rscript ~/bin/GenManPlot.R "GWASResults.lmm.assoc.txt"

"""


}



// nextflow workflow (input/output for each process)

workflow{

VC2Plink(params.SNPs)
GenPheno(params.Pheno, VC2Plink.out.cleanbed)
RelMat(GenPheno.out.gemmaFiles)
Gwas(GenPheno.out.gemmaFiles, RelMat.out.relFile)
ManPLot(Gwas.out.gResults)

}




