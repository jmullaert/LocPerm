source("functions.R")

# load genotypes
data = read.table("data.012")[,-1]	# Reads genotypes from plink 012 format
geno = aggregate_rareV_CAST(data)

# load phenotypes
pheno = read.table("pheno")$V1

# compute permutation list
pca = read.table("plink.eigenvec")[3:12]	# Reads PC (eigenvecs and eigenvals) from plink output format
eig = read.table("plink.eigenval")$V1[1:10]
permlist = generate_permlist(1000,pca,eig)

# Statistical test procedure
Q = test_stat_LR(pheno,geno)	# observed statistic
Qperm = sapply(1:nrow(permlist), function(j) test_stat_LR(pheno[permlist[j,]],geno)) # list of permuted statistics
get_pvalue_FE(Qperm,Q,bilateral=F) #Full empiric with LR statistic : unilateral p-value

# Other test with z transform
Q2 = test_stat_LRz(pheno,geno)
Qperm2 = sapply(1:nrow(permlist), function(j) test_stat_LRz(pheno[permlist[j,]],geno))
get_pvalue_FE(Qperm2,Q2,bilateral=T) #Full empiric with z-transformed LR statistic : bilateral p-value
get_pvalue_SE(Qperm2,Q2)			 #Semi empiric 