aggregate_rareV_CAST = function(data){
	apply(data,1,function(x) any(x>0)) 	# CAST statistic
}

generate_permlist = function(N,pca,eig){
	permute = function(r) {
		for (k in 1:length(L)){
			J = which(r[L[k],]==1)
			Z = J[r[L[J],k]==1]
			if (length(Z)==1) next;
			l = sample(Z,1)
			tmp = L[k]
			L[k] <<- L[l]
			L[l] <<- tmp
		}
	}

	rr = Reduce('+', lapply(1:10, function(k) eig[k]*outer(pca[,k],pca[,k],'-')^2))
	r30 = apply(rr,1,function(l) 1*(l<=sort(l)[30]))
	L = 1:nrow(pca)
	replicate(50,permute(r30))	# burn in
	do.call(rbind,lapply(1:N, function(i) {replicate(5,permute(r30)); L}))
}

LL = function(a,b,c,d){	# log-liklihood of logistic model with binary covariate
	r = 0
	if (a>0) r = r + a*log(a)
	if (b>0) r = r + b*log(b)
	if (c>0) r = r + c*log(c)
	if (d>0) r = r + d*log(d)
	if (a+b>0) r = r - (a+c)*log(a+c)
	if (c+d>0) r = r - (b+d)*log(b+d)
	r
}

test_stat_LR = function(pheno,geno){ # Returns LR statistic from the contingency table
	a = sum(geno*pheno)
	b = sum((1-geno)*pheno)
	c = sum(geno*(1-pheno))
	d = sum((1-geno)*(1-pheno))
	LL(a,b,c,d)
}

test_stat_LRz = function(pheno,geno){# Returns the z-transformed LR statistic from the contingency table
	N = length(pheno)
	m = mean(pheno)
	a = sum(geno*pheno)
	b = sum((1-geno)*pheno)
	c = sum(geno*(1-pheno))
	d = sum((1-geno)*(1-pheno))
	ll1 = LL(a,b,c,d)
	sign(a-(a+b)*(a+c)/N)*sqrt(2*(ll1-N*m*log(m/(1-m))-N*log(1-m)))
}

get_pvalue_FE = function(Qperm,Q,bilateral=F){	# Full empiric derivation of p-value with continuity correction
	p1 = mean(Qperm>=Q)
	p2 = mean(Qperm>Q)
	p = runif(1,min=p2,max=p1)
	if(bilateral) p = ifelse(p<0.5,2*p,2*(1-p))
	p
}

get_pvalue_SE = function(Qperm,Q){	# Semi empiric derivation of p-value (assumes Qperm is normally distributed)
	m = mean(Qperm)
	v = var(Qperm)
	pchisq((Q-m)^2/v,df=1,low=F)
}

