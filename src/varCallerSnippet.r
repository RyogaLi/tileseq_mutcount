
phredToProb <- function(phred) {
	if (phred == "#") return(1)
	phredNum <- as.integer(charToRaw(phred))-33
	10^(-phredNum/10)
}

#calculster posteriors for variant calls, taking into account the size of the cluster
# in which the calls were found
bayesBase <- function(calls, quals, templ, mutRate=0.0025, clusterSize=1) {
	#convert phred quals to error probabilities
	probs <- sapply(quals, function(phred) {
		#multiple qualities for the same base can occur for inserations or deletions
		#they are treated as multiplicative (i.e. the probability of all being errors)
		prod(sapply(toChars(phred), phredToProb))
	})
	#all possible hypothetical bases (=q) to test
	qs <- unique(c(templ,calls))
	#calculate priors for each base given the template and mutation rate
	#note: this could be sped up by keeping a hash/lookup table of all possible values around
	r <- mutRate
	m <- clusterSize
	logPriorOdds <- sapply(qs,
		function(q) if (q==templ) {#WT case
			log((1-r)^m)-log(1-(1-r)^m)
		} else {#Mutant case
			log(1-(1-r)^m)-log(3)-log(1-(1-(1-r)^m)/3)
		}
	)
	#log posterior odds = sum of log likelihood ratios plus prior
	k <- colSums(do.call(rbind,lapply(1:length(calls), function(i) {
		#error probablity
		p <- probs[[i]]
		#for each possible base: calculate the log likelihood ratio (LLR).
		#Given that illumina discretizes the qualities, a lookup table could speed this up too
		setNames(sapply(qs, function(q) {
			#if the basecall supports the hypothetical base
			if (calls[[i]] == q) {
				log(1-p)-log(p)+log(3)
			#otherwise, if the basecall contradicts the hypothetical base
			} else {
				log(p) - log(3) - log(1/3-p/9)
			}
		}),qs)
	}))) + logPriorOdds
	#convert to posterior probabilities
	return(exp(k)/(1+exp(k)))
}

#input "callClusters" is list of lists of dataframes. Each top layer element correponds to a
# read pair, containing a list clusters. The clusters themselves are datafraames that list proximal
# proposed variant calls with columns refbase, refpos, readbase1, qual1, readbase2, qual2
# in case of insertions, refbase is "-" and readbase and qual can contain more than one value,
# in case of deletions, readbase1/readbase2 is "-" and qual1/qual2 would contain the qualities of the neighbouring two base calls.
processClusters <- function(callClusters, mutrate=0.0025) {
	#calculate posterior probabilites for each call in each cluster
	out <- lapply(callClusters, function(clusters) {
		#iterate over each cluster for the current read pair
		lapply(clusters, function(cluster) {
			if (!is.data.frame(cluster)){
				return(cluster)
			}
			#cluster size
			csize <- nrow(cluster)
			#iterate over the proposed nucleotide change in the current cluster
			#  and calculate posteriors
			posteriors <- lapply(1:csize,function(i) with(cluster[i,],{
				#check if any calls are missing or "N"
				na1 <- is.na(readbase1) || readbase1=="N"
				na2 <- is.na(readbase2) || readbase2=="N"
				if (na1 && na2) {
					#if there's no usable data, it's WT
					# return(cbind(cluster,varcall=refbase,post=1))
					setNames(1,refbase)
				} else if (na1) {
					bayesBase(readbase2,qual2,refbase,mutrate,csize)
				} else if (na2) {
					bayesBase(readbase1,qual1,refbase,mutrate,csize)
				} else {
					bayesBase(
						c(readbase1,readbase2),
						c(qual1,qual2),
						refbase,mutrate,csize
					)
				}
			}))
			varcall <- sapply(posteriors,function(x)names(which.max(x)))
			callPosteriors <- mapply(function(posteriors,varcall)posteriors[[varcall]],posteriors,varcall)
			return(cbind(cluster,varcall=varcall,post=callPosteriors))
		})
	})
	return(out)
}
