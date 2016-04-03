### Copyright 2016 Mauricio Bosshard
### 
### Licensed under the Apache License, Version 2.0 (the "License");
### you may not use this file except in compliance with the License.
### You may obtain a copy of the License at
### 
###     http://www.apache.org/licenses/LICENSE-2.0
### 
### Unless required by applicable law or agreed to in writing, software
### distributed under the License is distributed on an "AS IS" BASIS,
### WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
### See the License for the specific language governing permissions and
### limitations under the License.
### 
### 


### 				INSTRUCTIONS
### 
### Loading the function into your workspace:
### 
### 	Simply copy-paste the contents of this file into the R console.
### 	Anternatively use the source() function or the menu provided by
###	your GUI.
### 
### Running the script:
### 
### 	Run the main function resample.ms.div(), for example
### 
### 	x<-resample.ms.div(
###		x = "example_genepop_file.txt",
###		ploidy=2,
###		motif_lengths = c(2,2,2,3),
###		iterations=1000,
###		fn=allelic.richness,
###		min.size=16,
### 	)
### 	x
### 
### Parameters:
### 
### 	x: a standard genepop file as input
###
### 	ploidy: if data are haploid it should be 1, if diploid 2, etc.
###
### 	motif_lengths: the lengths of the microsatellite motifs (e.g. if
### 		       you have AT, AG, AT, AAT, CG it should be 
###		       c(2,2,2,3,2) )
###
### 	iterations: number of iterations. There is a default of 1000.
###
### 	fn: the function used to calculate the diversity index, can be 
### 	    allelic.richness, allelic.range, He (expected heterozigocity) 
### 	    or MSD (mean squared distance). Advanced users may use a 
###	    custom function modelled after functions below.
###
### 	min.size: the minimum size (in alleles not individuals) 
###		  a pop has to have at a locus to be used (else NA is 
### 		  returned), defaults to 8 * ploidy
###
### 	resampling.size: the number of alleles of resampled populations,
### 			 defaults to min.size and is not recommended to
### 			 exceed this value, since results would still be
###			 biased. Only results with the same resampling 
###			 size are comparable.
###		
###	replace: if True (default) multiple random samples are drawn 
###		 using the frequency distribution from each original 
###		 sample. If False each original sample is reduced to
###		 resampling size multiple times. If unsure, use default.
### 	
### Citation:	
### 	
### 	Bosshard et al. in preparation (I will update this information 
### 	once paper is published)
### 	
	

# the 4 functions below take a vector of allele sizes and the 
# microsatellite motif length as an integer and return their diversity 
# index as a scalar

allelic.richness = function(x, motif_length){
	length(unique(x))
}

allelic.range = function(x, motif_length){
	(range(x)[2]-range(x)[1]) / motif_length
}

MSD = function(x, motif_length){
	Sum = 0
	count = 0
	for (i in 1:(length(x) - 1)){
		for (j in (i+1):length(x)){
			Sum = Sum + ((x[i]-x[j]) / motif_length)^2
			count = count + 1
		}
	}
	Sum / count
}

He = function(x, motif_length){
	props = vector("numeric",0)
	for (i in unique(x)){
		props = cbind(props, sum(x==i) / length(x))
	}
	1 - sum(props^2)
}

# Main function. Takes a genpop file, ploidy as scalar, the 
# microsatellite motif_length lengths as a numeric vector and optionally
# the number of iterations, the function to use, the minimum size 
# (in alleles not individuals) a pop has no have at a locus to be used 
# and the size to be resampled to; returns matrix of diversity indices
# with pops as rows and loci (and mean between loci) as columns

resample.ms.div = function(
	x, 
	ploidy, 
	motif_lengths, 
	iterations=1000, 
	fn=allelic.richness, 
	min.size=8 * ploidy, 
	resampling.size=min.size,
	replace=T
){
	# check if fn behaves correctly
	if(class( try( fn(c(2,2,4), 2) ) ) == "try-error"){
		stop('invalid "fn" argument')
	}
	if(!is.numeric(fn(c(2,2,4), 2))){
		stop('"fn" must return numeric value')
	}
	if(length(fn(c(2,2,4), 2)) != 1){
		stop('"fn" must return scalar')
	}
	
	# read in file and prepare data
	a = scan(x,'character',sep='\n')
	a = a[which(a!='')]
	a[length(a)+1] = 'pop'
	b = grep('Pop', a, ignore.case=T)
	tmp = a
	tmp[b] = gsub(' ', '', tmp[b])
	tmp[b] = gsub('\t', '', tmp[b])
	b = b[nchar(tmp[b]) == 3]
	popsizes = rep(0, length(b)-1)
	for (i in 1:length(popsizes)){
		popsizes[i] = b[i+1] - b[i] -1
	}
	if(!all(popsizes > 0)){
		stop('input file may not include empty pops')
	}
	pops = length(popsizes)
	loci = b[1] - 2
	locus_names = a[(1:loci) + 1]
	locus_names = gsub('\t', '', locus_names)

	
	alleles<-list()
	for (locus in 1:loci){
		alleles[[locus]] = list()
	}
	popnames = rep('', pops)
	for (pop in 1:pops){
		tmp = a[-c(1:b[pop], b[pop+1]:length(a))]
		tmp = strsplit(tmp,',')
		popnames[pop] = as.matrix(data.frame(tmp))[1,1]
		tmp = as.matrix(data.frame(tmp))[2,]
		mat = matrix('',ncol=loci,nrow=length(tmp))
		for (i in 1:length(tmp)){
			tmp2 = gsub('\t',' ',tmp[i])
			tmp2 = strsplit(tmp2,  ' ')[[1]]
			tmp2 = tmp2[tmp2 != '']
			if(length(tmp2) != loci){
				stop(paste('wrong number of loci in pop ',pop,' (ind ',i,')', sep=''))
			}
			mat[i,]	= tmp2
		} 
		for(locus in 1:loci){
			if(ploidy == 1) alleles[[locus]][[pop]] = mat[,locus]
			else {
				tmp = ''
				field_length = length(strsplit(mat[i,locus], '')[[1]])
				for (i in 1:nrow(mat)){
					all_code_len = round((1/ploidy) * field_length)
					for(j in 1:ploidy){
						indices = ((j-1) * all_code_len + 1):(j * all_code_len)
						tmp = c(tmp, paste(strsplit(mat[i,locus],'')[[1]][indices], collapse=''))
					}
				}
				tmp = tmp[-1]
				if(!all(nchar(tmp) == nchar(tmp[1]))){
					stop(paste('wrong allelle code in pop ',pop,', locus ',locus, sep=''))
				}
				tmp = as.numeric(tmp)
				tmp = tmp[tmp != 0]
				alleles[[locus]][[pop]] = tmp
			}
		}
	}
	popnames = gsub('\t', '', popnames)

	# resampling
	final.table = matrix(nrow=pops, ncol=loci)
	for (locus in 1:loci){
		for (pop in 1:pops){
			res = vector("numeric",iterations)
			tmp = alleles[[locus]][[pop]]
			if(length(tmp) < min.size) final.table[pop,locus] = NA
			else{
				for (it in 1:iterations){
					res[it] = fn(sample(tmp, resampling.size, replace), motif_lengths[locus])
				}
				final.table[pop,locus] = mean(res)
			}
			cat(paste("Locus ",locus," Pop ",pop,"\n"))
		}
	}
	
	# final output
	means = matrix(nrow=pops)
	for (pop in 1:pops){
		means[pop] = mean(final.table[pop,], na.rm=T)
	}
	final.table = cbind(final.table, means)
	rownames(final.table) = popnames
	colnames(final.table) = c(locus_names, 'Mean')
	final.table
}
