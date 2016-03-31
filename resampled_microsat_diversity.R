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
###		your GUI.
### 
### Running the script:
### 
### 	Run the main function resample.ms.div(), for example
### 
### 	x<-resample.ms.div(
###			x = "my_genepop_file.txt",
###			ploidy=2,
###			motif_lengths = c(2,2,2,3),
###			iterations=1000,
###			fn=allelic.richness,
###			min.size=16,
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
### 				   you have AT, AG, AT, AAT, CG it should be 
###				   c(2,2,2,3,2) )
###
### 	iterations: number of iterations. There is a default of 1000.
###
### 	fn: the function used to calculate the diversity index, can be 
### 	    allelic.richness, allelic.range, He (expected heterozigocity) 
### 	    or MSD (mean squared distance).
###
### 	min.size: the minimum size (in alleles not individuals) 
###		  a pop has no have at a locus to be used (else NA is 
### 		  returned), defaults to 8 * ploidy
###
### 	resampling.size: the number of alleles of resampled populations,
### 			 defaults to min.size and is not recommended to
### 			 exceed this value. Only results with the same 
### 			 resampling.size are comparable.
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
# and the size to be resampled to
resample.ms.div = function(
	x, 
	ploidy, 
	motif_lengths, 
	iterations=1000, 
	fn=allelic.richness, 
	min.size=8 * ploidy, 
	resampling.size=min.size
){
	a = scan(x,'character',sep='\n')
	a = a[which(a!='')]
	a[length(a)+1] = 'pop'
	b = grep('Pop', a, ignore.case=T)
	popsizes = rep(0, length(b)-1)
	for (i in 1:length(popsizes)){
		popsizes[i] = b[i+1] - b[i] -1
	}
	pops = length(popsizes)
	popnames = rep('', pops)
	loci = b[1] - 2
	
	alleles<-list()
	for (locus in 1:loci){
		alleles[[locus]] = list()
	}
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
			mat[i,]	= tmp2
		} 
		for(locus in 1:loci){
			if(ploidy == 1) alleles[[locus]][[pop]] = mat[,locus]
			else {
				tmp2 = ''
				field_length = length(strsplit(mat[i,locus], '')[[1]])
				for (i in 1:nrow(mat)){
					all_code_len = round((1/ploidy) * field_length)
					for(j in 1:ploidy){
						indices = ((j-1) * all_code_len + 1):(j * all_code_len)
						tmp2 = c(tmp2, paste(strsplit(mat[i,locus],'')[[1]][indices], collapse=''))
					}
				}
				tmp2 = tmp2[-1]
				tmp2 = as.numeric(tmp2)
				tmp2 = tmp2[tmp2 != 0]
				alleles[[locus]][[pop]] = tmp2
			}
		}
	}

	# resampling
	final.table = matrix(nrow=pops, ncol=loci)
	for (locus in 1:loci){
		for (pop in 1:pops){
			res = vector("numeric",iterations)
			tmp = alleles[[locus]][[pop]]
			if(length(tmp) < min.size) final.table[pop,locus] = NA
			else{
				for (it in 1:iterations){
					res[it] = fn(sample(tmp, resampling.size, T), motif_lengths[locus])
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
	colnames(final.table) = c(a[(1:loci) + 1], 'Mean')
	final.table
}