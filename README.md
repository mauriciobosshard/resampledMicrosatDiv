Simple script to calculate resampled genetic population diversity from 
microsatellite data.

In population genetics, molecular ecolgy, conservation genetics and other
fields genetic population diversity is useful to assess population 
history, for example past population bottlenecks leave a signature of
reduced population diversity even after recovery of population numbers.
However, estimated diversity depends heavily on sample size, difficulting
direct comparisons between samples. Suggested methods to address this 
problem include random resampling all populations to a specific size 
(e.g. Leberg 2002). 

This script is an implementation of this approach. It calculates diversity 
indices on multiple random samples of a specified size and reports their 
mean for each loci and in total, by population. It takes microsatellite 
data read from a standard genepop file. Available indices include allellic 
richness, allellic range, expected heterozigocity and mean squared 
distance between allelles. Custom functions may also be used. 

References:

Leberg PL 2002 Estimating allelic richness: Effects of sample size and 
bottlenecks. Molecular Ecology 11: 2445–2449
