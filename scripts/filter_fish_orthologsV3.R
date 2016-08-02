#Filter orthologs using pair-wise dS


#set working directory and load nec
setwd('~/git_Repositories/dNdS_for_grunts/dS_filter')
source('~/git_Repositories/dNdS_for_grunts/scripts/dS_filtering_functions.R')
require(plotrix)


#upload clustering packages
library(mixtools)
library(mclust)


#upload the overall ortholog table
odat = read.table("reciprocalOrthos_cov50_pctID75_rep1.0.txt", header = T)
colnames(odat) = c('carb', 'flav', 'macro')
head(odat)

#upload the pairwise dS values for carb
CUT = 0.75
carb = read.table('pair-wise_dNdS_carb.txt', header = T)
species = levels(carb$species)
format.names = c('H.flavolineatum', 'H.macrostomum')
do.bic(carb, format.names)
plot.ds(carb, CUT)
carb.class = filter.orthologs(carb, 1, 2, format.names, 'dS')
carb.class = carb.class[,c(1,2,6)]
colnames(carb.class) = c('carb', 'carb-pair', 'carb.class')
head(carb.class)



#Here cdat will represent the combined data
#This is where we will compile the passing/failing data for each of the 
#three pair-wise comparison datasets
cdat = merge(odat, carb.class, by = 'carb')
head(cdat)
nrow(cdat)
nrow(odat)



#upload the flav pairwise data
flav = read.table('pair-wise_dNdS_flav.txt', header = T)
head(flav)
species = levels(flav$species)
format.names = c('H.carbonarium', 'H.macrostomum')
do.bic(flav, format.names)
plot.ds(flav, CUT)
flav.class = filter.orthologs(flav, 1, 2, format.names, 'dS')
head(flav.class)
flav.class = flav.class[,c(1,2,6)]
colnames(flav.class) = c('flav', 'flav-pair', 'flav.class')
head(flav.class)


#merge up
head(cdat)
cdat = merge(cdat, flav.class, by = 'flav')
nrow(cdat)
head(cdat)


#upload the macro data (for the sake of beauty and symetry)
macro = read.table('pair-wise_dNdS_macro.txt', header = T)
head(macro)
species = levels(macro$species)
format.names = c('H.carbonarium', 'H.flavolineatum')
do.bic(macro, format.names)
plot.ds(macro, CUT)
macro.class = filter.orthologs(macro, 1, 2, format.names, 'dS')
head(macro.class)
macro.class = macro.class[,c(1,2,6)]
colnames(macro.class) = c('macro', 'macro-pair', 'macro.class')
head(macro.class)

#final merge
head(cdat)
cdat = merge(cdat, macro.class, by = 'macro')
nrow(cdat)
head(cdat)

#now we have a big redundant dataframe with all the ortholog groups and the class
#indicating which ones failed the dS test


#build a handle for each orhtologous group as the three contig names together and take the unique of that
o.sets = paste(cdat$carb, paste(cdat$flav, cdat$macro, sep = "_"), sep = "_")
cdat$o.sets = o.sets
head(cdat)
u.sets = unique(o.sets)
length(u.sets)


#this loop will run through the dataframe and output the orthologs, replacing dS failing sequences with NAs
#the final result should be the same length as the number of orthologous groups minus a few that all three failed the dS test
res = data.frame()
out.spp = c('carb', 'flav', 'macro')#note order here must be same as classes within for loop
for (set in u.sets){
	# print(set)
	sub = cdat[cdat$o.sets == set,]
	#if the subset lacks any 3s just add it to the results
	con = append(append(as.vector(sub$carb.class), as.vector(sub$flav.class)), as.vector(sub$macro.class))
	if (!3 %in% con){
		res = rbind(res, sub[1,1:3])
		next
	}
	#if it's all threes just skip it
	if (table(con)['3'] == length(con)){
		next
	}
	classes = sub[,c('carb.class', 'flav.class', 'macro.class')]
	tots = apply(classes, 2, sum)
	out = out.spp[tots == max(tots)]
	if (length(out) == 1){
		sub[1,out] <- NA
		res = rbind(res, sub[1,1:3])
	}
}
head(res)
nrow(res)
res2 = res[,c(3,2,1)] #make sure they are in alphabetical order
head(res2)

#write out the dS filtered results
write.table(res2, 'dS_filtered_orthologs.tsv', sep = "\t", row.names = F, quote = F)

