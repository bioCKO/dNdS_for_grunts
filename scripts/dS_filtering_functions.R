#dS_filtering_functions.R
#This contains the functions implemented in filter_fish_orthologsV3.R



#this funciton assigns dS classes for each pairwise comparison
#the two classes with the lowest mean dS are assumed to capture true orthologs
#the third component with the highest mean dS is assumed to capture false orthologs
#this function is run for each of the three sets of pairwise comparisons (comparisons to carb, flav, and macro)
#It outputs figures showing the histogram of dS values with the components traced over it
#These are not terribly informative because of scaling and a few huge dS outliers
#also reports how many ortholog calls were flagged as false (class 3)
filter.orthologs = function(dnds.dat, n.panel.rows, n.panel.columns, names, mut.type){
	diff.list = c()
	prop.removed.list = c()
	par(mfrow = c(n.panel.rows, n.panel.columns))
	filt.dnds.dat = data.frame()
	fail.df = data.frame()
	species = levels(dnds.dat$species)
	par(xpd = T)
	for (i in 1:length(species)){
		print("-----------------------------------")
		sp = species[i]
		sp.name = names[i]
		print(paste('Filtering orthologs for', sp))
		sub = dnds.dat[dnds.dat$species == sp,]
		x = na.omit(sub[,mut.type])
		mod1 = emmix(sub,mut.type,NULL,NULL,NULL,3)
		mus = mod1$mu
		mus = mus[order(mus)]
		d = max(mus) / mus[2]
		diff.list = append(diff.list, d)
		title(sp.name, line = 0.5)
		mod = Mclust(x, G = 3, modelNames = c('V'))
		sub$class = mod$class
		x.pos = min(sub[sub$class == 3, mut.type]) #use 3 normally
		points(x.pos, 0, cex = 2, pch = 17, col = 'black')
		sub2 = sub[sub$class < 3,]
		fail.sub = sub[sub$class == 3,]
		removed = nrow(sub) - nrow(sub2)
		prop.removed.list = append(prop.removed.list, removed/nrow(sub))
		print(paste(removed, "orthologs were flagged as likely paralogs and removed"))
		filt.dnds.dat = rbind(filt.dnds.dat, sub)
		fail.df = rbind(fail.df, fail.sub)
		title(paste(removed, 'orthologs removed'), line = -.5)
		title(paste(round(removed/nrow(sub)*100, 0), '% False Positive', sep = ""), line = -1.5)
	}
	par(mfrow = c(1,1))
	print('##########################################')
	print(paste('Mean increase from 3rd to 2nd component =', mean(diff.list)))
	print(paste('Mean proportion orthologs removed =', mean(prop.removed.list)))
	print('##########################################')
	return(filt.dnds.dat)
}




#funciton for building mixture models using mixtools
emmix = function(df, var, means,sigma,lambda,k){
  x = df[,var]
  mod <- normalmixEM(x, mu = means, sigma = sigma, lambda = lambda, k=k, arbvar=T)
  summary(mod)
  par(xpd = T)
  plot(mod, which = 2, breaks = 30, density = TRUE, cex.axis = 1.4, cex.lab = 1.5, cex.main = 1.5, title = "\n") 
  return(mod)
}


#function for running bic from Mclust
do.bic = function(dnds, format.names){
	species = levels(dnds$species)
	sub1 = dnds[dnds$species == species[1],]
	sub2 = dnds[dnds$species == species[2],]
	mod1 = Mclust(sub1$dS, G = 1:8, modelNames = c('V'))
	mod2 = Mclust(sub2$dS, G = 1:8, modelNames = c('V'))
	par(mfrow = c(1,2))
	plot(mod1, what = c("BIC"))
	title(main = format.names[1])
	plot(mod2, what = c("BIC"))
	title(main = format.names[2])
}


#function to plot the dS distirbution (zoomed in)
#and the component boundaries indicating where failed dS calls are
plot.ds = function(dnds, CUT){
	XLIM = c(0, CUT)
	species = levels(dnds$species)
	sub1 = dnds[dnds$species == species[1],]
	sub2 = dnds[dnds$species == species[2],]
	mod1 = Mclust(sub1$dS, G = 3, modelNames = c('V'))
	mod2 = Mclust(sub2$dS, G = 3, modelNames = c('V'))
	m1f = sub1$dS[mod1$classification == 3]
	m2f = sub2$dS[mod2$classification == 3]
	line1 = min(m1f)
	line2 = min(m2f)
	par(mfrow = c(2,2))
	hist(sub1$dS[sub1$dS < CUT], xlim = XLIM, breaks = 20, xlab = 'dS', main = format.names[1])
	points(x = line1, y = 0, pch = 17)
	hist(sub2$dS[sub2$dS < CUT], xlim = XLIM, breaks = 20, xlab = 'dS', main = format.names[2])
	points(x = line2, y = 0, pch = 17)
	plot(mod1, what = 'classification', xlim = XLIM, axes = F, xlab = 'dS', main = '\n')
	title(main = paste(length(m1f), 'false calls'))
	axis(1)
	# abline(v = line1, lty = 2)
	segments(x0 = line1, x1 = line1, y0 = -.2, y1 = 3.4, lty = 2)
	plot(mod2, what = 'classification', xlim = XLIM, axes = F, xlab = 'dS', main = '\n')
	title(main = paste(length(m2f), 'false calls'))
	segments(x0 = line2, x1 = line2, y0 = -.2, y1 = 3.4, lty = 2)
	axis(1)
}
