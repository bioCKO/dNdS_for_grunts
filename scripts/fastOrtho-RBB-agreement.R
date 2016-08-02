#build a histogram of the data output from cross_check_orthos_fish.py
setwd("~/git_Repositories/dNdS_for_grunts/fastOrtho_Results")
dat = read.table("agreementData.tsv", header = T)
head(dat)


hist(dat$percent, xlab = "Percent Agreement", main = '\n')