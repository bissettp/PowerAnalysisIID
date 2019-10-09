####################
## POWER ANALYSIS ##
####################

library(MASS)
library(RColorBrewer)

# variance between subjects and within subjects
between <- 1
within <- 0.5

# construct var-covarmatrix
# the design is 2x4, so 8 measures per subject
measures <- 8
Sigma <- array(within/between,dim=c(measures,measures))
diag(Sigma) <- 1

# prepare array for results
results <- data.frame(subjects=integer(),power_main=double(),power_int=double())

# effect size for both main and interaction effect
cohen <- 0.2

# loop over range of subjects
for(subs in 1:200){

	print(subs)
	mainsig <- c()
	intsig <- c()

	for (i in 1:10000){

		# simulate correlated noise
		noise <- mvrnorm(subs,rep(0,measures),Sigma)
		noise <- as.vector(t(noise))
		x1 <- rep(rep(c(0,1),each=4),subs) #between sub
		x2 <- rep(rep(c(0:3),times=2),subs) #within sub

		# generate data with effects + noise
		y <- cohen*x1 + cohen*x1*ifelse(x2==3,1,0) + noise

		# test for significance with alpha=0.05
		fit <- lm(y~factor(x1)+factor(x2)+factor(x1)*factor(x2))
		test <- anova(fit,test="F")
		mainsig[i] <- test[["Pr(>F)"]][1]<0.05 #True if significant
		intsig[i] <- test[["Pr(>F)"]][3]<0.05 #True if significant

	}

	results[subs,1] <- subs
	results[subs,2] <- mean(mainsig) #Percentage of detected main effects
	results[subs,3] <- mean(intsig) #Percentage of detected interaction effects

}

# write power table
write.table(results,"Power_table.txt",sep="\t",row.names=FALSE)

# minimal sample size for main effect for 1-beta = 0.8 or 0.95
ss80main <- min(results$subjects[results$power_main>0.8],na.rm=TRUE)
ss95main <- min(results$subjects[results$power_main>0.95],na.rm=TRUE)

# make figure
pdf("Power_figure.pdf",width=10,height=7)
	par(mar=c(6,6,2,1),oma=c(0,0,4,0))
	col <- brewer.pal(9,"Set1")
	plot(results$subjects,results$power_main,type="l",lwd=3,col=col[1],ylim=c(0,1),xlab="Subjects",ylab="Average power")
	lines(results$subjects,results$power_int,,lwd=3,col=col[2])
	legend(60,0.6,c("main effect, 2 levels","interaction, 2x4 levels","primary power target (80%)","secondary power target (95%)"),col=col[1:4],lwd=c(3,3,1,1),lty=c(1,1,2,2),bty="n")
	abline(h=0.8,lty=2,col=col[3])
	abline(h=0.95,lty=2,col=col[4])
	lines(x=rep(ss80main,2),y=c(0,0.8),col=col[3],lty=2)
	lines(x=rep(ss95main,2),y=c(0,0.95),col=col[4],lty=2)
	text(x=ss80main,y=-0.01,labels=ss80main,col=col[3])
	text(x=ss95main,y=-0.01,labels=ss95main,col=col[4])
	title("Power curves for within subjects 2x4 full factorial design \n (based on 10000 simulations)",outer=TRUE)
dev.off()
