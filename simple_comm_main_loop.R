### programmed by R. Dale, Jan. 2017
### built from the ESN demonstration code of
# Mantas Lukosevicius 2012, http://minds.jacobs-university.de/mantas

setwd('~/Dropbox/new.projects/multiModal/final-github/emergent-meaning')

seedRand = 62 # one shown in the paper
coefsTAll = c() # coefficients for topic, word, and sound across first 20 PCA components
coefsWAll = c()
coefsSAll = c()

# FOR LOOP FOR 100 SIMULATIONS
seeds = c(round(runif(100)*10000),100)
for (seedRand in seeds) { # so that figures reflect seed 100 (in paper)
  
  print(which(seeds==seedRand))
  set.seed(seedRand)
  
  system(paste('./data/sentgen data/grammarOfConvo.txt -e 10 -s ',seedRand,' -t data/training'))
  system(paste('./data/sentgen data/grammarOfConvo_t1.txt -e 1 -s ',seedRand,' -t data/test1'))
  system(paste('./data/sentgen data/grammarOfConvo_t2.txt -e 1 -s ',seedRand,' -t data/test2'))
  system('tail -22 data/test2.data > temp')
  system('cat data/test1.data temp > data/test.data')
  system('rm -f temp data/test1.data data/test2.data data/*.teach')
  
  source('simple_comm_single_run.R')
  
  coefsTAll = rbind(coefsTAll,coefsT)
  coefsWAll = rbind(coefsWAll,coefsW)
  coefsSAll = rbind(coefsSAll,coefsS)

} # for the seed for loop; comment if doing one seed

pdf(file="figures/pcs_by_predictors_ALL.pdf",height=6,width=5)
plot(colMeans(coefsTAll),lwd=3,col='gray',type='o',xlab='Principal Component',ylab='R-Squared')
text(4,2,'topics',col='gray',cex=1.5)
points(colMeans(coefsWAll),lwd=3,col=rgb(.5,.5,.5),type='o')
text(5.2,1.2,'words',col=rgb(.5,.5,.5),cex=1.5)
points(colMeans(coefsSAll),lwd=3,col='black',type='o')
text(10,.4,'sounds',col='black',cex=1.5)
dev.off()

pdf(file="figures/pcs_by_predictors_ALL.pdf",height=4,width=8)
par(mfrow=c(1,3))
for (i in 1:100) {
  if (i==1) {
    plot(coefsTAll[i,],col='gray',xlab='Principal Component',ylab='R-Squared',main='Topics',ylim=c(0,.8))
  } else {
    points(coefsTAll[i,],col='gray')
  }  
}
for (i in 1:100) {
  if (i==1) {
    plot(coefsWAll[i,],col=rgb(.5,.5,.5),xlab='Principal Component',ylab='R-Squared',main='Words',ylim=c(0,.8))
  } else {
    points(coefsWAll[i,],col=rgb(.5,.5,.5))
  }  
}
for (i in 1:100) {
  if (i==1) {
    plot(coefsSAll[i,],col='black',xlab='Principal Component',ylab='R-Squared',main='Sounds',ylim=c(0,.8))
  } else {
    points(coefsSAll[i,],col='black')
  }  
}
dev.off()




