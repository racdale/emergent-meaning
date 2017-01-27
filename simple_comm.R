### programmed by R. Dale, Jan. 2017
### built from the ESN demonstration code of
# Mantas Lukosevicius 2012, http://minds.jacobs-university.de/mantas

setwd('~/Dropbox/new.projects/multiModal')

# seedRand = 100 # one shown in the paper
coefsTAll = c() # coefficients for topic, word, and sound across first 20 PCA components
coefsWAll = c()
coefsSAll = c()
# 
for (seedRand in c(round(runif(100)*1000),100)) { # so that figures reflect seed 100 (in paper)
# 
print(seedRand)
set.seed(seedRand)

data = read.table('data/node_activations_13overlap.data',sep=',',skip=2)$V1 # uses sentgen output based on logic of convo (topics, words, phonemes)

dataIn = matrix(0,length(data),6)
dataIn[1:length(data)+(data[1:length(data)]-1)*length(data)] = 1
dataOut = dataIn[2:nrow(dataIn),] # it's prediction (then generation), so take one off
dataIn = dataIn[1:(nrow(dataIn)-1),] 

######
# generate the ESN reservoir 
inSize = 6 # 5 "phonemes" and one "end of turn" marker
outSize = 6
resSize = 500
a = 0.3 # leaking rate
l = length(data)-1 # since it's prediction, take one off

Win = matrix(runif(resSize*(1+inSize),-0.5,0.5),resSize) # input to reservoir
W = matrix(runif(resSize*resSize,-0.5,0.5),resSize) # reservoir recurrent connections
rhoW = abs(eigen(W,only.values=TRUE)$values[1]) # singular values of echo state W
W = W * 1.25 / rhoW # set spectral radius of 1 to ensure ES property (Jaeger, 2007)

X = matrix(0,1+inSize+resSize,l) # collect activations on reservoir + input + "on" state
Yt = dataOut[1:l,] # desired output

# run the reservoir with the data and collect X
x = rep(0,resSize)
for (t in 1:(l-1)){
  u = t(dataIn[t,])
  x = (1-a)*x + a*tanh( Win %*% c(1,u) + W %*% x )
  X[,t] = c(1,u,x)
}

# train the output
reg = 1e-8  # regularization coefficient
X_T = t(X)
Wout = t(Yt) %*% X_T %*% solve( X %*% X_T + reg*diag(1+inSize+resSize) )

# let's see how it does
Y = matrix(0,outSize,l)
for (t in 1:l){
  u = dataIn[t,]  
  x = (1-a)*x + a*tanh( Win %*% c(1,u) + W %*% x)
  y = tanh(Wout %*% c(1,u,x))
  Y[,t] = y
}

plot(Y[1,],dataOut[,1]) # predicting topic 1 only phoneme
plot(Y[4,],dataOut[,4]) # consistent phoneme (both topics)
plot(Y[5,],dataOut[,5]) # predicting topic 2 only phoneme
plot(Y[6,],dataOut[,6]) # predicting final word; should be hardest

###
### let's get the well-known test set, going from t1 to t2
data = read.table('data/node_activations_13overlap_test.data',sep=',',skip=2) # uses sentgen output based on logic of convo (topics, words, phonemes)
words = data$V3 # see data file
topics = data$V2

data = data$V1

dataIn = matrix(0,length(data),6)
dataIn[1:length(data)+(data[1:length(data)]-1)*length(data)] = 1
dataOut = dataIn[2:nrow(dataIn),] # it's prediction (then generation), so take one off
dataIn = dataIn[1:(nrow(dataIn)-1),] 
l = dim(dataIn)[1]

# let's see how it moves between "topic spaces"
Y = matrix(0,outSize,l)
xs = c()
for (t in 1:l){
  u = dataIn[t,]  
  x = (1-a)*x + a*tanh( Win %*% c(1,u) + W %*% x)
  xs = cbind(xs,x)
  y = tanh(Wout %*% c(1,u,x))
  Y[,t] = y
}

xsp = prcomp(t(xs)) # find components in activation space used to separate solutions

# find best predicting components for TOPIC
rsqs = c()
for (i in 1:min(resSize,l)) {
  rsqs = c(rsqs,summary(lm(topics[2:(l+1)]==1~xsp$x[,i]))$r.squared)
}
d1 = sort(rsqs,index=T,decreasing=T)$ix[1]
d2 = sort(rsqs,index=T,decreasing=T)$ix[2]

convLabels = c("a","b")
labSeq = convLabels[topics]

pdf(file="figures/convo_topic_transitions.pdf",height=6,width=5)
plotDat = xsp$x
plot(plotDat[,d1],plotDat[,d2],type='b',pch=dataOut+1,col='gray',main='Topic-Predictive Components',
     xlab=paste('Principal Component',d1),ylab=paste('Principal Component',d2))
points(plotDat[1,d1],plotDat[1,d2],type='p',cex=1.5,pch=15,col='green')
points(plotDat[l,d1],plotDat[l,d2],type='p',cex=1.5,pch=15,col='red')
choiceIxes = seq(from=cycleTime,to=nrow(plotDat),by=cycleTime)
text(plotDat[,d1],plotDat[,d2],lab=labSeq[2:(l+1)],cex=1)
dev.off()
td1 = d1 # save for later

# find the best components for WORD 1 vs. 2&3
rsqs = c()
for (i in 1:min(resSize,l)) {
  rsqs = c(rsqs,summary(lm(words[2:(l+1)]==1~xsp$x[,i]))$r.squared)
}
d1 = sort(rsqs,index=T,decreasing=T)$ix[1]
d2 = sort(rsqs,index=T,decreasing=T)$ix[2]

plotDat = xsp$x
plot(plotDat[,d1],plotDat[,d2],type='b',pch=13,cex=.5,col='gray',
     xlab=paste('Principal Component',d1),ylab=paste('Principal Component',d2))
points(plotDat[1,d1],plotDat[1,d2],type='p',cex=1.5,pch=15,col='green')
choiceIxes = seq(from=cycleTime,to=nrow(plotDat),by=cycleTime)
text(plotDat[,d1],plotDat[,d2],lab=words[2:(l+1)],cex=1)
if (d1==td1) {
  wd1 = d2
} else {
  wd1 = d1
}

convLabels = c("i","ii","iii")
labSeq = convLabels[words]

pdf(file="figures/convo_word_and_topic_transitions.pdf",height=6,width=5)
plotDat = xsp$x
plot(plotDat[,td1],plotDat[,wd1],type='b',pch=13,cex=.5,col='gray',main='Topic x Word Component',
     xlab=paste('Topic Principal Component',td1),ylab=paste('Word Principal Component',wd1))
points(plotDat[1,td1],plotDat[1,wd1],type='p',cex=1.5,pch=15,col='green')
points(plotDat[l,td1],plotDat[l,wd1],type='p',cex=1.5,pch=15,col='red')
choiceIxes = seq(from=cycleTime,to=nrow(plotDat),by=cycleTime)
text(plotDat[,td1],plotDat[,wd1],lab=labSeq[2:(l+1)],cex=1)
dev.off()

coefsT = c() # coefficients for topic, word, and sound across first 20 PCA components
coefsW = c()
coefsS = c()
for (i in 1:20) {
  # predict then get factor coefficient
  coefsT = c(coefsT,summary(lm(plotDat[,i]~as.factor(topics[2:(l+1)])))$r.sq)
  coefsW = c(coefsW,summary(lm(plotDat[,i]~as.factor(words[2:(l+1)])))$r.sq)
  coefsS = c(coefsS,summary(lm(plotDat[,i]~as.factor(data[2:(l+1)])))$r.sq)
}
# plot 'em together
pdf(file="figures/pcs_by_predictors.pdf",height=6,width=5)
plot(abs(coefsT),lwd=3,col='gray',type='o',xlab='Principal Component',ylab='R-Squared')
#text(4,2,'topics',col='gray',cex=1.5)
points(abs(coefsW),lwd=3,col=rgb(.5,.5,.5),type='o')
#text(5.2,1.2,'words',col=rgb(.5,.5,.5),cex=1.5)
points(abs(coefsS),lwd=3,col='black',type='o')
#text(10,.4,'sounds',col='black',cex=1.5)
dev.off()

coefsTAll = rbind(coefsTAll,coefsT)
coefsWAll = rbind(coefsWAll,coefsW)
coefsSAll = rbind(coefsSAll,coefsS)

} # for the seed for loop; comment if doing one seed

# if running seed loop, let's do all

#pdf(file="figures/pcs_by_predictors_ALL.pdf",height=6,width=5)
#plot(colMeans(coefsTAll),lwd=3,col='gray',type='o',xlab='Principal Component',ylab='R-Squared')
#text(4,2,'topics',col='gray',cex=1.5)
#points(colMeans(coefsWAll),lwd=3,col=rgb(.5,.5,.5),type='o')
#text(5.2,1.2,'words',col=rgb(.5,.5,.5),cex=1.5)
#points(colMeans(coefsSAll),lwd=3,col='black',type='o')
#text(10,.4,'sounds',col='black',cex=1.5)
#dev.off()



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




