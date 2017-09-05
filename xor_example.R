### programmed by R. Dale, Jan. 2017
### built from the ESN demonstration code of
# Mantas Lukosevicius 2012, http://minds.jacobs-university.de/mantas

setwd('~/Dropbox/new.projects/multiModal')

data = as.matrix(read.table('data/xor_input.txt')) # load sequences of XOR input output
dataIn = data[,1:2] # input (p1, p2)
dataOut = data[,3] # truth/falsity of XOR operator
l = dim(data)[1] # length of examples in text file

# generate the ESN reservoir
# global parameters of the model (terminology: Lukosevicius, 2012)
inSize = 2
outSize = 1
resSize = 20
a = 0.3 # leaking rate
cycleTime = 20 # duration to "settle" the ESN on an XOR pattern

set.seed(69)

Win = matrix(runif(resSize*(1+inSize),-0.5,0.5),resSize) # input to reservoir
W = matrix(runif(resSize*resSize,-0.5,0.5),resSize) # reservoir recurrent connections
rhoW = abs(eigen(W,only.values=TRUE)$values[1]) # singular values of echo state W
W = W * 1.25 / rhoW # set spectral radius of 1 to ensure ES property (Jaeger, 2007)

X = matrix(0,1+inSize+resSize,l) # collect activations on reservoir + input + "on" state
Yt = dataOut[1:l] # desired output

# run the reservoir with the data and collect X
x = rep(0,resSize)
for (t in 1:l){
  u = t(dataIn[t,])
  for (i in 1:cycleTime) {
    x = (1-a)*x + a*tanh( Win %*% c(1,u) + W %*% x )
  }
  X[,t] = c(1,u,x)
}

# train the output
# compute regression coefficients with reservoir activations and desired output
reg = 1e-8  # regularization coefficient
X_T = t(X)
Wout = t(Yt) %*% X_T %*% solve( X %*% X_T + reg*diag(1+inSize+resSize) )

# show trajectory of the xor solution; test with patterns again
Y = matrix(0,outSize,l)
xs = c()
for (t in 1:l){
  u = dataIn[t,]  
  for (i in 1:cycleTime) {
    x = (1-a)*x + a*tanh( Win %*% c(1,u) + W %*% x )
    xs = cbind(xs,x)
  }
  y = Wout %*% c(1,u,x)
  Y[,t] = tanh(y)
  # generative mode:
  u = y
  ## this would be a predictive mode:
  #u = data[l+t+1] 
}

cor.test(t(Y),dataOut[1:l]) # correlation with output expected
plot(t(Y),dataOut[1:l]) # plot the expected / observed

xsp = princomp(t(xs)) # find components in activation space used to separate solutions

# find best predicting components for this run
choiceIxes = seq(from=cycleTime,to=nrow(xsp$scores),by=cycleTime)
rsqs = c()
for (i in 1:resSize) {
  rsqs = c(rsqs,summary(lm(dataOut~xsp$scores[choiceIxes,i]))$r.squared)
}
d1 = sort(rsqs,index=T,decreasing=T)$ix[1]
d2 = sort(rsqs,index=T,decreasing=T)$ix[2]

# show some example trajectories in PC space that found best solution
#pdf(file='figures/xor_transitions.pdf',height=6,width=5)
cShown = 20
plotDat = xsp$scores[1:(cycleTime*cShown),]
plot(plotDat[,d1],plotDat[,d2],type='b',pch=13,cex=.5,col='gray',main='XOR Operator Components',
     xlab=paste('Principal Component',d1),ylab=paste('Principal Component',d2))
#points(plotDat[1,d1],plotDat[1,d2],type='p',cex=1.5,pch=15,col='green')
choiceIxes = seq(from=cycleTime,to=nrow(plotDat),by=cycleTime)
text(plotDat[choiceIxes,d1],plotDat[choiceIxes,d2],lab=dataOut[1:cShown],cex=2)
#dev.off()
#########



