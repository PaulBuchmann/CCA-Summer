cca.pca = function(x.pca,y.pca,tx=NULL,ty=NULL) {
## PERFORMS CANONICAL CORRELATION ANALYSIS ON X AND Y
## X AND Y ARE ASSUMED TO BE IN THE FOLLOWING FORMS:
## X = FX %*% EX^T ### AND ### Y = FY %*% EY^T
## WHERE FX AND FY HAVE COVARIANCE MATRICES = I
## FOR EXAMPLE: FROM PRINCIPAL COMPONENT ANALYSIS
## IF TX OR TY = NULL, THEN BOTH (TX,TY) SELECTED USING MIC
## INPUT:
#    X.PCA: LIST OUTPUT FROM EOF.LATLON[X-DATA]
#    Y.PCA: LIST OUTPUT FROM EOF.LATLON[Y-DATA]
#    TX: TRUNCATION FOR X (TX <= MX); IF NULL, TX IS SELECTED
#    TY: TRUNCATION FOR Y (TY <= MY); IF NULL, TY IS SELECTED
## OUTPUT:
#    MIC[MX,MY]: MUTUAL INFORMATION CRITERION
#    NMIN[1,2]: VALUES OF TX,TY THAT MINIMIZES MIC
#    CAN.COR[MIN(TX,TY)]: CANONICAL CORRELATIONS
#    RX[NTOT,MIN(TX,TY)]: CANONICAL VARIATES FOR X
#    RY[NTOT,MIN(TX,TY)]: CANONICAL VARIATES FOR Y
#    PX[SX  ,MIN(TX,TY)]: CANONICAL LOADING VECTORS FOR X
#    PY[SX  ,MIN(TX,TY)]: CANONICAL LOADING VECTORS FOR Y
#    QX.TILDE[TX,MIN(TX,TY)]: WEIGHTING VECTORS FOR X-FEATURES
#    QY.TILDE[TY,MIN(TX,TY)]: WEIGHTING VECTORS FOR Y-FEATURES
#    TX,TY: SELECTED VALUES OF TX AND TY

trend = 'detrended'
x = x.pca$pc
y = y.pca$pc
nsamp = dim(x)[1]
nx = dim(x)[2]
ny = dim(y)[2]

#### CALCULATE MIC
micc = mic.cca(x[,1:floor(nx/2)],y[,1:floor(ny/2)])
nmin = which(micc$mic == min(micc$mic,na.rm=TRUE),arr.ind=TRUE)
tx = nmin[1]
ty = nmin[2]
print(paste('Tx is',tx,'and Ty is',ty))

#### COMPUTE CANONICAL CORRELATIONS
sigmax_y = t(x[,1:tx])%*%y[,1:ty]/nsamp
CCs = svd(sigmax_y)$d
print(paste0('Canonical Correlations are ',round(CCs,5)))

#### CALCULATE THE CANONICAL VARIATES
u = svd(sigmax_y)$u
s = svd(sigmax_y)$d
v = svd(sigmax_y)$v

Rx = x[,1:tx]%*%u
Ry = y[,1:ty]%*%v

#### CALCULATE THE CANONICAL LOADING VECTORS
Px = x.pca$eof[,1:tx]%*%u
Py = y.pca$eof[,1:ty]%*%v

pdf(file='DJF_two_week_sep_loading_vectors.pdf')
par(mfrow=c(2,1))
for (nc in 1:min(tx,ty)){
    plot_latlon_v4(x.pca$lon,x.pca$lat,Px[,nc],shrinkdomain=TRUE,nbreaks=60,title.say=paste('DJF 2mT EOFs, Loadings for Component',nc,'Cor=',round(CCs[nc],3)))
    plot_latlon_v4(y.pca$lon,y.pca$lat,Py[,nc],shrinkdomain=TRUE,nbreaks=60,title.say=paste('DJF 2mT EOFs 2 weeks lagged, Loadings for Component',nc))
}
dev.off()

#### CALCULATE FRACTION OF VARIANCE EXPLAINED BY FIRST CRITICAL COMPONENT
var.x.tot = sum(x.pca$sval^2)/(nsamp-1)
exp.var.x = rep(NA,dim=min(tx,ty))
for (n in 1:min(tx,ty)) {
  exp.var.x[n] = sum(Px[!x.pca$lbad,n]^2*x.pca$weight[!x.pca$lbad]^2)/var.x.tot
}

var.y.tot = sum(y.pca$sval^2)/(nsamp-1)
exp.var.y = rep(NA,dim=min(tx,ty))
for (n in 1:min(tx,ty)) {
  exp.var.y[n] = sum(Py[!y.pca$lbad,n]^2*y.pca$weight[!y.pca$lbad]^2)/var.y.tot
}

print(paste('Explained variance from each Px is',round(exp.var.x,4)))
print(paste('Explained variance from each Py is',round(exp.var.y,4)))

#### CALCULATE 5% SIGNIFICANCE OF THE CORRELATIONS
CCs.bp = monte.carlo.cca(tx,ty,nsamp)

print('5% significance level of correlations:')
print(round(CCs.bp$crits,5))


}