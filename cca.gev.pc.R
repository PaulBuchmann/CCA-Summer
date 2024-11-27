cca.gev.pc = function(fx,fy,tx=NULL,ty=NULL) {
########################################################
## PERFORMS CANONICAL CORRELATION ANALYSIS ON X AND Y.
## X AND Y ARE ASSUMED TO BE IN THE FOLLOWING FORMS:
## X = FX %*% EX^T   ### AND ### Y = FY %*% EY^T
## IF TX AND& TY = NULL, THEN BOTH (TX,TY) SELECTED USING MIC
## INPUT:
#    FX: TIME SERIES OF X DATA
#    FY: TIME SERIES OF Y DATA
#    TX: TRUNCATION FOR X (TX <= MX); IF NULL, TX IS SELECTED
#    TY: TRUNCATION FOR Y (TY <= MY); IF NULL, TY IS SELECTED
## OUTPUT LIST:
#    MIC[MX,MY]: MUTUAL INFORMATION CRITERION
#    NMIN[1,2]: VALUES OF TX,TY THAT MINIMIZES MIC
#    CAN.COR[MIN(TX,TY)]: CANONICAL CORRELATIONS
#    RX[NTOT,MIN(TX,TY)]: CANONICAL VARIATES FOR X
#    RY[NTOT,MIN(TX,TY)]: CANONICAL VARIATES FOR Y
#    PX[SX  ,MIN(TX,TY)]: CANONICAL LOADING VECTORS FOR X
#    PY[SY  ,MIN(TX,TY)]: CANONICAL LOADING VECTORS FOR Y
#    QX.TILDE[TX,MIN(TX,TY)]: WEIGHTING VECTORS FOR X-FEATURES
#    QY.TILDE[TY,MIN(TX,TY)]: WEIGHTING VECTORS FOR Y-FEATURES
#    TX, TY: SELECTED VALUES OF TX AND TY 
########################################################

ntot = dim(fx)[1]
mx = dim(fx)[2]
my = dim(fy)[2]
if (ntot != dim(fy)[1]) stop('fx and fy have inconsistent time dimension')

lead.can.cor = array(NA,dim=c(mx,my))
micc  = array(NA,dim=c(mx,my))
for (ny in 1:my) for (nx in 1:mx) if (ntot-nx-ny-2 > 0) {
    penalty = (ntot+1)*((nx+ny)/(ntot-nx-ny-2)-nx/(ntot-nx-2)-ny/(ntot-ny-2))

    fx.mic = fx[,1:nx]
    dim(fx.mic) = c(ntot,nx)
    fy.mic = fy[,1:ny]
    dim(fy.mic) = c(ntot,ny)

    s.xx = cov(fx.mic)
    s.yy = cov(fy.mic)
    s.xy = cov(fx.mic,fy.mic)
    s.yx = cov(fy.mic,fx.mic)
    cors = gev(s.xy %*% solve(s.yy) %*% s.yx, s.xx)$lambda
    lead.can.cor[nx,ny] = cors[1]
    micc[nx,ny] = sum(log(1-cors)) + penalty

}
nmin = which(micc == min(micc,na.rm=TRUE),arr.ind=TRUE)
if (is.null(tx) & is.null(ty)){ 
   tx = nmin[1]
   ty = nmin[2]
}

fx.fin = fx[,1:tx]
dim(fx.fin) = c(ntot,tx)
fy.fin = fy[,1:ty]
dim(fy.fin) = c(ntot,ty)

s.xx = cov(fx.fin)
s.yy = cov(fy.fin)
s.xy = cov(fx.fin,fy.fin)
s.yx = cov(fy.fin,fx.fin)

out.x = gev(s.xy %*% solve(s.yy) %*% s.yx, s.xx)
out.y = gev(s.yx %*% solve(s.xx) %*% s.xy, s.yy)
can.cor = sqrt(out.x$lambda)[1:min(tx,ty)]
qx = out.x$q
qy = out.y$q
rx = fx.fin %*% qx
ry = fy.fin %*% qy
if (any(diag(cor(rx,ry))<0)){
   neg = which(diag(cor(rx,ry))<0)
   for (n in neg){
       qy[,n] = -qy[,n]
   }
   ry = fy.fin %*% qy
}

list(can.cor=can.cor,rx=rx,ry=ry,tx=tx,ty=ty,qx=qx,qy=qy,mic=micc,lead.can.cor=lead.can.cor)
}