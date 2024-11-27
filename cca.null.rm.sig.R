cca.null.rm.sig = function(tx,ty,nsamp,nseas,x.list,y.list,rm.sig=FALSE,cca.cor=2,sig=NULL,iseed=1,ntrials=1000,alpha=0.05) {

set.seed(iseed)
neof = max(tx,ty)
ncrit = min(tx,ty)
cca.samp = array(NA,dim=c(ntrials,ncrit))

if (cca.cor[1]==2) cca.cor = rep(0,ncrit)
cca.cor = c(0,cca.cor[-length(cca.cor)])

for (nc in 1:ncrit){
for (nt in 1:ntrials){
    field = array(rnorm(neof*nsamp),dim=c(nsamp,neof))
    n.per.seas = floor(nsamp/nseas)
    for (ns in 1:nseas) for (ne in 1:neof) {
    	si = (ns-1)*n.per.seas
	se = ns*n.per.seas
	if (se>=nsamp) se = nsamp-2
	### REMOVE THE LOCAL SEASONAL MEAN
	if (ns==nseas) {
	   field[(se+1):nsamp,ne] = field[(se+1):nsamp,ne]-mean(field[(se+1):nsamp,ne])
	} else {
    	   field[si:se,ne] = field[si:se,ne]-mean(field[si:se,ne])
    	}
    }
    ### REMOVE AN EXTRA SIGNAL, IF DESIRED
    if (rm.sig){
       for (ne in 1:neof){
       	   field[,ne] = residuals(lm(field[,ne]~sig))
       }
    }
    
    field.y = array(0,dim=c(nsamp,nc))
    for (n in 1:nc){
    	field.y[,n] = cca.cor[n]*field[,n]+rnorm(nsamp)*sqrt(1-cca.cor[n]^2)
    }

    x = field[x.list,1:tx]
    dim(x) = c(length(x.list),tx)
    y = array(0,dim=c(length(x.list),ty))
    for (n in 1:nc) y[,n] = field.y[x.list,n]
    if (nc+1<=ty) y[,(nc+1):ty] = field[y.list,(nc+1):ty]
    dim(y) = c(length(y.list),ty)

    s.xx = cov(x)
    s.yy = cov(y)
    s.xy = cov(x,y)
    s.yx = cov(y,x)

    out.x = gev(s.xy %*% solve(s.yy) %*% s.yx, s.xx)
    cca.samp[nt,nc] = sqrt(out.x$lambda[nc])
    }
}


cca.crit = rep(NA,ncrit)
for (nc in 1:ncrit) cca.crit[nc] = quantile(cca.samp[,nc],probs=1-alpha)

cca.crit

}