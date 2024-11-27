### THIS FILE TESTS IF OUT OF SAMPLE CORRELATION IS UNDERESTIMATED WHEN 
### OVERFITTING OCCURS
rm(list=ls())

dir.rlib = '/homes/pbuchman/thesis/functions/'
dir.save = '/homes/pbuchman/thesis/cv.CCA/'
dir.cca = '/homes/pbuchman/thesis/second.step/cca.data.era40/'
dir.post02 = '/homes/pbuchman/thesis/post2002/'

source(paste0(dir.rlib,'cca.gev.pc.R'))
source(paste0(dir.rlib,'gev.R'))

#########################################
### CONSTRUCT TIME SERIES
#########################################
nl = 270	#how many time steps in training data
nv = 120	#how many time steps in verification
npc = 20	#how many pcs
sig.l = seq(0,0.5,0.01)
nsig = length(sig.l)
tx = 3
ty = 7
ntrials = 500

print(paste('Tx =',tx))
print(paste('Ty =',ty))

can.cor.mean = rep(NA,nsig)
can.cor.sd = rep(NA,nsig)
oos.mean = rep(NA,nsig)
oos.sd = rep(NA,nsig)

for (sig in sig.l){
can.cor = rep(NA,ntrials)
oos.proj = rep(NA,ntrials)
in.sample = rep(NA,ntrials)
out.sample = rep(NA,ntrials)
set.seed(1)
for (nt in 1:ntrials){
    fx = array(rnorm(npc*(nl+nv)),dim=c(nl+nv,npc))
    fy = array(rnorm(npc*(nl+nv)),dim=c(nl+nv,npc))
    for (n in 1:npc){
    	fx[,n] = fx[,n]/sqrt(var(fx[,n]))
    	fy[,n] = fy[,n]/sqrt(var(fy[,n]))
    }
    fy[,1] = sig*fx[,1] + sqrt(1-sig^2)*rnorm(nl+nv)
    fy[,1] = fy[,1]/sqrt(var(fy[,1]))
    fx.t = fx[1:nl,]
    fy.t = fy[1:nl,]
    fx.oos = fx[(nl+1):(nl+nv),]
    fy.oos = fy[(nl+1):(nl+nv),]
    cca.out = cca.gev.pc(fx.t,fy.t,tx=tx,ty=ty)
    qx = cca.out$qx
    qy = cca.out$qy
    fx.proj = fx.oos[,1:tx] %*% qx
    fy.proj = fy.oos[,1:ty] %*% qy
    can.cor[nt] = cca.out$can.cor[1]
    oos.proj[nt] = cor(fx.proj[,1],fy.proj[,1])
    in.sample[nt] = cor(fx.t[,1],fy.t[,1])
    out.sample[nt] = cor(fx.oos[,1],fy.oos[,1])
}
print(sig)
print(paste(round(mean(can.cor),5),'   ',round(sd(can.cor),5)))
print(paste(round(mean(oos.proj),5),'   ',round(sd(oos.proj),5)))

num = which(sig==sig.l)
can.cor.mean[num] = mean(can.cor)
can.cor.sd[num] = sd(can.cor)
oos.mean[num] = mean(oos.proj)
oos.sd[num] = sd(oos.proj)
}

file.save = paste0(dir.save,'oos.montecarlo.tx',tx,'.ty',ty,'.RData')
oos.list = list(sig=sig.l,can.cor.mean=can.cor.mean,can.cor.sd=can.cor.sd,oos.mean=oos.mean,oos.sd=oos.sd)
save(oos.list,file=file.save)

