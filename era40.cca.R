### THIS FILE TESTS IF THERE IS A CORRELATION FROM ONE 2-WEEK MEAN TO THE NEXT
### THIS IS ONLY FOR ERA40 EOFs

rm(list=ls())

dir.save = '/homes/pbuchman/thesis/second.step/cca.data.era40/'
dir.data = '/project/statop/pbuchman/thesis.data/era40/daily/eofs/'
dir.rlib = '/homes/pbuchman/thesis/functions/'

source(paste0(dir.rlib,'cca.gev.R'))
source(paste0(dir.rlib,'gev.R'))
source(paste0(dir.rlib,'cca.null.rm.sig.R'))
source(paste0(dir.rlib,'rm.seas.cycle.R'))
source(paste0(dir.rlib,'rm.seas.mean.R'))

library(fields)
library(maps)
library(rworldmap)
library(ncdf4)
library(caret)
library(lubridate)

load(paste0(dir.data,'t2m.era40.eof.ave.14day.RData'))
#eof.14dayave
ndate = length(eof.14dayave$date.list)
date.list = eof.14dayave$date.list[1:(ndate-13)]

seas.mon = NULL
seas.mon[[1]] = c(12,1,2)
seas.mon[[2]] = c(3,4,5)
seas.mon[[3]] = c(6,7,8)
seas.mon[[4]] = c(9,10,11)
seas.name = c('djf','mam','jja','son')

for (seas in 1:4){
    neof = dim(eof.14dayave[[seas.name[seas]]]$eof)[2]
    mons = seas.mon[[seas]]
    seas.list = which(month(date.list)==mons[1] | month(date.list)==mons[2] | month(date.list)==mons[3])
    date.list.seas = date.list[seas.list]

    for (ne in 1:neof){
    	eof.14dayave[[seas.name[seas]]]$pc[,ne] = rm.seas.mean(eof.14dayave[[seas.name[seas]]]$pc[,ne],date.list,mons)
    }

x.list = which(day(date.list[seas.list])==1 | day(date.list[seas.list])==14)
y.list = x.list+14

if (any(x.list<0)) {
   day.rm = which(x.list<0)
   y.list = y.list[-day.rm]
   x.list = x.list[-day.rm]
}

day.rm = NULL
for (nd in 1:length(x.list)){
    if (date.list[seas.list][y.list[nd]]-date.list[seas.list][x.list[nd]]!=14) day.rm = c(day.rm,nd)
}
if (!is.null(day.rm)){
   x.list = x.list[-day.rm]
   y.list = y.list[-day.rm]
}

date.list.x = date.list.seas[x.list]
date.list.y = date.list.seas[y.list]

x.pca = eof.14dayave[[seas.name[seas]]]
x.pca$pc = x.pca$pc[x.list,]
y.pca = eof.14dayave[[seas.name[seas]]]
y.pca$pc = y.pca$pc[y.list,]
lat = x.pca$lat
lon = x.pca$lon
npc = length(eof.14dayave[[seas.name[seas]]]$pc)/neof

cca.out = cca.gev(x.pca,y.pca)
tx = cca.out$tx
ty = cca.out$ty
print(seas.name[seas])
print(paste('Correlations:',round(cca.out$can.cor,5)))
cca.crit = cca.null.rm.sig(tx,ty,npc,42,x.list,y.list,cca.cor=cca.out$can.cor)
print(paste('Crit:',round(cca.crit,5)))

cca.base.seas = list(cca.out=cca.out,cca.crit=cca.crit,date.list.x=date.list.x,date.list.y=date.list.y)
if (seas==1) {
   cca.base = list(cca.base.seas)
} else {
   cca.base[seas] = list(cca.base.seas)
}

} #end of seas in 1:4

names(cca.base) = seas.name
save(cca.base,file=paste0(dir.save,'cca.era40.week12.RData'))



