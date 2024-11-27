rm.seas.mean = function(pc,date.list,mons){

### MONS = C(12,1,3) FOR EXAMPLE
### DATE.LIST NEEDS TO BE THE FULL SET OF DATES SPANNING THE DATASET

library(lubridate)

seas.list = NULL
for (nm in 1:length(mons)){ 
    seas.list = c(seas.list,which(month(date.list)==mons[nm]))
}
seas.list = sort(seas.list)

if (length(seas.list)!=length(pc)) stop('In rm.seas.mean, pc length not consistent with season list - the field entered should already be only that season.')

seas.breaks = NULL
for (i in 1:(length(seas.list)-1)) {
    seas.breaks[i] = seas.list[i+1]-seas.list[i]
}
loc.breaks = which(seas.breaks>90)

### REMOVE LOCAL SEASONAL MEAN
prev.loc.breaks = 0
for (ny in 1:length(loc.breaks)){
    pc[(prev.loc.breaks+1):loc.breaks[ny]] = pc[(prev.loc.breaks+1):loc.breaks[ny]]-mean(pc[(prev.loc.breaks+1):loc.breaks[ny]])
    prev.loc.breaks = loc.breaks[ny]
}

# REMOVE SEASONAL MEAN OF LAST SEASON
npc = length(pc)
pc[(loc.breaks[ny]+1):npc] = pc[(loc.breaks[ny]+1):npc]-mean(pc[(loc.breaks[ny]+1):npc])

pc

}