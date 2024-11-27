rm.seas.cycle = function(field,nharm=3,npoly=3,include.poly=TRUE){


t = seq(1,length(field))
t.list = array(NA,dim=c(length(t),2*nharm+npoly))
for (j in seq(1,2*nharm,2)) {
    t.list[,j] = cos((j+1)*pi*t/365.25)
    t.list[,j+1] = sin((j+1)*pi*t/365.25)
}
if (include.poly){
   poly.curr = 1
   for (j in seq(2*nharm+1,2*nharm+npoly)){
    	t.list[,j] = t^poly.curr
	poly.curr = poly.curr + 1
    }}
    
anoms = residuals(lm(field~t.list,na.action=na.exclude))

anoms

}