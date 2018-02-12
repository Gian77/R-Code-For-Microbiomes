######################################################
##  *** plot a regression equation from lm() ***    ##
##        on a ggplot regression-line plot          ##
##				                    ##
## which I have modified from:                      ##
## http://stackoverflow.com/q/7549694/857416        ##
##                                                  ##
##              Gian MN Benucci, PhD                ##
##          Michigan State University               ##
##                benucci@msu.edu                   ##
##              February 2nd, 2018                  ##
######################################################

lm_eqn = function(m) {
  # Displays regression line equation and R^2 value on plot
  # Usage:
  # p + annotate("text", x=25, y=300, label=lm_eqn(lm(y ~ x, df)), parse=TRUE)
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(signif(summary(m)$r.squared), digits = 2, nsmall=2),
  				p = format(summary(m)$coef[2,4], digits = 2));
  
	if ((coef(m)[2] >= 0) & (summary(m)$coef[2,4]<=0.01)) {
		eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)<="0.01",l)
  		} else if ((coef(m)[2] >= 0) & (summary(m)$coef[2,4]>0.01)) {
		eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)~"="~p,l)
		} else if ((coef(m)[2] < 0) & (summary(m)$coef[2,4]<=0.01)) {
    	eq <- substitute(italic(y) == a - b * italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)<="0.01",l)
		} else if ((coef(m)[2] < 0) & (summary(m)$coef[2,4]>0.01)) {
		eq <- substitute(italic(y) == a - b * italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)~"="~p,l)  
	}
  
  as.character(as.expression(eq));                 
}



