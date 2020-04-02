## Code to reproduce Figures 2 & 3 of Tabi et al 2020. Species multidimensional effects explain idiosyncratic responses of communities to environmental change

# Load  packages --------------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
require(mvtnorm)
require(kernlab)
library(data.table)
library(grid)
library(scales)
library(matrixcalc)
library(cowplot)
library(grid)
library(gridExtra)

# Functions -------------------------------------------------------------------------------------------------------------------------------------------

# generate interaction matrix (A) where you can set the number of species, the level of variance and the connectance
generate_Interaction_matrix <- function(S, P, conne){
  Inte <- runif(S*S,  min = -1, max = 1 )*P 
  zeroes <- sample(c(rep.int(1,floor(S*S*conne)), rep.int(0,(S*S-floor(S*S*conne)))))
  Inte[which(zeroes==0)] <- 0
  Inte <- matrix(Inte, ncol = S, nrow = S)
  diag(Inte) <- -1
  Inte <- -Inte
  return(Inte)
}

# to calculate the omega (the size of the feasibility domain)
Omega <- function(A){
  n <- nrow(A)
  Sigma <-solve(t(A) %*% A)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- d[1]
  return(out) 
}

# to get the span of the column vectors of the A matrix
spanned_vectors <- function(A, num){
  G <- matrix(0, ncol=ncol(A), nrow=nrow(A))
  for(k in 1:num) G[,k] <- -A[,k]
  G
}


# to get the vector of carrying capacities that yield equal equilibrium abundances
parameterization_center <- function(A, num){
  G <- spanned_vectors(A, num)
  lambda <- rep.int(1, num)
  lambda <- lambda/sum(lambda)
  growth <- G %*% matrix(lambda, ncol=1) %>%
    as.vector()
}

# Simulations to reproduce Figure 2 -------------------------------------------------------------------------------------------------------------------

num <- 500    # number of communities
sp <- c(2,3)  # number of species
p <- c(0.1, 0.5, 0.9)  # level of variance (asymmetry)
C <- 1        # connectance

sim_res <- array( dim= c(num, 4, length(p), length(sp)), dimnames=list( 1:num, c("evenness", "theta", "omega", "phi"), p, sp  ) )


for(u in 1:length(sp)) {
  for(j in 1:length(p) ) {
    for(i in 1:num) {
      
      P <- p[j]
      S <- sp[u]
      
      # generate interaction matrix
      alpha <- generate_Interaction_matrix(S, P , C) 
      
      # calculate omega (size of the feasibility domain)
      sim_res[ i, "omega", j, u] <- tryCatch( Omega(alpha) ,error = function(e) {NA})
      
      # calculate phi (asymmetry of the feasibility domain)
      sv <- spanned_vectors(alpha,S) ### 
      sim_res[i, "phi", j, u] <- sd(sqrt(colSums(sv^2))) ###### 
      
      N <- exp(rnorm(S)*5*runif(1)) # equilibrium biomass
      K <- alpha %*% N  # calculate carrying capacities 
      Kc <- -parameterization_center(alpha,S)   # to calculate the centroid
      
      # calculate evenness
      pr <- N /sum(N)  
      sim_res[i, "evenness", j, u] <- -sum(pr*log10(pr))/log10(S)  
      
      # calculate theta (relative performance in isolation)
      sim_res[i, "theta", j, u] <- acos(sum(K*Kc)/(sqrt(sum(K^2))*sqrt(sum(Kc^2))))
      
    } 
  }
}



dd <- spread(reshape2::melt(sim_res ), Var2, value)

names(dd)[1:3] <-  c( "num", "P", "S")

# relative performance in isolation (theta) normalized by the size of the feasibility domain (omega)
dd$norm_theta <- dd$theta * (0.5-dd$omega)

dd$sp <- factor(dd$S)
levels(dd$sp) <- c("S = 2", "S = 3")
dd$p <- factor(dd$P)
levels(dd$p) <- c("low", "medium", "high")

# Spearman's rank correlation between evenness and theta(omega)
cors <- dd %>% group_by(  sp, p) %>% dplyr::summarise(cor = cor.test( norm_theta,evenness,method = 'spearman')$est )


fancy_scientific <- function(l) {
  
  ifelse( 0 < l & l < 0.01,  formatC(l, format = "e",  digits = 0), round(l,2) )
  
}

coord_x <- dd %>% group_by(sp, p) %>% dplyr::summarise(coord = range( norm_theta, na.rm=T)[2]*0.7 )

Fig2 <- ggplot()+
  geom_point(data=dd, aes(x=norm_theta, y=evenness, fill=factor(p), color=factor(p) ), size=2, pch=21, alpha=0.6)+
  facet_wrap(sp ~ p, scales = "free_x") +
  scale_x_continuous(labels = fancy_scientific ) +
  geom_text(data=cors, aes(label=paste( "rho == ", round(cor,2),   sep="")), x=rep(coord_x$coord,1), y=0.95, size=6, hjust=0, parse = TRUE,
            color="blue", fontface="italic" )+
  ylab("Evenness")+
  xlab(  expression( paste( "Relative Performance in Isolation (", theta[n], ")"  )  ))+
  scale_fill_manual(values=c("cyan3", "gold", "brown1") ) +
  scale_color_manual(values=c("darkblue", "orange", "darkred") ) +
  ylim(0,1) +
  theme_bw()+ theme(legend.direction = 'horizontal', legend.position = 'none') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1) )+
  theme(strip.background =element_rect(colour="white", fill="white") )+
  theme(strip.text = element_text(size = 16), axis.text = element_text(size=14) )+
  theme(plot.title = element_text(lineheight=.8, face="bold", size=20, hjust = 0.5), axis.title = element_text(size=20))


Fig2


# Figure 3 --------------------------------------------------------------------------------------------------------------------------------------


labs <- dd %>% group_by(sp,p) %>% dplyr::summarise(dist = round( quantile( omega, c(0.75), na.rm=T ) - quantile( omega, c(0.25), na.rm=T ), 2))

meds <- dd %>% group_by(sp,p) %>% dplyr::summarise( med = round( median( omega,  na.rm=T ), 3), 
                                                    u = round( quantile( omega, c(0.75), na.rm=T ), 3),
                                                    l = round( quantile( omega, c(0.25), na.rm=T ), 3) )

meds

dodge=position_dodge(width=0.5) 

omegas <- ggplot(data=meds, aes(x=p, y=med, color=p, shape=sp)) +   
  geom_point(  size=5, stroke=2, position=dodge )+
  geom_errorbar(data=meds, aes(x=p, y=med, ymax=u, ymin=l, color=p), width=0.1, size=1, position=dodge)+
  geom_text(data=meds, aes(label =labs$dist, y = u, x=p ), vjust = -.5, color="black", position=dodge, size=6) +
  scale_color_manual(values=c("cyan3", "gold", "brown1") ) +
  scale_shape_manual(values=c(1,2,5) ) +
  ylab(  expression( paste("Size of the Feasibility Domain"  ))  )+ xlab("")+
  theme_bw()+ theme(legend.direction = 'horizontal', legend.position = 'none') +
  ylim(0,0.35)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1) )+
  theme( axis.text = element_text(size=18, colour = "black"), axis.title = element_text(size=20))

omegas


labs2 <- dd %>% dplyr::group_by(sp,p) %>% dplyr::summarise(dist = round( quantile( theta, c(0.75), na.rm=T ) - quantile( theta, c(0.25), na.rm=T ), 2))

meds2 <- dd %>% group_by(sp,p) %>% dplyr::summarise( med = round( median( theta,  na.rm=T ), 3),
                                                     u = round( quantile( theta, c(0.75), na.rm=T ), 3), 
                                                     l = round( quantile( theta, c(0.25), na.rm=T ), 3) )
meds2

dodge=position_dodge(width=0.5) 

thetas <- ggplot(data=meds2, aes(x=p, y=med, color=p, shape=sp)) +   
  geom_point(  size=5, stroke=2, position=dodge )+
  geom_errorbar(data=meds2, aes(x=p, y=med, ymax=u, ymin=l, color=p), width=0.1, size=1, position=dodge)+
  geom_text(data=meds2, aes(label =labs2$dist, y = u, x=p ), vjust = -.5, color="black", position=dodge, size=6) +
  scale_shape_manual(values=c(1,2,5) ) +
  scale_color_manual(values=c("cyan3", "gold", "brown1") ) +
  ylab(  expression( paste( "Relative Performance in Isolation" ))  )+ xlab("")+
  theme_bw()+ theme(legend.direction = 'horizontal', legend.position = 'none') +
  ylim(0,1)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1) )+
  theme( axis.text = element_text(size=18, colour = "black"), axis.title = element_text(size=20))

thetas



labs3 <- dd %>% group_by(sp,p) %>% dplyr::summarise(dist = round( quantile( phi, c(0.75), na.rm=T ) - quantile( phi, c(0.25), na.rm=T ), 2))

meds3 <- dd %>% group_by(sp,p) %>% dplyr::summarise( med = round( median( phi,  na.rm=T ), 3),
                                                     u = round( quantile( phi, c(0.75), na.rm=T ), 3),
                                                     l = round( quantile( phi, c(0.25), na.rm=T ), 3) )

meds3

dodge=position_dodge(width=0.5) 

phis <- ggplot(data=meds3, aes(x=p, y=med, color=p, shape=sp)) +   
  geom_point(  size=5, stroke=2, position=dodge )+
  geom_errorbar(data=meds3, aes(x=p, y=med, ymax=u, ymin=l, color=p), width=0.1, size=1, position=dodge)+
  geom_text(data=meds3, aes(label =labs3$dist, y = u, x=p ), vjust = -.5, color="black", position=dodge, size=6) +
  scale_shape_manual(values=c(1,2,5) ) +
  scale_color_manual(values=c("cyan3", "gold", "brown1") ) +
  ylab(  expression( paste( "Asymmetry of the Feasibility Domain" ))  )+
  xlab("")+
  theme_bw()+ theme(legend.direction = 'horizontal', legend.position = 'none') +
  ylim(0,0.2)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1) )+
  theme( axis.text = element_text(size=18, colour = "black"), axis.title = element_text(size=20))

phis


head(dd)
dd$fsp <- as.factor(dd$sp)
levels(dd$fsp) <- c("2-species", "3-species")

L <- ggplot(dd) + 
  geom_point(aes(x=omega, y=theta, shape=fsp ), size=4, stroke=1, position=dodge )+
  theme_bw()+ 
  scale_shape_manual(values=c(1,2,5) ) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size=1) )+
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=20), legend.direction = "horizontal",
        legend.box.background = element_rect(colour = "black"))

legend <- get_legend(L)

grid <- plot_grid(  omegas, thetas , phis,
                    labels = c(  "A", "B", "C"), 
                    ncol=3, align="h",  label_size = 22  )

p_grid <- plot_grid( legend,grid, ncol = 1, rel_heights = c(.2, 1) )

x.grob <- textGrob( label=expression(paste("Level of Asymmetry in the Feasibility Domain"  )), 
                    gp=gpar(col="black", fontsize=20))

Fig3 <- grid.arrange(arrangeGrob(p_grid, bottom = x.grob))