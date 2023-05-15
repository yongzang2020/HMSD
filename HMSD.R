


conditional.of.mu <- function(mu,theta,k) {

	## theta is a row of the theta.mcmc matrix (n.ite by 5)
	if (mu<=1 & mu>=0)
		return((prod(theta))^(mu*k)*(prod(1-theta))^(-mu*k)*mu^(A-1)*(1-mu)^(B-1)/beta(mu*k,(1-mu)*k)^2)
	else return(0)
}

conditional.of.k <- function(k,mu,theta) {

	## theta is a row of the theta.mcmc matrix (n.ite by 5)
	if (k>=0)
		return((prod(theta))^(mu*k)*(prod(1-theta))^((1-mu)*k)*k^(S-1)*exp(-R*k)/beta(mu*k,(1-mu)*k)^2)
	else return(0)
}

gibbs.p00p10 <- function(n,y,A,B,S,R,n.mcmc,n.burnin) {

	theta.gibs <- matrix(0,n.mcmc,2)
	mu.gibs <- rep(0,n.mcmc)
	k.gibs <- rep(0,n.mcmc)

	## set initial value
	mu <- .5
	k <- 10

	n.accept.mu <- 0
	n.accept.k <- 0
	for (i in 1:n.mcmc) {

		theta <- rep(0,2)
		for (j in 1:2) {
			theta[j] <- rbeta(1,y[j]+mu*k,n[j]-y[j]+(1-mu)*k)
			if (theta[j]<.01) theta[j] <- .01
			if (theta[j]>.99) theta[j] <- .99

		}
		theta.gibs[i,] <- theta
		sd.prop <- .3
		mu.star <- rnorm(1,mu,sd.prop)
		alpha.mu <- min(1,conditional.of.mu(mu.star,theta,k)/conditional.of.mu(mu,theta,k))
		r <- runif(1,0,1)
		if (r <= alpha.mu) {
			mu <- mu.star
			n.accept.mu <- n.accept.mu+1
		}
		mu.gibs[i] <- mu

		sd.prop <- 14
		k.star <- rnorm(1,k,sd.prop)
		alpha.k <- min(1,conditional.of.k(k.star,mu,theta)/conditional.of.k(k,mu,theta))
		r <- runif(1,0,1)
		if (r <= alpha.k) {
			k <- k.star
			n.accept.k <- n.accept.k+1
		}
		k.gibs[i] <- k


	}
	
	list(theta.gibs=theta.gibs[(n.burnin+1):n.mcmc,])
}



integrand.0 <- function(y,alpha.A,beta.A,alpha.B,beta.B,margin) {
	dbeta(y,alpha.B,beta.B)*pbeta(y-margin,alpha.A,beta.A)
}
Prob.B.greater.A <- function(alpha.A,beta.A,alpha.B,beta.B,margin) { 
	integrate(integrand.0,0,1,alpha.A,beta.A,alpha.B,beta.B,margin)$value
}



integrand <- function(y,alpha.A,beta.A,alpha.B,beta.B,delta) {
	dbeta(y,alpha.B,beta.B)*(pbeta(y+delta,alpha.A,beta.A)-pbeta(y-delta,alpha.A,beta.A))
	
}


getQ <- function(alpha.A,beta.A,alpha.B,beta.B,delta) {
	integrate(integrand,0,1,alpha.A,beta.A,alpha.B,beta.B,delta)$value
}


assign.cohort <- function(marker,current.cohort.marker.status,P0,P1,true.prob) { ## marker=0 or 1
	n.marker.current.cohort <- sum(current.cohort.marker.status==marker)

	Pm <- ifelse(marker==0,P0,P1)

	trt.assign <- rbinom(n.marker.current.cohort,1,Pm)
		
	
	
	n.0 <- sum(trt.assign==0) ## number of patients in the current cohort with the given marker that are assigned to T=0
	n.1 <- sum(trt.assign==1) ## number of patients in the current cohort with the given marker that are assigned to T=1
	y.0 <- rbinom(1,n.0,true.prob[(marker+1),1])
	y.1 <- rbinom(1,n.1,true.prob[(marker+1),2])

	list(n.0=n.0,n.1=n.1,y.0=y.0,y.1=y.1)
}

get.P0P1 <- function(n.mt,n.1mt,n.0mt,alpha.mt,beta.mt,alpha00,beta00,margin) {
		

	
	p00p10 <- gibbs.p00p10(n.mt[,1],n.1mt[,1],A,B,S,R,n.mcmc,n.burnin)$theta.gibs
	p00 <- p00p10[,1]	
	p10 <- p00p10[,2]	
	

	p01 <- rbeta(n.post, alpha.mt[1,2]+n.1mt[1,2], beta.mt[1,2]+n.0mt[1,2])
	p11 <- rbeta(n.post, alpha.mt[2,2]+n.1mt[2,2], beta.mt[2,2]+n.0mt[2,2])

	P0 <- mean(p01>p00)
	P1 <- mean(p11>p10)

	#Q <- mean(abs(p00-p10)<delta)
	#Q <- getQ(alpha.00,beta.00,alpha.10,beta.10,delta)

	list(P0=P0,P1=P1)
}


each.cohort <- function(n.mt,n.1mt,n.0mt,P0,P1,stage,coh.size,prevalence,tau.L,tau.U,true.prob,alpha.mt,beta.mt,alpha00,beta00,margin) {

	current.cohort.marker.status <- rbinom(coh.size,1,prevalence)
	recomm <- c(-10,-10)

	## action=0 terminate both marker groups. action==1, terminate at most 1 marker group
	if (P0<tau.L & P1<tau.L) {
		recomm <- c(0,0)  ## futility for both marker
		action <- 0
	} else if (P0>tau.U & P1>tau.U) {
		recomm <- c(1,1)  ## superiority for both marker
		action <- 0
	} else if (P0>tau.U & P1<tau.L) {
		recomm <- c(1,0)  ## superiority for M=0 and futility for M=1
		action <- 0
	} else if (P0<tau.L & P1>tau.U) {
		recomm <- c(0,1)  ## futility for M=0 and superiority for M=1
		action <- 0
	} else {
		action <- 1
		if (P0<tau.L | P0>tau.U) {   ## marker0 is terminated 
			recomm[1] <- ifelse(P0<tau.L,0,1)

			assign.coh <- assign.cohort(1,current.cohort.marker.status,P0,P1,true.prob) 

			n.mt[2,1] <- n.mt[2,1] + assign.coh$n.0
			n.mt[2,2] <- n.mt[2,2] + assign.coh$n.1
			n.1mt[2,1] <- n.1mt[2,1]+assign.coh$y.0
			n.1mt[2,2] <- n.1mt[2,2]+assign.coh$y.1
		
		} else if (P1<tau.L | P1>tau.U) { ## marker 1 is terminated
			recomm[2] <- ifelse(P1<tau.L,0,1)
			
			assign.coh <- assign.cohort(0,current.cohort.marker.status,P0,P1,true.prob) 

			n.mt[1,1] <- n.mt[1,1] + assign.coh$n.0
			n.mt[1,2] <- n.mt[1,2] + assign.coh$n.1
			n.1mt[1,1] <- n.1mt[1,1]+assign.coh$y.0
			n.1mt[1,2] <- n.1mt[1,2]+assign.coh$y.1


		} else { ## both marker groups are not terminated
			assign.coh.1 <- assign.cohort(1,current.cohort.marker.status,P0,P1,true.prob) 

			n.mt[2,1] <- n.mt[2,1] + assign.coh.1$n.0
			n.mt[2,2] <- n.mt[2,2] + assign.coh.1$n.1
			n.1mt[2,1] <- n.1mt[2,1]+assign.coh.1$y.0
			n.1mt[2,2] <- n.1mt[2,2]+assign.coh.1$y.1

			assign.coh.0 <- assign.cohort(0,current.cohort.marker.status,P0,P1,true.prob) 

			n.mt[1,1] <- n.mt[1,1] + assign.coh.0$n.0
			n.mt[1,2] <- n.mt[1,2] + assign.coh.0$n.1
			n.1mt[1,1] <- n.1mt[1,1]+assign.coh.0$y.0
			n.1mt[1,2] <- n.1mt[1,2]+assign.coh.0$y.1

		}
		n.0mt <- n.mt-n.1mt

		if (stage==2) {
			get.P0P1.current <- get.P0P1(n.mt,n.1mt,n.0mt,alpha.mt,beta.mt,alpha00,beta00,margin)
			if (P0>tau.L & P0<tau.U)
				P0 <- get.P0P1.current$P0
			if (P1>tau.L & P1<tau.U)
				P1 <- get.P0P1.current$P1
		}

	}
	list(recomm=recomm,n.mt=n.mt,n.1mt=n.1mt,n.0mt=n.0mt,action=action,P0=P0,P1=P1)

}

set.seed(1)

main <- function(n.sim,true.prob,alpha.mt,beta.mt,alpha00,beta00,margin,prevalence,coh.size,tau,n1,n2,tau.L,tau.U) {
	final.recomm.sim <- matrix(-100,n.sim,2)
	alloc.sim <- matrix(0,nrow=2,ncol=2)
	n.total.sim <- rep(0,n.sim)
	n.early.stop.futility.marker0 <- 0
	n.early.stop.futility.marker1 <- 0
	n.early.stop.superiority.marker0 <- 0
	n.early.stop.superiority.marker1 <- 0
	overall.response.rate <- rep(0,n.sim)
	for (ite in 1:n.sim) {

		### stage I

		marker.status <- rbinom((n1*coh.size),1,prevalence)
		trt.ind <- rbinom((n1*coh.size),1,0.5)  ## equal randomization
		n.mt <- matrix(0,nrow=2,ncol=2) ## number of patients in each of the 4 arms
		n.1mt <- n.mt ## number of patients with Y==1 in each of the 4 arms
		for (i in 1:2)
			for (j in 1:2) {
				n.mt[i,j] <- sum(marker.status==(i-1) & trt.ind==(j-1))
				n.1mt[i,j] <- rbinom(1,n.mt[i,j],true.prob[i,j])
			}
		n.0mt <- n.mt-n.1mt

		get.P0P1.current <- get.P0P1(n.mt,n.1mt,n.0mt,alpha.mt,beta.mt,alpha00,beta00,margin)
		P0 <- get.P0P1.current$P0
		P1 <- get.P0P1.current$P1

		n <- n1*coh.size

	
		### stage II
		for (coh.ind in 1:n2) {
			coh.current <- each.cohort(n.mt,n.1mt,n.0mt,P0,P1,stage=2,coh.size,prevalence,tau.L,tau.U,true.prob,alpha.mt,beta.mt,alpha00,beta00,margin)
			action <- coh.current$action
			final.recomm <- coh.current$recomm

			n.mt <- coh.current$n.mt
			n.1mt <- coh.current$n.1mt
			n.0mt <- coh.current$n.0mt

			if (action==0) {
				final.recomm.sim[ite,] <- final.recomm
				if (coh.ind<n2) {
					n.early.stop.futility.marker0 <- n.early.stop.futility.marker0 + ifelse(final.recomm[1]==0,1,0)
					n.early.stop.futility.marker1 <- n.early.stop.futility.marker1 + ifelse(final.recomm[2]==0,1,0)
					n.early.stop.superiority.marker0 <- n.early.stop.superiority.marker0 + ifelse(final.recomm[1]==1,1,0)
					n.early.stop.superiority.marker1 <- n.early.stop.superiority.marker1 + ifelse(final.recomm[2]==1,1,0)
				}
				break
			} else {
				P0 <- coh.current$P0
				P1 <- coh.current$P1
				n <- n+coh.size
				if (coh.ind==(n2-1)) {
					n.early.stop.futility.marker0 <- n.early.stop.futility.marker0 + ifelse(final.recomm[1]==0,1,0)
					n.early.stop.futility.marker1 <- n.early.stop.futility.marker1 + ifelse(final.recomm[2]==0,1,0)
					n.early.stop.superiority.marker0 <- n.early.stop.superiority.marker0 + ifelse(final.recomm[1]==1,1,0)
					n.early.stop.superiority.marker1 <- n.early.stop.superiority.marker1 + ifelse(final.recomm[2]==1,1,0)

				}
			}
		}
		
		
		if (action!=0) {
			if (final.recomm[1] == -10) final.recomm[1] <- ifelse(P0>tau,1,0)
			if (final.recomm[2] == -10) final.recomm[2] <- ifelse(P1>tau,1,0)
			final.recomm.sim[ite,] <- final.recomm
		}
		
		alloc.sim <- alloc.sim + n.mt
		n.total.sim[ite] <- n
	
		overall.response.rate[ite] <- sum(n.1mt)/sum(n.mt)
	}

	list(recomm=final.recomm.sim,final.recomm=apply(final.recomm.sim,2,mean),alloc=alloc.sim/n.sim,alloc.mean=sum(alloc.sim/n.sim),n=mean(n.total.sim),earlyStop.futility.marker0=n.early.stop.futility.marker0/n.sim,
	earlyStop.futility.marker1=n.early.stop.futility.marker1/n.sim,earlyStop.sup.marker0=n.early.stop.superiority.marker0/n.sim,earlyStop.sup.marker1=n.early.stop.superiority.marker1/n.sim,response=mean(overall.response.rate))
}


### truth
true.prob <- matrix(c(.3,.1,.3,.15),nrow=2,byrow=T)
### p00 p01
### p10 p11

### prior
alpha.mt <- matrix(c(.3,.3,.3,.3),nrow=2,byrow=T)
beta.mt <- matrix(c(.7,.7,.7,.7),nrow=2,byrow=T)
alpha00 <- .3;beta00 <- .7


margin <- .05
prevalence <- .5  ## prevalence of marker positive
coh.size <- 5
#delta <- .05
#eps <- .05
tau.L <- 0.1
tau.U <- .9
tau <- .8
#delta1 <- .25 # cutoff to assign all patients to control
#delta2 <- .7 # cutoff to assign all patients to experimental arm

n1 <- 4 # number of cohorts in stage I
n2 <- 20 # number of cohorts in stage II

## parameterize the beta distribution in terms of mean mu and sample size k. mu~beta(A,B), k~gamma(S,R)
A <- 2
B <- 2
S <- 1   ## shape parameter of k
R <- .1  ## rate parameter of k

n.mcmc <- 5200
n.burnin <- 200
n.post <- 5000

set.seed(1)
main(200,true.prob,alpha.mt,beta.mt,alpha00,beta00,margin,prevalence,coh.size,tau,n1,n2,tau.L,tau.U)












