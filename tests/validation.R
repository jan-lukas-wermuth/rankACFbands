library(rankACFbands)
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label)
  checks <<- checks + 1L
}
near <- function(x, y, label, tolerance = 1e-9) {
  check(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance)), label)
}
fails <- function(expr, label) check(inherits(tryCatch({force(expr); NULL}, error = identity), "error"), label)
mid <- function(q, z) as.numeric(q > z) + 0.5 * as.numeric(q == z)
naive <- function(x, H) {
  T <- length(x)
  G <- function(q) mean(mid(q, x))
  p <- function(q) mean(x == q)
  rho <- tau <- nu <- numeric(H + 1)
  pr <- list()
  for (h in 0:H) {
    u <- x[seq_len(T-h)]; v <- x[seq_len(T-h)+h]; n <- length(u)
    pairs <- combn(n, 2)
    signs <- sign(u[pairs[1,]]-u[pairs[2,]]) * sign(v[pairs[1,]]-v[pairs[2,]])
    tau[h+1] <- mean(signs); nu[h+1] <- mean(signs == 0)
    rho[h+1] <- 12/T * sum((vapply(u,G,0)-0.5)*(vapply(v,G,0)-0.5))
    GH <- function(a,b) mean(mid(a,u)*mid(b,v))
    if (h == 0) {
      pr[["rho.0"]] <- 1 - vapply(x,p,0)^2 - rho[1]
      pr[["tau.0"]] <- 1 - vapply(x,p,0) - tau[1]
      pr[["nu.0"]] <- vapply(x,p,0) - nu[1]
    } else {
      pr[[paste0("tau.",h)]] <- vapply(seq_len(n), function(i) 4*GH(u[i],v[i])-2*(G(u[i])+G(v[i]))+1-tau[h+1],0)
      pr[[paste0("nu.",h)]] <- vapply(seq_len(n), function(i) p(u[i])+p(v[i])-mean(u==u[i]&v==v[i])-nu[h+1],0)
      pr[[paste0("rho.",h)]] <- vapply(seq_len(n),function(i) {
        g1 <- mean(vapply(x,function(y) GH(u[i],y),0))
        g2 <- mean(vapply(x,function(y) GH(y,v[i]),0))
        4*(g1+g2+G(u[i])*G(v[i])-G(u[i])-G(v[i]))+1-rho[h+1]
      },0)
    }
  }
  list(rho=rho,tau=tau,nu=nu,projections=pr)
}
set.seed(902)
series <- list(c(1,1,3,2,4,2,1,3), rnorm(19), sample(0:4,25,TRUE),
               c(-3,0,0,0.2,8,0.2,-3,2,2,1), seq_len(12))
for (x in series) {
  H <- 3; ref <- naive(x,H)
  expected <- list(rho=ref$rho,tau=ref$tau,rho_b=ref$rho/ref$rho[1],
                   tau_b=ref$tau/ref$tau[1],gamma=ref$tau/(1-ref$nu))
  for (m in names(expected)) {
    fit <- rank_acf(x,H,m)
    near(fit$acf,expected[[m]],paste("literal estimator",m))
    near(rank_acf(-x,H,m)$acf,fit$acf,paste("decreasing transform",m))
    near(rank_acf(exp(x/10),H,m)$acf,fit$acf,paste("increasing transform",m))
  }
  pr <- rank_projections(x,H)
  for (nm in names(pr)) near(pr[[nm]],ref$projections[[nm]],paste("literal projection",nm))
  near(rank_acf(x,H,"rho_b")$acf,as.numeric(stats::acf(rank(x),lag.max=H,plot=FALSE)$acf),"grade ACF matches Pearson ACF of full ranks")
}
# Exact combinatorial check of tau-a and gamma, distinct from cor()'s tau-b.
x <- c(1,1,2,3,4,4,2,5)
for(h in 1:3) {
  u<-head(x,-h); v<-tail(x,-h); n<-length(u)
  tx<-1-sum(table(u)*(table(u)-1))/(n*(n-1))
  ty<-1-sum(table(v)*(table(v)-1))/(n*(n-1))
  near(rank_acf(x,3,"tau")$acf[h+1],cor(u,v,method="kendall")*sqrt(tx*ty),"Kendall tau-a vs R tau-b")
}
# Independently evaluate every HAC entry with literal double loops.
x <- c(1,2,1,3,2,4,1,2,4,3,2,1,3,4,2,3,1,4)
for(m in c("rho","tau","rho_b","tau_b","gamma")) {
  v <- rank_hac(x,3,m,bandwidth=2,psd="none")
  ref <- naive(x,3); labels<-colnames(v$jacobian)
  expected<-matrix(0,length(labels),length(labels))
  for(a in seq_along(labels)) for(b in seq_along(labels)) {
    ha<-as.integer(sub(".*\\.","",labels[a])); hb<-as.integer(sub(".*\\.","",labels[b]))
    ra<-if(startsWith(labels[a],"rho"))3 else 2; rb<-if(startsWith(labels[b],"rho"))3 else 2
    for(s in -2:2) for(t in seq_len(length(x)-abs(s)-max(ha,hb))) {
      aa<-ref$projections[[labels[a]]]; bb<-ref$projections[[labels[b]]]
      expected[a,b]<-expected[a,b]+ra*rb/length(x)*(1-abs(s)/3)*
        if(s>=0) aa[t]*bb[t+s] else bb[t]*aa[t-s]
    }
  }
  near(v$primitive_covariance,expected,paste("literal HAC",m))
  near(v$raw_covariance,v$jacobian%*%expected%*%t(v$jacobian),paste("delta transform",m))
  near(vcov(v),v$covariance/length(x),"vcov scaling")
}
# Numerical derivatives check Jacobian signs and ordering.
for(m in c("rho_b","tau_b","gamma")) {
  v<-rank_hac(x,3,m,bandwidth=0,psd="none")
  raw<-rank_acf(x,3,m)$raw
  f<-if(m=="gamma") function(z) z[c(1,3,5)]/(1-z[c(2,4,6)]) else function(z) z[-1]/z[1]
  z<-if(m=="gamma") as.vector(rbind(raw$tau[-1],raw$nu[-1])) else raw[[if(m=="rho_b")"rho" else "tau"]]
  J<-sapply(seq_along(z),function(i){plus<-minus<-z;plus[i]<-plus[i]+1e-6;minus[i]<-minus[i]-1e-6;(f(plus)-f(minus))/2e-6})
  near(v$jacobian,J,paste("numerical Jacobian",m),1e-7)
}
# IID variances, Box-Pierce scaling, exact independent critical values.
ties<-rank_ties(x)
for(m in c("rho","tau","rho_b","tau_b","gamma")) {
  b<-rank_acf_sigbands(x,3,m,plot=FALSE)
  near(b$critical,qnorm((1+0.95^(1/3))/2),"exact Sidak")
  near(b$table$se,rep(sqrt(rank_iid_variance(x,m)/length(x)),3),"iid SE")
  q<-rank_box_test(x,3,m)
  near(unname(q$statistic),length(x)*sum(coef(rank_acf(x,3,m))^2)/rank_iid_variance(x,m),"Box-Pierce uses variance")
}
near(rank_iid_variance(x,"rho_b"),1,"grade iid variance one")
near(supt_critical(diag(1)),qnorm(.975),"H=1")
near(supt_critical(diag(3),alternative="greater"),qnorm(.95^(1/3)),"one-sided iid critical")
# Correlated integration, RNG preservation, singular PSD, one-sided semantics.
set.seed(123); before<-.Random.seed
v<-matrix(c(1,.5,.5,1),2)
a<-supt_critical(v,maxpts=1000000,seed=9)
check(identical(before,.Random.seed),"local RNG restoration")
near(a,supt_critical(v,maxpts=1000000,seed=9),"seed reproducibility")
check(supt_critical(v,alternative="greater",maxpts=1000000,seed=9)<=a,"one-sided threshold")
near(supt_critical(matrix(1,2,2),maxpts=1000000,seed=2),qnorm(.975),"singular perfect correlation",.04)
set.seed(33); x<-as.numeric(arima.sim(list(ar=.3),n=180))
for(m in c("rho","tau","rho_b","tau_b","gamma")) {
  b<-rank_acf_confbands(x,3,m,alternative="greater",maxpts=1000000,plot=FALSE)
  check(all(is.infinite(b$table$upper)&b$table$upper>0),"lower confidence band")
  near(b$table$lower,b$table$estimate-b$critical*b$table$se,"lower confidence formula")
  b2<-rank_acf_confbands(x,3,m,alternative="less",maxpts=1000000,plot=FALSE)
  check(all(is.infinite(b2$table$lower)&b2$table$lower<0),"upper confidence band")
  near(b$critical,b2$critical,"Gaussian sign symmetry")
  s<-rank_acf_sigbands(x,3,m,alternative="greater",plot=FALSE)
  check(all(is.infinite(s$table$lower)),"upper significance threshold")
  near(s$table$upper,s$critical*s$table$se,"significance threshold formula")
}
pt<-rank_acf_sigbands(x,3,type="pointwise",plot=FALSE)
check(is.na(pt$p.value)&&is.na(pt$reject),"no false global pointwise decision")
# Inputs and degeneracy must fail clearly.
fails(rank_acf(c(1,NA,3,4)),"missing values")
fails(rank_acf(c(1,Inf,3,4)),"infinite values")
fails(rank_acf(rep(2,10)),"constant normalized series")
fails(rank_acf(1:10,lag.max=9),"insufficient pairs")
fails(rank_acf(1:10,lag.max=2.5),"fractional lag")
fails(rank_acf(factor(c("a","b","a"))),"unordered factor")
fails(supt_critical(diag(c(1,0))),"zero variance")
fails(supt_critical(matrix(c(1,2,2,1),2)),"indefinite covariance")
fails(rank_box_test(x,2,null=.2),"iid nonzero null")
near(rank_acf(ordered(c("low","high","low","mid","high"),levels=c("low","mid","high")),2)$acf,
     rank_acf(c(1,3,1,2,3),2)$acf,"ordinal input")
# Plot and S3 smoke tests, including one lag and one-sided bounds.
pdf(file=tempfile(fileext=".pdf"))
plot(rank_acf(x,1)); plot(rank_acf_sigbands(x,1,alternative="less",plot=FALSE))
plot(rank_acf_confbands(x,2,maxpts=1000000,plot=FALSE))
dev.off()
check(nrow(as.data.frame(rank_acf(x,2)))==3,"ACF data frame")
cat("Validated",checks,"mathematical, input, inference and interface checks.\n")
