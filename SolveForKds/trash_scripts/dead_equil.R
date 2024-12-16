# ModelFunction <- function(pars) {
#   # Split up the parameters into the kd and background parameters.
#   kds  <- 10^pars[1 : nrow(sXc)]
#   names(kds) <- rownames(sXc)
#   bgs  <- 10^pars[length(kds) + 1]
#   stock.ago <- 10^pars[length(kds) + 2]
#   c.agos <- dil.agos * stock.ago
#   # c.agos.old <<- c.agos
#   c.bounds <- as.matrix(
#     sapply(c.agos, function(c.ago) {
#       return(BoundRNA(kds, c.I.tots, c.ago))
#     }
#   ))
#   # print("c.bounds.old")
#   # print(c.bounds)
#   # c.bounds.old <<- c.bounds
#   # Get the amount of background binding by subtracting the bound from the
#   # total sites in each experiment, normalizing. Must transpose to multiply
#   # each column.
#   c.frees <- c.totals - c.bounds
#   c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
#   c.all <- c.bounds + c.bgs
#   # print("c.all.old")
#   # print(c.all)
#   # c.all.old <<- c.all
#   c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
#   # c.final.old <<- c.final
#   return(c.final)
# }

# ModelLikelihood <- function(pars){
#   model <- ModelFunction(pars)
#     # lc.final <- log(model+1)
#   # ldata <- log(data+1)
#   model_norm <- t(t(model) / colSums(model))
#   sumofsquares <- -sum(data*log(model_norm))
#   return(sumofsquares)
# }

# Gradient <- function(pars) {
#   # Split up the parameters into the kd and background parameters.
#   kds  <- 10^pars[1 : nrow(sXc)]
#   names(kds) <- rownames(sXc)
#   B  <- 10^pars[length(kds) + 1]
#   stock.ago <- 10^pars[length(kds) + 2]
#   c.agos <- stock.ago * dil.agos
#    f.mat <- matrix(
#       rep(sapply(c.agos, function(c.ago) {
#            return(FreeAgo(kds, c.I.tots, c.ago))
#          }
#          ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)

#   f.jvec <- sapply(c.agos, function(c.ago) {
#            return(FreeAgo(kds, c.I.tots, c.ago))})

#   A.mat <- matrix(c.agos, nrow=nrow(data), ncol=ncol(data),
#                   byrow=TRUE)

#   A.jvec <- c.agos

#   K.mat <- matrix(kds, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
#   K.ivec <- kds

#   l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
#   l.ivec <- c.I.tots

#   R.mat <- matrix(colSums(data), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
#   R.jvec <- colSums(data) 
#   L <- 100

#   # Get the amount of background binding by subtracting the bound from the
#   # total sites in each experiment, normalizing. Must transpose to multiply
#   # each column.

#   time.init <- proc.time()
#   c.vec_num <- -(
#                  (l.ivec %*% t(R.jvec * f.jvec * (f.jvec + L - A.jvec)) +
#                  (l.ivec * K.ivec * B)  %*% t(R.jvec))
#                 ) 
#   C1.jvec <- L - B - 2 * A.jvec
#   C2.jvec <- A.jvec^2 + (B - L) * A.jvec - L * B


#   c.vec_dem <- t(
#                 (f.jvec^3 + f.jvec^2 * C1.jvec + f.jvec * C2.jvec) +
#                 t(K.ivec %*% t(f.jvec^2 + f.jvec * C1.jvec + C2.jvec))
#                 )^-1
#   c.final <- c.vec_num * c.vec_dem

#   time.init <- proc.time()
#   dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)
#   time.new <- proc.time()

#   dF.dK.mat <- sweep(f.mat * l.mat * (f.mat + K.mat)^(-2),MARGIN=2,dF.base, "*")

#   C1.mat <- L - B - 2*A.mat
#   C2.mat <- A.mat^2 + (B - L)*A.mat - L*B

#   # The d (each model point) d Free derivative:)
#   dc.ai.dF <- R.mat * l.mat * (
#     (
#       f.mat^4
#     ) + (
#       2 * (L - A.mat) * f.mat^3
#     ) + (
#       ((L - A.mat)^2 + K.mat * (4 * B + A.mat)) * f.mat^2
#     ) + (
#       2 * K.mat * (A.mat * (L - B) + B * (K.mat - 2 * L) - (A.mat + B)^2) * f.mat
#     ) + (
#       K.mat * (A.mat^3 + 2 * A.mat^2 * (B - L) - B * (B - L) * (K.mat + L) +
#                A.mat * (B^2 - 2 * B * K.mat - 3 * B * L + L^2))
#     )
#   ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

#   time.init <- proc.time()


#   dc.ai.dKj_specific <- R.mat * l.mat * (
#     (
#       f.mat^4
#     ) + (
#       (2 * C1.mat + A.mat) * f.mat^3
#     ) + (
#       (C2.mat + (C1.mat + A.mat) * C1.mat) * f.mat^2
#     ) + (
#       (C2.mat * (C1.mat + A.mat)) * f.mat
#     )
#   ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

#   dc.ai.dA_specific <- R.mat * l.mat * (
#     (
#       -f.mat^3
#     ) + (
#       (2 * (A.mat - L)) * f.mat^2
#     ) + (
#       ((A.mat - L) * C1.mat - 2 * B * K.mat + C2.mat) * f.mat
#     ) + (
#       - C1.mat * B * K.mat
#     )
#   ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

#   dc.ai.dB <- R.mat * l.mat * (
#     (
#       -f.mat^3
#     ) + (
#       (2 * (A.mat - L) - K.mat) * f.mat^2
#     ) + (
#       ((L - A.mat) * (A.mat - L) - K.mat * (B + C1.mat)) * f.mat
#     ) + (
#       + A.mat * B * K.mat - B * K.mat * L - K.mat * C2.mat
#     )
#   ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

#   residuals <- -data/c.final

#   grad_derivs <- (log(10)*kds *
#                   (colSums(residuals*dc.ai.dF) %*%
#                    t(dF.dK.mat
#                   ) + rowSums(dc.ai.dKj_specific*residuals)))
  
#   dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
#                     MARGIN=2,dil.agos, "*")

#   colnames(c.final) <- colnames(data)

#   gradient_with_Ago <- log(10)*stock.ago*sum(residuals * dc.ai.dA)

#   gradient_with_bg <- log(10)*B*sum(residuals * dc.ai.dB)
#   gradient_all <- c(grad_derivs, gradient_with_bg, gradient_with_Ago)
#   names(gradient_all) <- names(pars)
#   return(gradient_all)
# }
