#check for required packages

list.of.packages <- c("plspm", "turner", "amap")
new.packages <- list.of.packages[
			!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


#modified from plspm/R/res.clus.r
#takes in a plspm object as input

resid <-  function (pls, Y = NULL) {

    if (class(pls) != "plspm") 
        stop("\n'res.clus()' requires a 'plspm' object")

    if (any(pls$model$specs$modes != "A")) 
        stop("\nSorry, REBUS only works for mode 'A'")

    if (!pls$model$specs$scaled) 
        stop("\nSorry, REBUS only works with scaled='TRUE'")

    test_dataset(Y, pls$data, pls$model$gens$obs)
    IDM <- pls$model$IDM
    blocks <- pls$model$blocks
    blocklist = indexify(blocks)

    if (!is.null(pls$data)) {
        DM = pls$data
        dataset = TRUE
    } else {
        dataset = FALSE
        DM = get_manifests(Y, blocks)
    }

    lvs = nrow(IDM)
    lvs.names = rownames(IDM)
    mvs = pls$model$gen$mvs
    X = get_data_scaled(DM, TRUE)
    Y.lvs <- pls$scores
    loads <- pls$outer_model$loading
    Path <- pls$path_coefs
    endo <- rowSums(IDM)
    endo[endo != 0] <- 1
    outer_residuals = DM
    inner_residuals = Y.lvs[, endo == 1]

    for (j in 1:lvs) {
        X.hat = Y.lvs[, j] %*% t(loads[blocklist == j])
        outer_residuals[, blocklist == j] = X[, blocklist == 
            j] - X.hat
    }

    if (sum(endo) != 1) 
        Y.hat <- Y.lvs %*% t(Path[endo == 1, ])

    if (sum(endo) == 1) 
        Y.hat = Y.lvs %*% Path[endo == 1, ]

    inner_residuals = Y.lvs[, endo == 1] - Y.hat
    res = cbind(outer_residuals, inner_residuals)

    return(res)

}
