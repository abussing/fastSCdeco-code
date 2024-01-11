
source(file = "/work/abussing/fast_sc_deco.R")



chunk2 <- function(x, n) {
  split(x, cut(seq_along(x), n, labels = FALSE))
}



metadat <- readRDS("/work/abussing/metadat_4471")
countmat <- readRDS("/work/abussing/countmat_4471")


xgene <- "COL12A1"

x_use <- mean(countmat[xgene, ])

###############  create gene pairs


gene_pairs <- t(combn(rownames(countmat)[-which(rownames(countmat) == xgene)], 2))


from_where <- 1
to_where <- 583740


# Let's break up this gene_pairs into subsets
nsubs <- 100

gene_idcs_list <- chunk2(x = from_where:to_where, n = nsubs)

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
taskid <- as.numeric(taskid)

gene_pairs <- matrix(gene_pairs[gene_idcs_list[[taskid]], ], ncol = 2)



############### for each gene pair, fit the model

countr <- 0
problemcountr <- 0
printcountr <- 0
storeres <- list()
problempairs <- list()


for (i in 1:nrow(gene_pairs)) {
  y1gene <- gene_pairs[i, 1]
  y2gene <- gene_pairs[i, 2]

  print(c(y1gene, y2gene))

  moddf <- data.frame(t(countmat[c(y1gene, y2gene, xgene), ]))

  colnames(moddf) <- c("y1", "y2", "x")

  moddf$celltype_longname <- metadat[metadat$X %in% rownames(moddf), "celltype_minor"]
  moddf$celltype <- 1 * (moddf$celltype_longname == "CAFs myCAF-like")

  moddf$celltype_times_x <- moddf$celltype * moddf$x

  eq1 <- y1 ~ x + celltype + celltype_times_x
  eq2 <- y2 ~ x + celltype + celltype_times_x
  eq3 <- ~1
  eq4 <- ~1
  eq5 <- ~ x + celltype + celltype_times_x
  eqlist <- list(eq1, eq2, eq3, eq4, eq5)


  gjrm_out <- gjrm(formula = eqlist, data = moddf, copula = "N", model = "B", margins = c("NBI", "NBI"), gamlssfit = TRUE, rinit = 50, rmax = 10000)


  if (!gjrm_out$fit$converged) {
    problemcountr <- problemcountr + 1
    problempairs[[problemcountr]] <- c(y1gene, y2gene)
  } else {
    Vb <- gjrm_out$Vb

    upper <- gjrm_out$coefficients + 1.96 * sqrt(diag(Vb))
    lower <- gjrm_out$coefficients - 1.96 * sqrt(diag(Vb))

    CImat <- cbind(lower, gjrm_out$coefficients, upper)

    numparams <- length(gjrm_out$coefficients)

    # tau1 and tau3 portions of the Vb matrix
    Vb_cut <- Vb[c(numparams - 4, numparams - 2), c(numparams - 4, numparams - 2)]
    sum_se <- sqrt(matrix(c(1, 1), nrow = 1) %*% Vb_cut %*% matrix(c(1, 1), ncol = 1))

    sum_lower <- sum(gjrm_out$coefficients[c(numparams - 4, numparams - 2)]) - 1.96 * sum_se
    sum_upper <- sum(gjrm_out$coefficients[c(numparams - 4, numparams - 2)]) + 1.96 * sum_se

    CImat <- rbind(c(sum_lower, sum(gjrm_out$coefficients[c(numparams - 4, numparams - 2)]), sum_upper), CImat)

    if (sign(lower[numparams - 2]) == sign(upper[numparams - 2])) {
      countr <- countr + 1

      print("okay")
      print("okay okay")
      print("okat okay okay")
      storeres[[countr]] <- list(c(y1gene, y2gene), CImat)
      print("okay okay okay")
      print("okay okay")
      print("okay")
      if (gjrm_out$fit$e.v < 0) {
        storeres[[countr]] <- append(storeres[[countr]], "hessian is not positive semidefinite")
      }

      if (gjrm_out$fit$gradi > 10) {
        storeres[[countr]] <- append(storeres[[countr]], "maximum gradient value is > 10")
      }
    }
  }

  if (printcountr < countr) {
    printcountr <- printcountr + 1
    print(storeres[[countr]])
  }
}


saveRDS(storeres, paste("/work/abussing/meeting_2024-01-12/RDS files/allgenes_from", from_where, "_to", to_where, " ", taskid, sep="" ))

saveRDS(problempairs, paste("/work/abussing/meeting_2024-01-12/RDS files/problemgenes_from", from_where, "_to", to_where, " ", taskid, sep=""))
