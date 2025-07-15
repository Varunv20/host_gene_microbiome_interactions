library(shiny)
library(bslib)
library(dplyr)
library(rsconnect)

## Lasso tutorial
## Sambhawa Priya

## Goal
## Our aim in this tutorial is to show how to run our lasso analysis on a small set of genes (~2-3 genes)
## to identify host gene-taxa associations for the demo dataset.


## Install and import libraries
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("glmnet","data.table","hdi","stabs") 
check.packages(packages)

############## Functions ###############

## Please execute these functions

load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
  genes <- as.matrix(genes)
  
}

load_microbiome_abnd <- function(filename){
  microbes <- read.table(filename,sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
  microbes <- as.matrix(microbes)
  
}

estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  
  
  ## Fit a lasso object
  lasso.fit = glmnet(x,y_i,alpha = 1) ## this is same as cv.fit$glmnet.fit from loocv code below.
  beta <- as.vector(coef(lasso.fit, s = bestlambda)) ## This gives coefficients of fitted model, not predicted coeff.
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  ## predicted coefficients, same as coefficient of fitted model lasso. Either one is fine.
  # beta = predict(lasso.fit,s=bestlambda, type="coef")
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1))
    sigma.flag = 0
  } else{
    sigma = 1 ## conservative option
    # sigma = ss_res ## lenient option
    sigma.flag = 2
  }
  
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}

fit.cv.lasso <- function(x, y_i, kfold){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  
  ## glmnet CV
  cv.fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE, standardize = T)  
  lambdas = data.frame(cv.fit$lambda,cv.fit$cvm)
  
  ## get best lambda -- lambda that gives min cvm
  bestlambda <- cv.fit$lambda.min
  bestlambda_index <- which(cv.fit$lambda == bestlambda)
  
  ## Get R^2 of final model
  final_model <- cv.fit$glmnet.fit
  r_sqr_final_model <- cv.fit$glmnet.fit$dev.ratio[bestlambda_index]
  
  ## Get adjusted R^2
  r_sqr_final_adj <- adj_r_squared(r_sqr_final_model, n = nrow(x), 
                                   p = sum(as.vector(coef(cv.fit$glmnet.fit, 
                                                          s = cv.fit$lambda.min)) > 0))
  
  return(list(bestlambda = bestlambda, r.sqr = r_sqr_final_model, 
              r.sqr.adj = r_sqr_final_adj
  ))
}
## Sparse CCA tutorial
## Sambhawa Priya

## Install and import libraries
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("PMA","data.table") 
check.packages(packages)

################# Functions ################

load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F, stringsAsFactors = F)
  genes <- as.matrix(genes)
  
}

load_microbiome_abnd <- function(filename){
  microbes <- read.table(filename,sep="\t",head=T, row.names = 1, check.names =F, stringsAsFactors = F)
  microbes <- as.matrix(microbes)
  
}

run_sparseCCA <- function(X, Z, CCA.K, penaltyX, penaltyZ, vInit=NULL, outputFile=NULL){
  CCA.out <-  CCA(X,Z,typex="standard",typez="standard",K=CCA.K,
                  penaltyx=penaltyX,penaltyz=penaltyZ,
                  v=vInit)
  if(!is.null(outputFile)){
    sink(outputFile)
    print(CCA.out)
    sink()
  }
  
  ## add rownames to output factors
  rownames(CCA.out$u) <- colnames(X)
  rownames(CCA.out$v) <- colnames(Z)
  ## compute contribution of selected features to each of the samples.
  CCA_var_genes <- X %*% CCA.out$u ## canonical variance for genes 
  CCA_var_microbes <- Z %*% CCA.out$v ## canonical variance for microbes
  
  return(list(CCA.out, CCA_var_genes, CCA_var_microbes))
  
}

get_avg_features <- function(cca_cov, CCA.K){
  num_features <- 0
  for(k in 1:CCA.K){
    num_features <- num_features + length(which(cca_cov[,k]!=0))
  }
  avg_features <- num_features/CCA.K
}

save_CCA_components <- function(CCA.out, CCA.K, dirname){
  ## Print canonical covariates in files 
  for(i in CCA.K){
    print(paste0("Writing significant component = ", i))
    selected_X <- which(CCA.out$u[,i]!=0) 
    selected_X <- rownames(CCA.out$u)[selected_X]
    coeff_X <- unname(CCA.out$u[selected_X,i])
    selected_Z <- which(CCA.out$v[,i]!=0)
    selected_Z <- rownames(CCA.out$v)[selected_Z]
    coeff_Z <- unname(CCA.out$v[selected_Z,i])
    ## Make all vectors of same length to avoid repetition of elements from shorter vectors.
    n <- max(length(selected_X), length(selected_Z))
    length(selected_X) <- n                      
    length(selected_Z) <- n
    length(coeff_X) <- n
    length(coeff_Z) <- n
    selected_XZ <- as.data.frame(cbind(gene = selected_X, gene_coeff = coeff_X,
                                       taxa = selected_Z, taxa_coeff = coeff_Z))
    write.table(selected_XZ, file=paste0(dirname,"gene_taxa_component_",i,".txt"), sep = "\t", col.names = NA)
  }
  
}

tune_params_grid_search <- function( X, Y){
  
  scoreXcv <- c()
  scoreYcv <- c()
  penaltyX <- seq(0.05,0.4,length=10)
  penaltyY <- seq(0.05,0.4,length=10)
  corr_demo <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
  num_samples <- nrow(genes)
  start_time <- Sys.time()
  for( i in 1:length(penaltyX)){
    for(j in 1:length(penaltyY)){
      # print(paste0("Index: i = ",i,", j =", j)); flush.console()
      for(k in 1:num_samples){
        
        print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()
        #compute weights with sample k held out:
        #Default niter = 15 edited to 5 to speed this up.
        res <- CCA(X[-k,],Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F, standardize = T)
        ## Compute scores for k'th sample for first pair of canonical variables
        ## Take weight of features (res$u and res$v) computed using all except 
        ## the kth sample and multiple it by values for the kth sample in the 
        ## feature matrix X and Y.
        scoreXcv[k] <- X[k,]%*%res$u ## single value
        scoreYcv[k] <- Y[k,]%*%res$v ## single value
      }
      ## correlation between scores for X and Y for all held out samples.
      corr_demo[i,j] = cor(scoreXcv,scoreYcv) 
    }
  }
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste0("Time elapsed for param tuning = ", time_elapsed)); flush.console()
  
  row.names(corr_demo) <- as.character(penaltyX)
  colnames(corr_demo) <- as.character(penaltyY)
  
  corr_demo_df <- as.data.frame(corr_demo)
  rownames(corr_demo_df)
  colnames(corr_demo_df)
  
  ##identify best penalty parameters
  # find index with max absolute corr
  bestpenalty <- which(abs(corr_demo) == max(abs(corr_demo)), arr.ind = TRUE)
  bestpenalty
  bestpenaltyX <- penaltyX[bestpenalty[1]]
  bestpenaltyY <- penaltyY[bestpenalty[2]]
  
  return (c(bestpenaltyX,bestpenaltyY))
}

test_significance_LOOCV <- function(X, Y, bestpenaltyX, bestpenaltyY, num_components){
  cca.k = num_components
  scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
  scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
  corr_pval <- c()
  corr_r <- c()
  for(i in 2:nrow(X)){ #n = no. of samples
    
    print("hello")
    print(dim(X))
    print(dim(Y))
    print(dim(X[-i,]))
    print(dim(Y[-i,]))
    
    #compute weights with sample i held out:
    res <- CCA(X[-i,],Y[-i,], typex="standard",typez="standard", penaltyx=bestpenaltyX, penaltyz=bestpenaltyY, K=cca.k)
 ## default niter = 15 which is spit out when trace = T (default)
    ###compute scores for i'th sample for each component (pair of canonical variables),typex="standard",typez="standard",K=CCA.K,
    for(j in 1:cca.k){
      #print(paste0("i = ", i," K = ", j)); flush.console()
      scoresXcv[i,j] <- X[i,]%*%res$u[,j]
      scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
    }
  }
  ## Test for each components
  for(j in 1:cca.k){
    corr <- cor.test(scoresXcv[,j],scoresYcv[,j]) ## Pearson correlation.
    corr_pval[j] <- corr$p.value
  }
  corr_pval
}



## functions to compute R2
r_squared <- function(y, yhat) {
  ybar <- mean(y)
  ## Total SS
  ss_tot <- sum((y - ybar)^2)
  ## Residual SS
  ss_res <- sum((y - yhat)^2)
  ## R^2 = 1 - ss_res/ ss_tot
  1 - (ss_res / ss_tot)
}
## Function for Adjusted R^2
## n sample size, p number of prameters
adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}

read_demo_output <- function(dir, bestpenaltyX, bestpenaltyY) {
  file_path <- paste0(dir,"/CCA_demo_output_",bestpenaltyX,"_",bestpenaltyY,".txt")
  file_size <- file.info(file_path)$size 
  entire_text <- readChar(file_path, file_size) 
  entire_text
}

run_SparseCCA_ALL <- function(gene_input, microbe_input, output_dir, bestpenaltyX, bestpenaltyY, k) {

  gene_input <- gene_input$datapath
  microbe_input <- microbe_input$datapath
  genes <- load_gene_expr(paste0(gene_input))
  # 
  microbes <- load_microbiome_abnd(paste0(microbe_input))
  # ## Ensure same sampleIDs in both genes and microbes data before sparse CCA
  stopifnot(all(rownames(genes) == rownames(microbes)))
  
  print(dim(genes))
  print(dim(microbes))
  # 
  # ## SKIP if using pre-computed values above
  # ## select tuning parameters
  # # bestPenalty <- tune_params_grid_search(genes,microbes)
  # # bestpenaltyX <- bestPenalty[1]
  # # bestpenaltyY <- bestPenalty[2]
  # 
  # #### Run sparse CCA
  # 
  # ## Set the number of desired components
  cca.k = k
  outputfile <- paste0(output_dir,"CCA_demo_output_",bestpenaltyX,"_",bestpenaltyY,".txt")
  print(outputfile)
  # 
  # ## Run sparse CCA using selected tuning param using permutation search
  cca <- run_sparseCCA(genes, microbes, cca.k, bestpenaltyX, bestpenaltyY,
                       outputFile=outputfile)
  # 
  # ## average number of genes and microbes in resulting components
  avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
  
  # 
  avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
  # 
  # #### Test significance of components using LOOCV
  CCA_pval <- test_significance_LOOCV(genes, microbes, bestpenaltyX, bestpenaltyY, cca.k)
  # 
  length(which(CCA_pval < 0.1))
  which(CCA_pval < 0.1)
  # 
  CCA_padj <- p.adjust(CCA_pval, method = "BH")
  
  # 
  length(which(CCA_padj < 0.1))
  which(CCA_padj < 0.1)
  
  # #### Output significant components
  sig_cutoff <- 0.1
  sig <- which(CCA_padj < sig_cutoff)
  dirname <- paste0(output_dir,"/demo_gene_taxa_components/")
  # ## This will return FALSE if the directory already exists or is uncreatable,
  # ## and TRUE if it didn't exist but was succesfully created.
  ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
  save_CCA_components(cca[[1]],sig,dirname)
  out_text = read_demo_output(output_dir,bestpenaltyX, bestpenaltyY)
  return(out_text)
  
}

run_Lasso_ALL <- function(genes_path, microbes_path, gene_num) {
  genes_path <- genes_path$datapath
  microbes_path <- microbes_path$datapath
  
  genes <- load_gene_expr(paste0(genes_path))
  dim(genes) #[1]    44 3
  # 
  # ## load microbiome data (note, here we load microbiome data with sex covariate)
  microbes <- load_microbiome_abnd(paste0(microbes_path))
  dim(microbes) #[1]  44 236 -- 235 taxa + 1 sex covariate
  # 
  # 
  ## Ensure same sampleIDs in both genes and microbes data before sparse CCA
  stopifnot(all(rownames(genes) == rownames(microbes)))
  # 
  y <- genes #response
  x <- microbes #predictors
  # 
  # 
  # ############ Fit lasso model and test inference using HDI ############
  # 
  # ## We are going to test three genes: WNT5A, RIPK3, and SMAP2 for their association with microbes
  # 
  # ## Extract expression of first gene in the matrix
  df_list <- list()
  
  for (i in 1:dim(y)[2]) {
    print(length(y))
    print(dim(y))
    y_i <- y[,i]
    gene_name <- colnames(y)[i]
    # 
    # ## Make sure y_i is numeric before model fitting
    stopifnot(class(y_i) == "numeric")
    # 
    # ## Fit lasso CV model
    fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i))
    bestlambda <- fit.model$bestlambda
    r.sqr <- fit.model$r.sqr ## note this will give us R^2 for the gene's final model fit using bestLambda
    ## This R^2 reflects final model R^2 for this gene using all the microbes in the model,
    # ## and does not correspond to each gene-microbe pair.
    # 
    # ## Estimate sigma and betainit using the estimated LOOCV lambda.
    # ## Sigma is the standard deviation of the error term or noise.
    sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
    sigma <- sigma.myfun$sigmahat
    beta <- as.vector(sigma.myfun$betahat)[-1] ## remove intercept term
    sigma.flag <- sigma.myfun$sigmaflag
    # 
    # ## Inference using lasso projection method, also known as the de-sparsified Lasso,
    # ## using an asymptotic gaussian approximation to the distribution of the estimator.
    lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
    # ## A few lines of log messages appear here along with a warning about substituting sigma value (standard deviation of error term or noise)
    # ## because we substituted value of sigma using our computation above.
    # # Warning message:
    # #   Overriding the error variance estimate with your own value.
    # 
    # ## get 95% confidence interval (CI)
    lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
    # 
    # ## prep lasso output dataframe
    lasso.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)),
                           taxa = names(lasso.proj.fit$pval.corr),
                           r.sqr = r.sqr,
                           pval = lasso.proj.fit$pval,
                           ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                           row.names=NULL)
    
    # 
    # ## sort by p-value
    lasso.df <- lasso.df[order(lasso.df$pval),]
    head(lasso.df)
    # 
    # ################# Stability selection #################
    # 
    # ## set a seed for replicability
    set.seed(0511)
    # 
    # ## perform stability selection using glmnet lasso
    stab.glmnet <- stabsel(x = x, y = y_i,
                           fitfun = glmnet.lasso, cutoff = 0.6,
                           PFER = 1)
    # 
    taxa.selected <- names(stab.glmnet$selected)
    if(length(taxa.selected) == 0) taxa.selected <-"None"
    
    
    stabsel.df <- data.frame("gene" = gene_name, "taxa" = taxa.selected)
    if(taxa.selected == "none"){
      stabsel.df$stability_selected = "no"
    }else stabsel.df$stability_selected = "yes"
    
    head(stabsel.df)
    
    # ################ Merge output of lasso+hdi and stabsel #################
    # 
    overlap_lasso_stabsel <- merge(lasso.df,stabsel.df, by = c("gene","taxa"))
    df_list[[i]] = overlap_lasso_stabsel
  }
  bind_rows
  return(bind_rows(df_list))
}

sparseCCA_page <- page_fluid(layout_columns( 
  card( 
    card_header("Data Input"),
    p("Gene Expression Data"),
    fileInput("gen_ex", "Upload file"),
    p("Microbiome Data"),
    fileInput("microbiome", "Upload file"),
    textInput( 
      "cca_file_output", 
      "Output Folder:"
      ), 
  ), 
  card( 
    card_header("Variables"),
    numericInput("bestpenaltyX", "bestpenaltyX", 0.05),
    numericInput("bestpenaltyY", "bestpenaltyY", 0.3222),
    numericInput("k", "Number Of Components", 10),
    numericInput("sig_cutoff", "sig_cutoff", 0.1),
  ), 
  ),
  actionButton("runSparse", "Run Sparse CCA"),
  verbatimTextOutput("SCCA")
)

lasso_page <- page_fluid(layout_columns( 
  card( 
    card_header("Data Input"),
    p("Gene Expression Data"),
    fileInput("gen_ex", "Upload file"),
    p("Microbiome Data"),
    fileInput("microbiome", "Upload file"),
  ), 
  card( 
    card_header("Variables"),
  ), 
),
actionButton("runLasso", "Run Lasso"),
tableOutput("LASSO") 
)
ui <- page_fluid(
  titlePanel("Host gene-microbiome associations"),
  navset_pill( 
    nav_panel("Lasso", lasso_page), 
    nav_panel("Sparse CCA", sparseCCA_page), 
   
  ), 
  
)

server <- function(input, output) {
  output$filepath <- renderPrint({
    req(input$file)  # Wait for file to be uploaded
    input$file$datapath  # This is the temp file path on the server
  })
  output$SCCA <- renderText({
    result_text <- run_SparseCCA_ALL(input$gen_ex, input$microbiome, input$cca_file_output, input$bestpenaltyX, input$bestpenaltyY, input$k)
    return(result_text)
  }) |> bindEvent(input$runSparse) 
  output$LASSO <- renderTable({
    result_text <- run_Lasso_ALL(input$gen_ex, input$microbiome, input$num)
    return(result_text)
  }, striped=TRUE) |> bindEvent(input$runLasso) 
}

shinyApp(ui, server)

