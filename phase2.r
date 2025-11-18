library(dplyr)
library(caret)
library(LogicReg)
library(Metrics)
library(foreach)
library(parallel)
library(tidyverse)
library(combinat)

simulate_data <- function() {
    #Simulate Raw Data
    x <- rgamma(10000, 2, 3)
    x <- matrix(x, nrow = 1000, ncol = 10) #each row is a person, each column is a feature

    #Normalize Raw Data
    x_sums <- rowSums(x)
    x_norm <- as.data.frame(x / x_sums)

    ### Create a phi vector (true cutoff)
    true_phi <- rep(0.05, 10) #simulate phi (true cutoff). If an element in x is lower than respective element in phi, then it is 0. If an element in x is higher than respective element in phi, then it is 1. 

    z <- mapply(function(x_col, phi) {
    # Apply z_to_x to each element of the column z_col with the corresponding phi
    sapply(x_col, function(x_element) x_to_z(x_element, phi))
    }, as.data.frame(x_norm), true_phi)

    #### Multiply betas with Zs
    w <- matrix(NA, 1000, 2)
    w[,1] <- z[,1] & z[,2] 
    w[,2] <- z[,3] | z[,4]
    # w[,3] <- z[,5] & z[,6]
    # w[,4] <- z[,7] | z[,8]
    # w[,5] <- z[,9] & z[,10]

    ### Generate regression coefficients (Betas)
    betas <- c(-5,5)

    dot_prod <- w %*% betas #dot product of z (confirmed that this works)

    ### Inverse Logit (go from \# to probability)
    p <- exp(dot_prod)/(1+exp(dot_prod))

    ### Go From Probability to Y
    y = rbinom(1000,1, p = p)

    df <- data.frame(cbind(x_norm,y))

    return(df)
}

logic_cv <- function(a, df, phi_test, folds_inner, trees, leaves) {
    train_cv <- df[-folds_inner[[a]],] #training folds
    test_cv <- df[folds_inner[[a]],] #test fold
    
    z_train <- mapply(function(z_col, phi_col) {
    # Apply z_to_x to each element of the column z_col with the corresponding phi
    sapply(z_col, function(z_element) x_to_z(z_element, phi_col))
    }, as.data.frame(train_cv[,1:10]), phi_test)
    
    z_test <- mapply(function(z_col, phi_col) {
    # Apply z_to_x to each element of the column z_col with the corresponding phi
    sapply(z_col, function(z_element) x_to_z(z_element, phi_col))
    }, as.data.frame(test_cv[,1:10]), phi_test)
    
    capture.output(mod_greed <- logreg(
    train_cv$y,
    z_train,
    type = 3,
    select = 6, ntrees = trees, nleaves = leaves), file = NULL) #type = logistic regression, select = greedy
    
    #PICK LAST MODEL (because it is associated with the trees/leaves we are testing)
    num_models <- length(mod_greed$alltrees) #finds the number of Models it made. We want to pick the last Model because it has the lowest score
    
    y_pred <- predict(mod_greed, newbin = z_test[,1:10])[[num_models]] #this returns the predicted values (probabilities) from the test data
    
    y_truth <- test_cv$y #real y 

    return(Metrics::logLoss(y_truth, y_pred)) #log loss on particular fold
}

pick_neighbor_to_j <- function(df, i,j, phis){
    ##PICKING A NEIGHBOR OF J
          phi_new <- phis[[i]] #make copy of old phi 
          xj <- df[,j] #jth column of x
          diff <- abs(phi_new[j]-xj) #get the differences between element of phi and each element in xj
          diffs_filt <- diff[diff > 0 & diff <= 0.30] #find the differences within threshold of 0.30

          if (length(diffs_filt) == 0) {
            rand_diff <- min(diff[diff > 0])  #Pick the smallest nonzero difference if there aren't any differences in the threshold 
          } 
          else 
          {
            rand_diff <- sample(diffs_filt, 1)  #Pick a difference within the threshold 
          }

          index_rand_diff <- which(diff == rand_diff) #find which element in difference vector has that random difference
          #there might be more than one index that has that random difference
          if (length(index_rand_diff) == 1) {
            candidate <- xj[index_rand_diff] 
          }
          else 
          {
            sampled <- sample(index_rand_diff,1) #take a random sample of those indices 
            candidate <- xj[sampled]
          }

          phi_new[j] <- candidate #changing jth element of phi_new to candidate

          return(phi_new)
}

sim_anneal_of_phi_vector_for_tree_leaf_combo_cv <- function(a, df, outer_folds, trees, leaves) { 
    training <- df[-outer_folds[[a]],] #training folds
    hold_out <- df[outer_folds[[a]],] #test fold

    inform <- paste("TESTING", trees, "TREE and", leaves, "leaves")
    print(inform)

    #initialize phi start such that each element is an actual x from their respective column in the phi vector
    make_phi_start <- function(x) {
        return(sample(x,1))
    }

    phi_start <- unname(apply(training[1:10], 2, make_phi_start))

    #START THE SIMULATED ANNEALING PROCESS WITH THE LOGIC REGRESSION
    #Initializing temperature, sigma, and i
    temp <- 0.5
    #if change in log loss is 0.01. after 10 iterations, the probability of accepting something worse is 0.50. 
    i <- 1

    log_loss_tracker <- c() #this stores the log losses
    phis <- list()

    phis[[1]] <- phi_start #initialize the phi start 

    set.seed(2024) #setting seed so we have the same folds every time
    folds_inner <- caret::createFolds(training$y, k = 10) #creating folds
    set.seed(NULL)

    patience <- 10 #Number of iterations needed to break the while loop
    tol <- 0.001
    switchy = TRUE

    while (switchy == TRUE) { #when probability gets smaller, or do like an epsilon check
      #add for loop for each temp settings run a bunch of iterations (allows more iterations at these larger temperatures and give you space to search) more iterations at the hotter temperature
        for (j in 1:10) { #for each element in phi 
        print(paste("i", i))
        print(paste("j", j)) #prints which iteration we are on 
        print(paste("temp",temp)) #prints the temperature we are on 

        phi_new <- pick_neighbor_to_j(training, i, j, phis)

        log_loss_phi_new <- mean(unlist(mclapply(1:10, function(a) {logic_cv(a, training, phi_new, folds_inner, trees, leaves)}, mc.cores = 10))) #test the current proposal phi 
        log_loss_phi_old <- mean(unlist(mclapply(1:10, function(a) {logic_cv(a, training, phis[[i]], folds_inner, trees, leaves)}, mc.cores = 10))) #test the previous phi
  
        if (log_loss_phi_new < log_loss_phi_old)  {
          #if the new loss is lower than previous loss
          phis[[i]][j] <- phi_new[j] #change the phi with new j (accept)
        }
        else {
          #if new loss is higher than or equal to previous loss
          p <- exp(-(log_loss_phi_new - log_loss_phi_old) / temp) #plug in new phi into p function. p is smaller when temp is higher? 
          unif_num <- runif(1)
          if (unif_num < p) {
            #if randomly generated number from uniform distribution  is less than p. when the temp is higher, p is also higher. so higher probability of being accepted. 
            phis[[i]][j] <- phi_new[j] #change the phi at j with new j (accept)
          } else {
            phis[[i]][j] <- phis[[i]][j] #we reject the new phi
          }
        }
        print(phis[[i]]) #prints phi after update of j element
        }

      print(paste("logloss",log_loss_tracker[length(log_loss_tracker)])) #prints the log loss 

      log_loss_post_j_updates <- mean(unlist(mclapply(1:10, function(a) {logic_cv(a, training, phis[[i]], folds_inner, trees, leaves)}, mc.cores = 10))) #calculating the log loss after updating all j elements in phi 

      log_loss_tracker <- append(log_loss_tracker,log_loss_post_j_updates)

      range_last_patience_losses <- range(log_loss_tracker %>% tail(patience)) #getting last 10 iterations
      diff_max_min <- range_last_patience_losses[2] - range_last_patience_losses[1]

      if ((i > 10) && diff_max_min <= tol)  { #if all consecutive differences are less than the tolerance
        switchy <- FALSE
      }

      #COOLING SCHEDULE
      temp <- 0.97 * temp #change temperature

      #UPDATE PHIS
      phis[[i+1]] <- phis[[i]] #updates the ith element of phi with correct phi 
      i <- i + 1 #next i
    } #END SIMULATED ANNEALING OF PHI VECTOR - NOW YOU HAVE A FINAL PHI VECTOR 

      final_phi <- phis[[length(phis)]] #find final phi

      #binarize the training with that final phi vector
      z_train_final_phi <- mapply(function(z_col, phi_col) {
      # Apply z_to_x to each element of the column z_col with the corresponding phi
      sapply(z_col, function(z_element) x_to_z(z_element, phi_col))
      }, as.data.frame(training[,1:10]), final_phi)

      #fit a final logistic regresison model on this dataframe
      capture.output(final_logic_model <- logreg(
      training$y,
      z_train_final_phi,
      type = 3,
      select = 6, ntrees = trees, nleaves = leaves), file = NULL) #type = logistic regression, select = greedy

      #PICK LAST MODEL (because it is associated with the trees/leaves we are testing)
      num_models <- length(final_logic_model$alltrees) #finds the number of Models it made. We want to pick the last Model because it has the lowest score

      #BINARIZE HOLD OUT WITH FINAL PHI
      z_hold_out_final_phi <- mapply(function(z_col, phi_col) {
      # Apply z_to_x to each element of the column z_col with the corresponding phi
      sapply(z_col, function(z_element) x_to_z(z_element, phi_col))
      }, as.data.frame(hold_out[,1:10]), final_phi)
      
      hold_out_preds <- predict(final_logic_model, newbin = z_hold_out_final_phi[,1:10])[[num_models]] #this returns the predicted values (probabilities) from the hold out set
      
      hold_out_truth <- hold_out$y #real y 
      
      #compare predictions to truth
      hold_out_loss <- Metrics::logLoss(hold_out_truth, hold_out_preds)

      logic_expr <- final_logic_model$alltrees[[num_models]]
      logic_expr <- capture.output(logic_expr)

      return(c(hold_out_loss, logic_expr))
} 

sim_anneal_of_phi_vector_for_tree_leaf_combo_no_cv <- function(df, winner_tree, winner_leaf) { 
      inform <- paste("BUILDING A FINAL MODEL USING", winner_tree, "TREE and", winner_leaf, "leaves")
      print(inform)

      #initialize phi start such that each element is an actual x from their respective column in the phi vector
      make_phi_start <- function(x) {
          return(sample(x,1))
      }
      phi_start <- unname(apply(df[1:10], 2, make_phi_start))

      #START THE SIMULATED ANNEALING PROCESS WITH THE LOGIC REGRESSION
      #Initializing temperature, sigma, and i
      temp <- 0.5
      #if change in log loss is 0.01. after 10 iterations, the probability of accepting something worse is 0.50. 
      i <- 1

      log_loss_tracker <- c() #this stores the log losses
      phis <- list()

      phis[[1]] <- phi_start #initialize the phi start 

      set.seed(2024) #setting seed so we have the same folds every time
      folds_inner <- caret::createFolds(df$y, k = 10) #creating folds
      set.seed(NULL)

      patience <- 10 #Number of iterations needed to break the while loop
      tol <- 0.001
      switchy = TRUE

      while (switchy == TRUE) { #when probability gets smaller, or do like an epsilon check
        #add for loop for each temp settings run a bunch of iterations (allows more iterations at these larger temperatures and give you space to search) more iterations at the hotter temperature
          for (j in 1:10) { #for each element in phi 
          print(paste("i", i))
          print(paste("j", j)) #prints which iteration we are on 
          print(paste("temp",temp)) #prints the temperature we are on 

          phi_new <- pick_neighbor_to_j(df, i, j, phis)

          log_loss_phi_new <- mean(unlist(mclapply(1:10, function(a) {logic_cv(a, df, phi_new, folds_inner, winner_tree, winner_leaf)}, mc.cores = 10))) #test the current proposal phi 
          log_loss_phi_old <- mean(unlist(mclapply(1:10, function(a) {logic_cv(a, df, phis[[i]], folds_inner, winner_tree, winner_leaf)}, mc.cores = 10))) #test the previous phi
    
          if (log_loss_phi_new < log_loss_phi_old)  {
            #if the new loss is lower than previous loss
            phis[[i]][j] <- phi_new[j] #change the phi with new j (accept)
          }
          else {
            #if new loss is higher than or equal to previous loss
            p <- exp(-(log_loss_phi_new - log_loss_phi_old) / temp) #plug in new phi into p function. p is smaller when temp is higher? 
            unif_num <- runif(1)
            if (unif_num < p) {
              #if randomly generated number from uniform distribution  is less than p. when the temp is higher, p is also higher. so higher probability of being accepted. 
              phis[[i]][j] <- phi_new[j] #change the phi at j with new j (accept)
            } else {
              phis[[i]][j] <- phis[[i]][j] #we reject the new phi
            }
          }
          print(phis[[i]]) #prints phi after update of j element
          }

        print(paste("logloss",log_loss_tracker[length(log_loss_tracker)])) #prints the log loss 

        log_loss_post_j_updates <- mean(unlist(mclapply(1:10, function(a) {logic_cv(a, df, phis[[i]], folds_inner, winner_tree, winner_leaf)}, mc.cores = 10))) #calculating the log loss after updating all j elements in phi 

        log_loss_tracker <- append(log_loss_tracker,log_loss_post_j_updates)

        range_last_patience_losses <- range(log_loss_tracker %>% tail(patience)) #getting last 10 iterations
        diff_max_min <- range_last_patience_losses[2] - range_last_patience_losses[1]

        if ((i > 10) && diff_max_min <= tol)  { #if all consecutive differences are less than the tolerance
          switchy <- FALSE
        }

        #COOLING SCHEDULE
        temp <- 0.97 * temp #change temperature

        #UPDATE PHIS
        phis[[i+1]] <- phis[[i]] #updates the ith element of phi with correct phi 
        i <- i + 1 #next i
      } #END SIMULATED ANNEALING OF PHI VECTOR - NOW YOU HAVE A FINAL PHI VECTOR 

      final_phi <- phis[[length(phis)]] #find final phi

      #binarize the entire df with that final phi vector
      z_train_entire_data <- mapply(function(z_col, phi_col) {
      # Apply z_to_x to each element of the column z_col with the corresponding phi
      sapply(z_col, function(z_element) x_to_z(z_element, phi_col))
      }, as.data.frame(df[,1:10]), final_phi)

      #fit a final logic regression model on this dataframe
      capture.output(final_logic_model <- logreg(
      df$y,
      z_train_entire_data,
      type = 3,
      select = 6, ntrees = winner_tree, nleaves = winner_leaf), file = NULL) #type = logistic regression, select = greedy

      #PICK LAST MODEL (because it is associated with the trees/leaves we are testing)
      num_models <- length(final_logic_model$alltrees) #finds the number of Models it made. We want to pick the last Model because it has the lowest score
      
      train_preds <- predict(final_logic_model, newbin = z_train_entire_data[,1:10])[[num_models]] #this returns the predicted values (probabilities) from the test data
      
      train_truth <- df$y #real y 
      
      #compare predictions to truth
      training_loss <- Metrics::logLoss(train_truth, train_preds)

      logic_expr <- final_logic_model$alltrees[[num_models]]
      logic_expr <- capture.output(logic_expr)

      return(list(winner_tree, winner_leaf, logic_expr, final_phi)) #log loss on particular fold
    } 

#OUTPUT OF THIS IS A LIST OF T/L AND THEIR CV LOSSES
run_tree_leaf_combo <- function(tl_row, df){
  set.seed(2024) #setting seed so we have the same folds every time
  outer_folds <- caret::createFolds(df$y, k = 2) #creating folds
  set.seed(NULL)

  tree <- tree_leaves_combos[tl_row,][1,1]
  leaf <- tree_leaves_combos[tl_row,][1,2]
  
  inform <- paste("TESTING", tree, "TREE and", leaf, "leaves")
  print(inform)

  output <- mclapply(1:2, function (a) {sim_anneal_of_phi_vector_for_tree_leaf_combo_cv(a, df, outer_folds, tree, leaf)}, mc.cores = 2) 

  results_list <- list(tree, leaf, output[1], output[2])

  return(results_list)

}
start.time <- Sys.time()

args <- commandArgs(trailingOnly=TRUE)
dataset_id <- as.numeric(args[1])

#getting the tree and leaf combinations
t <- c(1,2,3,4,5)
l <- c(2,4,6)
tree_leaves_combos <- expand.grid(t,l) %>% filter(Var2 >= Var1) %>% arrange(Var1) 

#Defining a very important binarizing function
x_to_z <- function(x,phi){
#if x > phi --> return 1
# if x < phi --> return 0
ifelse(x >= phi, 1, 0)
}

#SIMULATING DATA
df <- simulate_data()

#RUNNING FUNCTIONS 
#Testing tree and leaf parameters and reporting their cv losses
results_list_per_run <- mclapply(1:11, function (tl_row) {run_tree_leaf_combo(tl_row, df)}, mc.cores = 11)
results_list_per_run_df <- data.frame(matrix(unlist(results_list_per_run), nrow=length(results_list_per_run), byrow=TRUE))
colnames(results_list_per_run_df) <- c("Tree", "Leaf", "HoldOut1", "LE1", "HoldOut2", "LE2")
results_list_per_run_df$HoldOut1 <- as.numeric(results_list_per_run_df$HoldOut1)
results_list_per_run_df$HoldOut2 <- as.numeric(results_list_per_run_df$HoldOut2)

#FIND T/L WITH LOWEST HOLD OUT LOSS
winner_ms <- results_list_per_run_df %>% mutate("AverageLoss" = rowMeans(select(results_list_per_run_df, starts_with("H")), na.rm = TRUE))  %>% arrange(AverageLoss) %>% head(1) %>% select(Tree, Leaf)
tree_winner <- as.numeric(winner_ms[1,1])
leaf_winner <- as.numeric(winner_ms[1,2])

#TRAIN LOGIC REG MODEL WITH WINNER T/L ON THE ENTIRE DATA SET 
final_logic_exp_and_phi <- sim_anneal_of_phi_vector_for_tree_leaf_combo_no_cv(df, tree_winner, leaf_winner)

cat("Saving results for dataset", dataset_id, "\n")

#SAVING RESULTS
saveRDS(results_list_per_run_df, sprintf("/home/shanbhagn123/phase2/MS_evals_phase2/find_MS_%03d.rds", dataset_id))
saveRDS(final_logic_exp_and_phi, sprintf("/home/shanbhagn123/phase2/winner_LEs_phase2/winner_LE_%03d.rds", dataset_id))
end.time <- Sys.time()
print(end.time-start.time)
