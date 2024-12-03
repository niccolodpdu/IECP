#' ITEA: Iterative equilibration of cell-type asymmetry for enhanced deconvolution
#' 
#' Cell-type asymmetry in expression profiles can substantially impact the outcomes of reference-free deconvolution. 
#' To address this, we proposed ITEA which aims to enhance the accuracy of proportion inference under the influence of 
#' cell-type asymmetry. ITEA iteratively identifies Consistently Expressed Genes (CEGs) across the cell-type expression 
#' profiles in S after each round of deconvolution and uses these CEGs to normalize the bulk expression matrix X.
#' 
#' @name ITEA
#' @title Iterative equilibration of cell-type asymmetry
#' @author Dongping Du
#' @description Main function of ITEA
#' @param X Input matrix, with samples on columns and features on rows.
#' @param transpose Boolean, default False. Some deconvolution methods require the input to have features on columns; if set 
#' to True, the input matrix for deconvolution will be transposed before being passed to the function.
#' @param K The number of sources/cell types/sub-types
#' @param iCEG_thres Threshold of the cosine similarity of potential iCEGs with the reference vector (1,1...1,1). 
#' Features with a higher cosine value than the 'iCEG_thres' will be used as iCEGs for normalizing X
#' @param iteration The number of iteration for ITEA The function will stop once the number of iterations reaches this number
#' @param removal Optional, number of rounds to be removed. Remove the results of first several rounds from the final output.
#' @param deconvo_func The deconvolution function used in ITEA.
#' @param deconvo_param A list of parameters used in the deconvolution function specified by the user. Please note that the
#' input matrix should be the first element of the list and its value should be set identical with X.
#' @param feat_sele_func Optional, a function used for feature selection. If both 'feat_sele_func' and 'feat_sele_param' are
#' provided, ITEA will commence feature selection prior to deconvolution in each round.
#' @param feat_sele_param Optional, a list of parameters used in the feature selection function. Please note that the input 
#' matrix should be the first element of the list and its value should be set identical with X.
#' @param A_deconvo The address of A matrix (mixing proportion matrix) from the deconvolution output. It is used for extracting
#' A matrix from that output with the help of the 'extract_element' function
#' @param S_deconvo Optional. The address of S matrix (source matrix) from the deconvolution output. If not provided, ITEA
#' will use NNLS to retrieve the source matrix S from the input matrix X and the proportion matrix A
#' 
#' @returns A list of 1. A matrix (mixing proportion matrix) from each iteration of ITEA 2. S matrix (source matrix) from each
#' iteration of ITEA 3. The index of the optimal output from the iterative search, based on the minimum reconstruction error 
#' 4. The index of the optimal output from the iterative search, based on the maximum number of iCEGs found in X and S
#' 
#' @export

ITEA<-function(X=NULL, transpose=FALSE, K, iCEG_thres=0.98, iteration=20, removal=NULL, deconvo_func=NULL, 
                deconvo_param=NULL,feat_sele_func=NULL,feat_sele_param=NULL,A_deconvo=NULL, S_deconvo=NULL){
  
  # Scale to (0,1]
  input<-X/max(X)
  
  if (transpose){
    deconvo_param[[1]] <- t(input)
  } else {
    deconvo_param[[1]] <- input
  }
  
  # Feature selection, if a feature selection function & its parameters are provided:
  if (!is.null(feat_sele_func) & !is.null(feat_sele_param)){
    
    # Update the input
    feat_sele_param[[1]]<-input
    
    # Pass the parameters to the function
    feat<-do.call(feat_sele_func,feat_sele_param)
    
    # 'input_short' as the selected features
    input_short<-input[feat,]
    input_short<-input_short/max(input_short)
    
    # Pass selected features into deconvolution function
    deconvo_param[[1]] <- input_short
  }
  
  # Pass the parameters to the function
  result<-do.call(deconvo_func,deconvo_param)
  
  # Find A by the user-specified A address 
  A<-extract_element(result,A_deconvo)
  
  # Make sure A is slim and tall
  if (dim(A)[1]<dim(A)[2]){
    A<-t(A)
  }
  
  # If the address of S is not specified, use the following to get S
  if (is.null(S_deconvo)){
    S <- tryCatch({
      t(NMF::.fcnnls(A, t(X))$coef) # Faster way, but if not working:
    }, error = function(msg) {
      return(t(apply(X, 1, function(x) (nnls::nnls(A,x))$x))) # Use the safer way
    })
  } else {
    # Find S by the user-specified S address 
    S<-extract_element(result,S_deconvo)
  }
  
  # Make sure S is slim and tall too, for the purpose of cosine calculation
  if (dim(S)[1]<dim(S)[2]){
    S<-t(S)
  }
  
  # Calculate the cosine values between S column vectors (features) with the reference vectors (1,1...1,1)
  cos_S<-cos_iCEG(S)
  
  # Find the index of iCEGs based on cosine values
  index_S_iCEG<-which(cos_S>iCEG_thres)
  
  # Initiate the  output
  num_iceg<-NULL
  error<-NULL
  A_list<-list(NULL)
  S_list<-list(NULL)
  
  #### Iterations Start ####
  
  for (q in 1:iteration){
    # Use the index of iCEGs to locate the iCEGs in X
    iCEG_common<-X[index_S_iCEG,]
    
    # Commence totalcount normalization on iCEGs across all samples
    iCEG_norm <- totalcount(iCEG_common)
    
    # Retrieve normalization factors (on samples)
    factor_cosbin <- iCEG_norm$norm_factor
    
    # Initiate X_norm, the output of normalization
    X_norm<-X
    
    for (i in 1:dim(X)[2]) {
      # Apply the normalization factors on each sample
      X_norm[, i] <- X[, i]/factor_cosbin[i]
    }
    
    # Scaling
    input<-X_norm/max(X_norm)
    
    sum(is.na(input))
    # Update the input for deconvolution
    if (transpose){
      deconvo_param[[1]] <- t(input)
    } else {
      deconvo_param[[1]] <- input
    }
    
    # Feature selection, same as before
    if (!is.null(feat_sele_func) & !is.null(feat_sele_param)){
      feat_sele_param[[1]] <- input                           # Update the input
      feat<-do.call(feat_sele_func,feat_sele_param)           # Feature selection
      input_short<-input[feat,]                               # Selecting & scaling
      input_short<-input_short/max(input_short)
      deconvo_param[[1]] <- input_short                       # Pass selected features into deconvolution function
    }
    
    # Pass the parameters to the function
    result<-do.call(deconvo_func,deconvo_param)
    
    # Same as before, retrieve A and S
    A<-extract_element(result,A_deconvo)
    if (dim(A)[1]<dim(A)[2]){
      A<-t(A)
    }
    
    if (is.null(S_deconvo)){
      S <- tryCatch({
        t(NMF::.fcnnls(A, t(input))$coef)
      }, error = function(msg) {
        return(t(apply(input, 1, function(x) (nnls::nnls(A,x))$x)))
      })
    } else {
      S<-extract_element(result,S_deconvo)
    }
    
    # Calculate the cosine values, get the iCEG index
    cos_S<-cos_iCEG(S)
    index_S_iCEG<-which(cos_S>iCEG_thres)
    
    # Record the absolute error in each round
    error_round<-sum(abs(X_norm-t(A%*%t(S)))^2)
    error<-c(error,error_round)
    
    # Record the number of iCEG in each round
    num_iceg<-c(num_iceg,length(index_S_iCEG))
    
    # Record the A and S in each round
    A_list[[q]]<-A
    S_list[[q]]<-S
    
  }
  
  # If the optional removal of first several round(s) is true
  if (is.null(removal)==FALSE){
    # Remove the first several round(s)
    error<-error[(removal+1):iteration]
    num_iceg<-num_iceg[(removal+1):iteration]
    A_list<-A_list[(removal+1):iteration]
    S_list<-S_list[(removal+1):iteration]
  }
  
  # Optimal result, based on error
  output_index_min_error<-which.min(error)
  
  # Optimal result, based on number of iCEG
  output_index_max_iceg<-which.max(num_iceg)
  
  # Combine the outputs
  output<-list(A_list,S_list,output_index_min_error,output_index_max_iceg)
  
  # Rename the outputs
  names(output)<-c("A list","S list","Min error index","Max ICEG index")
  
  return(output)
}
