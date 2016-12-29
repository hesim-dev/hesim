#' Simulate state transitions in markov cohort model
#'
#' This function simulates state transitions in a Markov cohort model. Allows for the simulation of multiple
#' cohorts. The model is fully deterministic and is appropriate for estimating population averages
#' in each cohort.
#'
#' @param z0 Matrix of starting distribution. One row for each cohort and one column
#' for each model state including abosrbing (i.e. death) states. If mortadj = TRUE, then
#' z0 should have one column for each model state in pmat and an additional column for the
#' state of death.
#'
#' @param ncycles Vector with the ith element indicating the number cycles to run the model for
#' the ith row in z0.
#'
#' @param pmat Array of distinct transition probability matrices.
#'
#' @param p_index Index vector denoting the probability matrix to apply to each cohort. If NULL,
#' the 3rd dimension of pmat must equal number of rows in z0.
#'
#'@param mortadj Logical. If TRUE, then the a death state is added to the probability matrix
#' and probabilities are adjusted accordingly. If FALSE, then probability matrix remains as is.
#'
#' @param mortprob Array of distinct mortality probability matrices. Each slice of array contains
#' a matrix with rows denoting model cycle (or patient age) and columns representing model state.
#'  Each element in the matrix is the probability of mortality for a particular group.
#' For example, a particular element might contain the probability of mortality for males
#' (3rd dimension of array) at age 55 (1st dimension of array) in model state 1 (2nd dimension
#' of array).
#'
#' @param mortprob_index Index vector denoting the mortality probability matrix to apply to each
#' cohort. If NULL, the 3rd dimension of mortprob must equal the number of rows in z0.
#'
#' @details The code is written in c++ to minmize simulation time.
#'
#' @return Matrix with the number of simulated individuals in each model state for each cohort
#' during each cycle.
#'
#' @export
markov_cohort_trans <- function(z0, ncycles, pmat, pmat_index = NULL,
                                mortadj = FALSE, mortprob = NULL, mortprob_index = NULL){
  # pmat
  if (is.null(pmat_index)){
    if(dim(pmat)[3] != nrow(z0)){
      stop("If pmat_index = NULL, then the number of matrices in pmat must equal number of
           rows in z0.")
    }
  }

  # mortality adjustment
  if(mortadj == TRUE){
    if(is.null(mortprob_index)){
      if(dim(mortprob)[3] != nrow(z0)){
        stop("If mortprob_index = NULL, then the number of matrices in mortprob must equal number
              of rows in z0.")
      }
    }
  }

  # c++ function
  return(markov_cohort_transC(z0, ncycles, pmat, pmat_index,
                              mortadj, mortprob, mortprob_index))
}

#' @export
array_convert <- function(x){
  if (is.vector(x)){
    x <- array(x, dim = c(length(x), 1, 1))
  }
  if (is.matrix(x)){
    x <- array(x, dim = c(nrow(x), ncol(x), 1))
  }
  return(x)
}

