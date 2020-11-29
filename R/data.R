# 4-state partitioned survival model -------------------------------------------
#'  Example data for a 4-state partitioned survival model
#'
#' A collection of example datasets containing simulated survival, costs, and utility data for a 4-state
#' partitioned survival model.
#'
#' 
#' @format 
#' A list containing the following elements:
#' \itemize{
#' \item{Survival}{ A data frame containing patient information and time to 3 separate survival endpoints.}
#' \item{Costs}{A list of data frames. The first data frame contains medical cost data and the
#' second data frame contains drug cost data.}
#' }
#' 
#' @section Survival data:
#' The survival data frame contains a list of 3 survival curves, each containing the following columns.
#' \describe{
#' \item{female}{An indicator variable equal to 1 if the patient is female and 0 otherwise.}
#' \item{age}{The age of the patient in years.}
#' \item{strategy_id}{The id of the treatment strategy used.}
#' \item{endpoint1_time}{Follow up time with right censored data to survival endpoint 1.}
#' \item{endpoint1_status}{A status indicator for survival endpoint 1 equal to 0 if alive and 1 if dead.}
#' \item{endpoint2_time}{Follow up time with right censored data to survival endpoint 2.}
#' \item{endpoint2_status}{A status indicator for survival endpoint 2 equal to 0 if alive and 1 if dead.}
#' \item{endpoint3_time}{Follow up time with right censored data to survival endpoint 3.}
#' \item{endpoint3_status}{A status indicator for survival endpoint 3 equal to 0 if alive and 1 if dead.}
#' }
#' 
#' @section Cost data:
#' The cost list contains two data frames.The first data frame contains data on the
#' medical costs by patient and health state, and contains the following columns:
#' \describe{
#'  \item{patient_id}{An integer denoting the id of the patient.}
#'   \item{female}{An indicator variable equal to 1 if the patient is female and 0 otherwise.}
#'   \item{state_name}{A categorical variable denoting the three possible health states.}
#'   \item{costs}{Annualized medical costs.}
#'  }
#'  
#' The second data frame contains data on the drug costs associated with each treatment strategy.
#' \describe{
#'  \item{strategy_id}{The id of each treatment strategy.}
#'  \item{costs}{Annualized drug costs.}
#'  }
"psm4_exdata"

# 3-state reversible multi-state model -----------------------------------------
#' Example data for a reversible 3-state multi-state model
#'
#' Example multi-state data for parameterizing a continuous time state
#' transition model. Costs and utility 
#' are also included to facilitate cost-effectiveness analysis. 
#' @format 
#' A list containing the following elements:
#' \itemize{
#' \item{transitions}{ A data frame containing the times at which patient transitions between health states based
#' on the [prothr][mstate::prothr] dataset from the [mstate][mstate::mstate] package.}
#' \item{costs}{ A list of data frames. The first data frame contains summary medical cost estimates and the
#' second data frame contains drug cost data.}
#' \item{utility}{ A data frame of summary utility estimates.}
#' }
#' 
#' @section Transitions data:
#' The data frame has the following columns:
#' \describe{
#' \item{strategy_id}{Treatment strategy identification number.}
#' \item{patient_id}{Patient identification number.}
#' \item{age}{Patient age (in years).}
#' \item{female}{1 if a patient is female; 0 if male.}
#' \item{from}{Starting state.}
#' \item{to}{Receiving state.}
#' \item{trans}{Transition number.}
#' \item{Tstart}{Starting time.}
#' \item{Tstop}{Transition time.}
#' \item{years}{Elapsed years between `Tstart` and `Tstop`.}
#' \item{status}{Status variable; 1=transition, 0=censored.}
#' }
#' 
#' @section Cost data:
#' The cost list contains two data frames. The first data frame contains 
#' data on the drug costs associated with each treatment strategy.
#' \describe{
#'  \item{strategy_id}{The treatment strategy identification number.}
#'  \item{costs}{Annualized drug costs.}
#'  }
#' 
#' The second data frame contains summary data on
#' medical costs by health state, and contains the following columns:
#' \describe{
#'   \item{state_id}{The health state identification number.}
#'   \item{mean}{Mean costs.}
#'   \item{se}{Standard error of medical costs.}
#'  }
#'  
#' @section Utility data:
#' The data frame has the following columns:
#' \describe{
#' \item{state_id}{The health state identification number.}
#' \item{mean}{Mean utility}
#' \item{se}{Standard error of utility}
#' }
#' 
"mstate3_exdata"

# 3-state multinomial model ----------------------------------------------------
#' Example data for a 3-state multinomial model
#'
#' Example discrete time health state transitions data simulated using 
#' multinomial logistic regression. Costs and utility 
#' are also included to facilitate cost-effectiveness analysis. 
#'
#' @format 
#' A list containing the following elements:
#' \itemize{
#' \item{transitions}{ A data frame containing patient transitions between health
#' states at discrete time intervals (i.e., on a yearly basis).}
#' \item{costs}{ A list of data frames. The first data frame contains 
#' drug cost data and the second contains summary medical cost estimates.}
#' \item{utility}{ A data frame of summary utility estimates.}
#' }
#' 
#' @section Transitions data:
#' The data frame has the following columns:
#' \describe{
#' \item{patient_id}{Patient identification number.}
#' \item{strategy_id}{Treatment strategy identification number.}
#' \item{strategy_name}{Treatment strategy name.}
#' \item{age}{Patient age (in years).}
#' \item{age_cat}{A factor variable with 3 age groups: (i) age less than 40, 
#' (ii) age at least 40 and less than 60, and (iii) age at least 60.}
#' \item{female}{1 if a patient is female; 0 if male.}
#' \item{year}{The year since the start of data collection with the first year equal to 1.}
#' \item{state_from}{State making a transition from.}
#' \item{state_to}{State making a transition to.}
#' \item{year_cat}{Factor variable for year with 3 categories: (i) year 3 and below, (ii)
#' year between 3 and 6, and (iii) year 7 and above.}
#' }
#' 
#' @section Cost data:
#' The cost list contains two data frames. The first data frame contains 
#' data on the drug costs associated with each treatment strategy.
#' \describe{
#'  \item{strategy_id}{The treatment strategy identification number.}
#'   \item{strategy_name}{The treatment strategy name.}
#'   \item{costs}{Annualized drug costs.}
#'  }
#' 
#' The second data frame contains summary data on
#' medical costs by health state, and contains the following columns:
#' \describe{
#'   \item{state_id}{The health state identification number.}
#'   \item{state_name}{The name of the health state.}
#'   \item{mean}{Mean medical costs.}
#'   \item{se}{Standard error of medical costs.}
#'  }
#'  
#' @section Utility data:
#' The data frame has the following columns:
#' \describe{
#' \item{state_id}{The health state identification number.}
#' \item{state_name}{The name of the health state.}
#' \item{mean}{Mean utility}
#' \item{se}{Standard error of utility.}
#' }
"multinom3_exdata"

# 3-state oncology model -------------------------------------------------------
#' Multi-state oncology data for 3-state model
#'
#' Simulated 3-state dataset in oncology with three health states
#' (Stable, Progression, and Death) and three possible transitions (Stable ->
#' Progression, Stable -> Death, and Progression -> Death). 
#'
#' @format 
#' A `data.table` with the following columns:
#' \describe{
#' \item{from}{Health state making a transition from.}
#' \item{to}{Health state making a transition to.}
#' \item{female}{1 if a patient is female; 0 if male.}
#' \item{age}{Patient age (in years).}
#' \item{patient_id}{Patient identification number.}
#' \item{time_start}{Starting time.}
#' \item{time_stop}{Stopping time.}
#' \item{time}{Elapsed years between `time_start` and `time_stop`.}
#' \item{status}{Status indicator: 1=transition, 0=censored.}
#' \item{transition_id}{Integer denoting transition: 1 = Stable -> Progression,
#' 2 = Stable -> Death, 3 = Progression -> Death.}
#' \item{strategy_name}{Standard of care (SOC) or new treatment (New).}
#' }
#' 
#' @examples 
#' head(onc3)
"onc3"

#' Convert multi-state data to PFS and OS data
#' 
#' Convert a multi-state dataset with irreversible transitions containing 3 health 
#' states to a dataset with one row per patient and progression-free survival (PFS)
#' and overall survival (OS) time-to-event outcomes. 
#' 
#' @param data A multi-state dataset.
#' @param patient_vars Character vector of the names of patient specific variables.
#' @param transition Character string with the name of the variable identifying
#' a transition. The transition variable should be integer valued with values
#' 1, 2, and 3 for the Stable -> Progression, Stable -> Death, and 
#' Progression -> Death transitions, respectively. 
#' @param status Character string with the name of the status variable (1 = event,
#' 0 = censored).
#' @param time_stop Character string with the name of the stopping time variable
#' (i.e., time patient transitions from state \eqn{r} to state \eqn{s}).
#' 
#' @examples 
#' as_pfs_os(onc3, patient_vars = c("patient_id", "female", "age")) 
#' 
#' @return A `data.table` with one row per patient containing each variable in 
#' `patient_vars`  as well as a time variable and status indicator for both 
#' PFS (`pfs_status`, `pfs_time`) and OS (`os_time`, `os_status`). 
#' @export
as_pfs_os <- function(data, patient_vars, status = "status", time_stop = "time_stop",
                      transition = "transition_id") {
  data <- as.data.table(data)
  
  # Checks for possible transitions
  unique_transvar <- unique(data[[transition]])
  if (!all(unique_transvar %in%  c(1L, 2L, 3L))) {
    stop(paste0("'", transition, "'", " should be a vector with unique values c(1, 2, 3)."))
  }
  min_transvar_msg <- 
  if (length(unique_transvar) == 2) {
    if (!isTRUE(all.equal(unique_transvar, c(1L, 2L)))) {
      stop(paste0("If '", transition, "'", " contains 2 values, they should be 1 and 2."))
    }
  }
  if (length(unique_transvar) == 1) {
    stop(paste0("'", transition, "'", " should contain at a minimum values 1 and 2."))
  }

  # Cast wide
  f <- as.formula(
    paste0(paste(patient_vars, collapse=" + "),
          "~", 
          transition)
  )
  data <- dcast(data, 
                formula = f,
                value.var = c(status, time_stop))
  
  # Helper functions so with strings passed as data.table columns
  timev <- function(i) paste0(time_stop, "_", i)
  gtimev <- function(i) data[[timev(i)]]
  statusv <- function(i) paste0(status, "_", i)
  gstatusv <- function(i) data[[statusv(i)]]
  
  # Compute PFS endpoint
  set(data, j = "pfs_time", 
      value = pmin(data[[timev(1)]], data[[timev(2)]]))
  set(data, j = "pfs_status", 
      value = pmax(data[[statusv(1)]], data[[statusv(2)]]))

  # Works even if no one progresses
  if (is.null(gstatusv(3))) data[, (statusv(3)) := NA]
  if (is.null(gtimev(3))) data[, (timev(3)) := NA]
  
  # Compute OS endpoint
  data[, os_status := fcase(
    gstatusv(2) == 1 | gstatusv(3) == 1, 1L,
    gstatusv(2) == 0 & gstatusv(3) == 0, 0L,
    gstatusv(2) == 0 & is.na(gstatusv(3)), 0L
    )
  ]
  
  data[, os_time := fcase(
    gstatusv(2) == 1, gtimev(2), # Observed death 1 -> 3
    gstatusv(2) == 0 & is.na(gstatusv(3)), gtimev(2), # Right censored 1 -> 3
    gstatusv(2) == 0 & gstatusv(3) == 1, gtimev(3), # Observed death 2 -> 3
    gstatusv(2) == 0 & gstatusv(3) == 0,  gtimev(3) # Right censored 2 -> 3
    )
  ]
  
  # Clean and return
  data[, c(statusv(1:3), timev(1:3)) := NULL]
  return(data[, ])
} 