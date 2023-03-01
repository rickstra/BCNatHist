#' @title MatchDataVariabes
#' @description Match a data set to the core variables used for estimation and prediction. 
#'
#' @param data Dataframe containing the data. One entry per individual.
#' @param base_variables List linking variable names of \code{data} to variables used in the prediction: 
#' \describe{
#'   \item{entry}{The age at the start of the prediction.}
#'   \item{exit}{The age at end of prediction follow-up.}
#'   \item{is_case}{BC case status (case = 1).}
#'   \item{mode}{Mode of detection (screen-detected = 1).}
#'   \item{size}{tumor size (mm).}
#'   \item{scr_hist}{Screening histories (each entry a list of ages screened).}
#' }
#' @param include_prevalent_screen Indicates if a screening performed at the time of entry should be included in the follow-up. 
#'
#' @return \code{data.frame} with the following variables matched from the input data:
#' \describe{
#'   \item{$entry}{The age at the start of the prediction.}
#'   \item{$exit}{The age at end of prediction follow-up.}
#'   \item{$case}{BC case status (case = 1).}
#'   \item{$mode}{Mode of detection (screen-detected = 1).}
#'   \item{$v}{tumor volume (mm3).}
#'   \item{$scr}{Screening histories (each entry a list of ages screened).}
#'   \item{$e_scr}{Screening histories before study entry (each entry a list of ages screened).} }
#' @examples
#' 
#' @export
MatchDataVariables <- function(data, 
                        base_variables = list(
                          entry = "entry", exit = "exit", is_case = "case", 
                          mode = "mode", size = "size", scr_hist = "scr"
                        ), include_prevalent_screen = TRUE) {
  `%<%` <- ifelse(include_prevalent_screen,`<`, `<=`)
  out <- data.frame(
    entry = data[, ifelse(is.null(base_variables$entry), "entry", base_variables$entry)],
    exit = data[, ifelse(is.null(base_variables$exit), "exit", base_variables$exit)],
    case = data[, ifelse(is.null(base_variables$is_case), "case", base_variables$is_case)],
    mode = data[, ifelse(is.null(base_variables$mode), "mode", base_variables$mode)],
    v = (pi / 6) * data[, ifelse(is.null(base_variables$size), "size", base_variables$size)] ^ 3
  )
  out$scr <- data[, ifelse(is.null(base_variables$scr_hist), "scr", base_variables$scr_hist)]
  out$e_scr <- lapply(1:nrow(data), function(i) {
    data$scr[[i]][data$scr[[i]] %<% data$entry[i]] # Note: include entry screen?
  })
  out
}