#' extract_and_merge
#' 
#' @details an function for pulling out parameters from a list object that has been created from obj$report(). Used to create a data.frame
#' with multiple models.
#' @param par_label character string of the parameter you want to extract
#' @param reports A list of lists. Each element is obj$report(), this will allow for multiple models
#' @param report_labels when creating the the data frame this will be the label representing each element in the list
#' @importFrom reshape2 melt
#' @return data.frame of quantities for each model
#' @export
#' @examples
#' \dontrun{
#'  n_sims = 10
#'  sim_ls = list()
#'  for(i in 1:n_sims)
#'    sim_ls[[i]] = obj$report()
#'  reports = list(sim_ls)  
#'  extract_and_merge("sigma", reports, report_labels = "sim")  
#' }
#' 
extract_and_merge <- function(par_label, reports, report_labels) {
  if(class(reports) != "list")
    stop("reports needs to be a list")
  
  n_lists = length(reports)
  if(!is.null(report_labels)) {
    if(length(report_labels) != n_lists)
      stop("if you supply report_labels there needs to be one for each element in reports")
  }
  # can only deal with scalar or vector parameters
  temp_par = get(par_label, reports[[1]][[1]])
  par_length = length(temp_par)
  df = NULL;
  for(i in 1:n_lists) {
    if(par_length == 1) {
      this_est = Reduce(c, lapply(reports[[i]], FUN = function(x) {get(par_label,x)}))
      lab = rep(report_labels[i], length(this_est))
      df = rbind(df, data.frame(lab = lab, value = this_est))
    } else {
      this_est = Reduce(cbind, lapply(reports[[i]], FUN = function(x) {get(par_label,x)}))
      # matrix n_par x n_sim
      rownames(this_est) = as.character(1:nrow(this_est))
      melted = reshape2::melt(this_est)
      #table(melted$Var1)
      melted$lab = report_labels[i]
      df = rbind(df, data.frame(lab = report_labels[i], value = melted$value, element = melted$Var1))
    }
  }
  return(df)
}
