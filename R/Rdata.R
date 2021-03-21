#' an example data set used in a CPUE analyst. Locations have been randomised for sensitivity so there is no 
#' sensitive information. The purpose is to demonstrate package functionality
#'
#' A dataset containing tlocation and measured variables from fishing activities for a species.
#'
#' @format A data frame with 4562 rows and 17 variables:
#' \describe{
#'   \item{primary_method}{gear type}
#'   \item{target_species}{species code that was reported as targetted}
#'   \item{fishing_duration}{species code that was reported as targetted}
#'   \item{effort_depth}{fishing depth}
#'   \item{effort_height}{height of net}
#'   \item{effort_width}{trawl door width}
#'   \item{start_latitude}{start latitude position decimal degrees}
#'   \item{start_longitude}{start longitude position decimal degrees}
#'   \item{vessel_key}{vessel identifier}
#'   \item{BAR_catch}{catch in kgs}
#'   \item{start_stats_area_code}{statistical area (discrete spatial variable)}

#'   ...
#' }
"cpue_df"

