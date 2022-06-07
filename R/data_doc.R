#' Spinal bone mineral density measurements (spnbmd)
#'
#' Relative spinal bone mineral density measurements on 261 North
#' American adolescents.  Each value is the difference in spnbmd
#' taken on two consecutive visits, divided by the average. The age is
#' the average age over the two visits.
#'
#' @format Variables:
#' \describe{
#'   \item{idnum}{identifies the child, and hence the repeat measurements}
#'   \item{age}{average age of child when measurements were taken}
#'   \item{gender}{male or female}
#'   \item{spnbmd}{Relative Spinal bone mineral density measurement}
#' }
#' @source \url{http://web.stanford.edu/~hastie/ElemStatLearn/}
#' @references 
#' Bachrach LK, Hastie T, Wang M-C, Narasimhan B, Marcus R. Bone Mineral 
#' Acquisition in Healthy Asian, Hispanic, Black and Caucasian Youth. A 
#' Longitudinal Study. J Clin Endocrinol Metab (1999) 84, 4702-12.
"spnbmd"


#' Biomechanics Data (Biomechanics_Data)
#'
#' Paired sample of n=11 subjects. 
#' 
#' Each subject performed one running stride with an extra cushioned 
#' running shoe and one running stride with a normally cushioned 
#' running shoe. 
#' 
#' Measured outcome: 
#' Torque (N m / kg) at ankle joint in sagittal pane during the stance phase
#' standardized with respect to the bodyweight of the subjects.
#'
#' @format Variables:
#' \describe{
#'   \item{Stance_Phase}{in percent}
#'   \item{Extra_Cushioned_1:Extra_Cushioned_11}{Torque curves (N m / kg) when running with extra cushioned running shoes.}
#'   \item{Normal_Cushioned_1:Normal_Cushioned_11}{Torque curves (N m / kg) when running with normal cushioned running shoes.}
#' }
"Biomechanics"