#' Spinal bone mineral density measurements (Fragments)
#'
#' Relative spinal bone mineral density (spnbmd) measurements on 261 North
#' American adolescents.  Each value is the difference in spnbmd
#' taken on two consecutive visits, divided by the average. The age is
#' the average age over the two visits.
#'
#' @format Variables:
#' \describe{
#'   \item{Y_f.i}{spnbmd measurements of female subjects i=1,...,140}
#'   \item{X_f.i}{ages of the female subjects at time of measurements}
#'   \item{Y_m.i}{spnbmd measurements of male subjects i=1,...,113}
#'   \item{X_m.i}{ages of the male subjects at time of measurements}
#' }
#' @source \url{http://web.stanford.edu/~hastie/ElemStatLearn/}
#' @references 
#' Bachrach LK, Hastie T, Wang M-C, Narasimhan B, Marcus R. Bone Mineral 
#' Acquisition in Healthy Asian, Hispanic, Black and Caucasian Youth. A 
#' Longitudinal Study. J Clin Endocrinol Metab (1999) 84, 4702-12.
"Fragments"


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