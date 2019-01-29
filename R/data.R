#' Pharmacokinetics of Acetaminophen
#'
#' The APAP dataset contains a human experiment of pharmacokinetic data
#' for acetaminophen, and it's major metabolites.
#'
#' @references
#' T. J. Zurlinden and B. Reisfeld, 2016,
#' Physiologically based modeling of the pharmacokinetics of acetaminophen
#' and its major metabolites in humans using a Bayesian population approach,
#' \emph{European Journal of Drug Metabolism and Pharmacokinetics}, 79(4), 1-26.
#'
#' @format A data frame containing the following columns:
#' \itemize{
#'   \item Wt: averaged weight of the study population (kg).
#'   \item Dose: given oral dose of acetaminophen administered (mg/kg).
#'   \item Time: time after drug administration (hr).
#'   \item APAP: concentration of acetaminophen in the sample (mg/L).
#'   \item AG: concentration of acetaminophen-glucuronide in the sample (mg/L).
#'   \item AS: concentration of acetaminophen-sulfate in the sample (mg/L).
#'   \item Study: Sourced reference.
#' }
"APAP"
