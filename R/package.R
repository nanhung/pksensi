#' About \pkg{pksensi} package
#'
#' Applying a global sensitivity analysis approach to reduce parameter dimensionality in
#' physiologically based kinetic modeling and evaluate the robustness of the algorithm under the given sampling number.
#'
#' The "extended Fourier amplitude sensitivity testing (eFAST)", a variance-based sensitivity analysis
#' method is used to estimate the parameter impact on model output
#' (Saltelli et al. 1999).
#' The eFAST is the effective algorithm to determine the influential parameter in
#' physiologically-based pharmacokinetic model calibration (Hsieh et al. 2018).
#' The eFAST algorithm is sourced from \pkg{sensitivity} package
#' but implemented the random-phase shift to evaluating the robustness of sensitivity measurement under the given sample size.
#'
#' @references
#' Saltelli, A., Tarantola, S., & Chan, K. S. (1999).
#' A quantitative model-independent method for global sensitivity analysis of model output.
#' \emph{Technometrics}, 41, 39-56.
#'
#' Hsieh, N. H., Reisfeld, B., Bois, F. Y., & Chiu, W. A. (2018).
#' Applying a global sensitivity analysis workflow to improve the computational efficiencies
#' in physiologically-based pharmacokinetic modeling.
#' \emph{Frontiers in Pharmacology}, 9, 588.
#'
#' @name about-pksensi
#' @aliases pksensi-package
NULL

