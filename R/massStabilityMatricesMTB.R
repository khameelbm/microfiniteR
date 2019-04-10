#' Mass matrix of a modified couple stress Timoshenko beam.
#'
#' Forms the mass matrix of a microscale Timoshenko beam (MTB)
#'
#' @param youngmod            Young's modulus.
#' @param shearmod            Shear modulus.
#' @param momentinertia       Moment of inertia (e.g. bh^3/12 for a rectangle).
#' @param area                Cross-sectional area
#' @param shearfactor         Shear correction factor (constant k).
#' @param poissonratio        Poisson's ratio.
#' @param length              Element's length.
#' @param rho                 Mass density (kg/cubic meter)
#'
#' @return                    Mass matrix involving both translational and
#'                            rotary mass matrices.
#' @export
#'
FormMassMTB <- function(youngmod, shearmod,
                        momentinertia, area,
                        shearfactor, poissonratio,
                        length, rho){

  L <- length
  ele_DOF <- 4
  phi <- (12 * youngmod * momentinertia) / (shearfactor * area* shearmod * L^2)
  beta3 <- (rho * area * length) / (420*(1 + phi)^2)
  beta4 <- (rho * momentinertia)/(30*L*(1 + phi)^2)

  m1 <- (156 + 294*phi + 140*phi^2)
  m2 <- (22 + 38.5*phi + 17.5*phi^2)
  m3 <- (54 + 126*phi + 70*phi^2)
  m4 <- (13 + 31.5*phi + 17.5*phi^2)
  m5 <- (4 + 7*phi + 3.5*phi^2)
  m6 <- (3 + 7*phi + 3.5*phi^2)
  m7 <- 37
  m8 <- 3 - 15*phi
  m9 <- 4 + 5*phi + 10*phi^2
  m10 <- 1 + 5*phi - 5*phi^2

  row1 <- c(m1, L*m2, m3, -L*m4)
  row2 <- c(L*m2, L^2*m5, L*m4, -L^2*m6)
  row3 <- c(m3, L*m4, m1, -L*m2)
  row4 <- c(-L*m4, -L^2*m6, -L*m2, L^2*m5)

  row21 <- c(m7, L*m8, -m7, L*m8)
  row22 <- c(L*m8, L^2*m9, -L*m8, -L^2*m10)
  row23 <- c(-m7, -L*m8, m7, -L*m8)
  row24 <- c(L*m8, -L^2*m10, -L*m8, L^2*m9)

  tran_massmatrix <- beta3 * matrix(c(row1,row2,row3,row4),
                                    nrow = ele_DOF,byrow = T)

  rota_massmatrix <- beta4 * matrix(c(row21,row22,row23,row24),
                                    nrow = ele_DOF,byrow = T)

  totalmass <- tran_massmatrix+rota_massmatrix
  return(totalmass)


}



#' Stability matrix of a modified couple stress Timoshenko beam.
#'
#' Forms the stability matrix of a microscale Timoshenko beam (MTB)
#'
#' @param youngmod            Young's modulus.
#' @param shearmod            Shear modulus.
#' @param momentinertia       Moment of inertia (e.g. bh^3/12 for a rectangle).
#' @param area                Cross-sectional area
#' @param shearfactor         Shear correction factor (constant k).
#' @param poissonratio        Poisson's ratio.
#' @param length              Element's length.
#'
#' @return                   Stability matrix for a microscale Timoshenko beam
#' @export
#'
FormStabilityMTB <- function(youngmod, shearmod,
                             momentinertia, area,
                             shearfactor, poissonratio, length){

  L <- length
  ele_DOF <- 4
  phi <- (12 * youngmod * momentinertia) / (shearfactor * area* shearmod * L^2)
  gamma <- 1/ (10*L*(1 + phi)^2)

  s1 <- 2*(6 + 10*phi + 5*phi^2)
  s2 <- L
  s3 <- (1/6)*(8 + 10*phi + 5*phi^2)
  s4 <- (1/6)*(2 + 10*phi + 5*phi^2)


  row1 <- c(s1, s2, -s1, s2)
  row2 <- c(s2, (L^2)*s3, -s2, -(L^2)*s4)
  row3 <- c(-s1, -s2, s1, -s2)
  row4 <- c(s2, -s4*(L^2), -s2, s3*(L^2))

  stabilitymatrix <- gamma * matrix(c(row1,row2,row3,row4),
                                    nrow = ele_DOF, byrow = T)
  return(stabilitymatrix)

}


#' Buckling loads (Timoshenko beam element)
#'
#' This function computes the buckling loads, but returns only the critical buckling load.
#'
#' @param reducedS      Reduced stability matrix obtained by applying
#'                      boundary conditions on the global stability matrix
#' @param reducedK      Reduced mass matrix obtained by applying
#'                      boundary conditions on the global stiffness matrix
#'
#' @return              Critical buckling load
#' @export

FindCriticalLoad <- function(reducedM, reducedK)
{
  stabilityinv <- solve(reducedS)
  productSK <- stabilityinv %*% reducedK
  criticalloads <- eigen(productSK)
  sortedloads <- sort(criticalloads$values)[1]
  return(sortedloads)
}


#' Natural frequencies (timoshenko beam element in radians)
#'
#' This function computes the natural frequencies, in radians per seconds.
#'
#' @param reducedM      Reduced mass matrix obtained by applying
#'                      boundary condition on the global mass matrix
#' @param reducedK      Reduced mass matrix obtained by applying
#'                      boundary condition on the global stiffness matrix
#'
#' @return              Natural frequencies in radian/s
#' @export
#'
FindFrequenciesRad <- function(reducedM, reducedK)
{
  massinv <- solve(reducedM);
  productMK <- massinv%*%reducedK
  syseigen <- eigen(productMK)
  sortedfreq <- sort(syseigen$values)

  return(sqrt(sortedfreq))
}

