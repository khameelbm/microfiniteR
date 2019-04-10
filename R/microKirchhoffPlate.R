
#' Stiffness matrix of a modified couple stress Kirchhoff plate element.
#'
#' Forms the stiffness matrix of a microscale Kirchhoff plate element (MKP)
#'
#' @param youngmod                Young's modulus
#' @param poissonratio            Poisson's ratio
#' @param edge_a                  Length along x-axis
#' @param edge_b                  Length along y-axis
#' @param thickness               Plate thinckness
#' @param lengthscale             Material length scale parameter
#'
#' @return                        Stiffness matrix (combination of the
#'                                classical stiffness matrix of a Timoshenko beam with
#'                                the effect arising from the modified
#'                                couple stress theory).
#' @export
FormStiffnessMKplate<- function(youngmod, poissonratio, edge_a,
                                edge_b, thickness,lengthscale){
    ele_DOF <- 12
    a <- edge_a
    b <- edge_b
    h <- thickness
    v <- poissonratio
    r <- a/b
    alpha1 <- youngmod * h^3 / (12* 30 * b^2 * r^3* (1 - v^2))
    alpha2 <- (youngmod * lengthscale^2 * h) / (30 * b^2 * r^3* (1 + v))

    alpha <- (alpha1 + alpha2)

  r1 <- c(30 + 21*r^2 + 30*r^4 - 6*r^2*v,
       3*b*r^2 + 30*b*r^4 + 12*b*r^2*v,
       -30*b*r - 3*b*r^3 - 12*b*r^3*v,
       -30 - 21*r^2 + 15*r^4 + 6*r^2*v,
       -3*b*r^2 + 15*b*r^4 - 12*b*r^2*v,
       -30*b*r - 3*b*r^3 + 3*b*r^3*v,
       -15 + 21*r^2 - 15*r^4 - 6*r^2*v,
       -3*b*r^2 + 15*b*r^4 + 3*b*r^2*v,
       -15*b*r + 3*b*r^3 - 3*b*r^3*v,
       15 - 21*r^2 - 30*r^4 + 6*r^2*v,
       3*b*r^2 + 30*b*r^4 - 3*b*r^2*v,
       -15*b*r + 3*b*r^3 + 12*b*r^3*v
       )
  r2 <- c(3*b*r^2 + 30*b*r^4 + 12*b*r^2*v,
         8*b^2*r^2 + 40*b^2*r^4 - 8*b^2*r^2*v,
         -30*b^2*r^3*v,
         -3*b*r^2 + 15*b*r^4 - 12*b*r^2*v,
         -8*b^2*r^2 + 20*b^2*r^4 + 8*b^2*r^2*v,
         0,
         3*b*r^2 - 15*b*r^4 - 3*b*r^2*v,
         2*b^2*r^2 + 10*b^2*r^4 - 2*b^2*r^2*v,
         0,
         -3*b*r^2 - 30*b*r^4 + 3*b*r^2*v,
         -2*b^2*r^2 + 20*b^2*r^4 + 2*b^2*r^2*v,
         0)

  r3 <- c(-30*b*r - 3*b*r^3 - 12*b*r^3*v,
          -30*b^2*r^3*v,
          40*b^2*r^2 + 8*b^2*r^4 - 8*b^2*r^4*v,
          30*b*r + 3*b*r^3 - 3*b*r^3*v,
          0,
          20*b^2*r^2 - 2*b^2*r^4 + 2*b^2*r^4*v,
          15*b*r - 3*b*r^3 + 3*b*r^3*v,
          0,
          10*b^2*r^2 + 2*b^2*r^4 - 2*b^2*r^4*v,
          -15*b*r + 3*b*r^3 + 12*b*r^3*v,
          0,
          20*b^2*r^2 - 8*b^2*r^4 + 8*b^2*r^4*v)

  r4 <- c(-30 - 21*r^2 + 15*r^4 + 6*r^2*v,
          -3*b*r^2 + 15*b*r^4 - 12*b*r^2*v,
          30*b*r + 3*b*r^3 - 3*b*r^3*v,
          30 + 21*r^2 + 30*r^4 - 6*r^2*v,
          3*b*r^2 + 30*b*r^4 + 12*b*r^2*v,
          30*b*r + 3*b*r^3 + 12*b*r^3*v,
          15 - 21*r^2 - 30*r^4 + 6*r^2*v,
          3*b*r^2 + 30*b*r^4 - 3*b*r^2*v,
          15*b*r - 3*b*r^3 - 12*b*r^3*v,
          -15 + 21*r^2 - 15*r^4 - 6*r^2*v,
          -3*b*r^2 + 15*b*r^4 + 3*b*r^2*v,
          15*b*r - 3*b*r^3 + 3*b*r^3*v)
  r5 <- c(-3*b*r**2 + 15*b*r**4 - 12*b*r**2*v,
          -8*b**2*r**2 + 20*b**2*r**4 + 8*b**2*r**2*v,
          0,
          3*b*r**2 + 30*b*r**4 + 12*b*r**2*v,
          8*b**2*r**2 + 40*b**2*r**4 - 8*b**2*r**2*v,
          30*b**2*r**3*v,
          -3*b*r**2 - 30*b*r**4 + 3*b*r**2*v,
          -2*b**2*r**2 + 20*b**2*r**4 + 2*b**2*r**2*v,
          0, 3*b*r**2 - 15*b*r**4 - 3*b*r**2*v,
          2*b**2*r**2 + 10*b**2*r**4 - 2*b**2*r**2*v,
          0)
  r6 <- c(-30*b*r - 3*b*r**3 + 3*b*r**3*v, 0,
          20*b**2*r**2 - 2*b**2*r**4 + 2*b**2*r**4*v,
          30*b*r + 3*b*r**3 + 12*b*r**3*v,
          30*b**2*r**3*v,
          40*b**2*r**2 + 8*b**2*r**4 - 8*b**2*r**4*v,
          15*b*r - 3*b*r**3 - 12*b*r**3*v,
          0,
          20*b**2*r**2 - 8*b**2*r**4 + 8*b**2*r**4*v,
          -15*b*r + 3*b*r**3 - 3*b*r**3*v, 0,
          10*b**2*r**2 + 2*b**2*r**4 - 2*b**2*r**4*v)

  r7 <- c(-15 + 21*r**2 - 15*r**4 - 6*r**2*v,
          3*b*r**2 - 15*b*r**4 - 3*b*r**2*v,
          15*b*r - 3*b*r**3 + 3*b*r**3*v,
          15 - 21*r**2 - 30*r**4 + 6*r**2*v,
          -3*b*r**2 - 30*b*r**4 + 3*b*r**2*v,
          15*b*r - 3*b*r**3 - 12*b*r**3*v,
          30 + 21*r**2 + 30*r**4 - 6*r**2*v,
          -3*b*r**2 - 30*b*r**4 - 12*b*r**2*v,
          30*b*r + 3*b*r**3 + 12*b*r**3*v,
          -30 - 21*r**2 + 15*r**4 + 6*r**2*v,
          3*b*r**2 - 15*b*r**4 + 12*b*r**2*v,
          30*b*r + 3*b*r**3 - 3*b*r**3*v)
  r8 <- c(-3*b*r**2 + 15*b*r**4 + 3*b*r**2*v,
          2*b**2*r**2 + 10*b**2*r**4 - 2*b**2*r**2*v,
          0, 3*b*r**2 + 30*b*r**4 - 3*b*r**2*v,
          -2*b**2*r**2 + 20*b**2*r**4 + 2*b**2*r**2*v,
          0, -3*b*r**2 - 30*b*r**4 - 12*b*r**2*v,
          8*b**2*r**2 + 40*b**2*r**4 - 8*b**2*r**2*v,
          -30*b**2*r**3*v,
          3*b*r**2 - 15*b*r**4 + 12*b*r**2*v,
          -8*b**2*r**2 + 20*b**2*r**4 + 8*b**2*r**2*v,
          0)

  r9 <- c(-15*b*r + 3*b*r**3 - 3*b*r**3*v,
          0, 10*b**2*r**2 + 2*b**2*r**4 - 2*b**2*r**4*v,
          15*b*r - 3*b*r**3 - 12*b*r**3*v, 0,
          20*b**2*r**2 - 8*b**2*r**4 + 8*b**2*r**4*v,
          30*b*r + 3*b*r**3 + 12*b*r**3*v,
          -30*b**2*r**3*v,
          40*b**2*r**2 + 8*b**2*r**4 - 8*b**2*r**4*v,
          -30*b*r - 3*b*r**3 + 3*b*r**3*v, 0,
          20*b**2*r**2 - 2*b**2*r**4 + 2*b**2*r**4*v)

  r10 <- c(15 - 21*r**2 - 30*r**4 + 6*r**2*v,
           -3*b*r**2 - 30*b*r**4 + 3*b*r**2*v,
           -15*b*r + 3*b*r**3 + 12*b*r**3*v,
           -15 + 21*r**2 - 15*r**4 - 6*r**2*v,
           3*b*r**2 - 15*b*r**4 - 3*b*r**2*v,
           -15*b*r + 3*b*r**3 - 3*b*r**3*v,
           -30 - 21*r**2 + 15*r**4 + 6*r**2*v,
           3*b*r**2 - 15*b*r**4 + 12*b*r**2*v,
           -30*b*r - 3*b*r**3 + 3*b*r**3*v,
           30 + 21*r**2 + 30*r**4 - 6*r**2*v,
           -3*b*r**2 - 30*b*r**4 - 12*b*r**2*v,
           -30*b*r - 3*b*r**3 - 12*b*r**3*v)

  r11 <- c(3*b*r**2 + 30*b*r**4 - 3*b*r**2*v,
           -2*b**2*r**2 + 20*b**2*r**4 + 2*b**2*r**2*v,
           0, -3*b*r**2 + 15*b*r**4 + 3*b*r**2*v,
           2*b**2*r**2 + 10*b**2*r**4 - 2*b**2*r**2*v,
           0, 3*b*r**2 - 15*b*r**4 + 12*b*r**2*v,
           -8*b**2*r**2 + 20*b**2*r**4 + 8*b**2*r**2*v,
           0, -3*b*r**2 - 30*b*r**4 - 12*b*r**2*v,
           8*b**2*r**2 + 40*b**2*r**4 - 8*b**2*r**2*v,
           30*b**2*r**3*v)

  r12 <- c(-15*b*r + 3*b*r**3 + 12*b*r**3*v,
           0,
           20*b**2*r**2 - 8*b**2*r**4 + 8*b**2*r**4*v,
           15*b*r - 3*b*r**3 + 3*b*r**3*v,
           0,
           10*b**2*r**2 + 2*b**2*r**4 - 2*b**2*r**4*v,
           30*b*r + 3*b*r**3 - 3*b*r**3*v,
           0, 20*b**2*r**2 - 2*b**2*r**4 + 2*b**2*r**4*v,
           -30*b*r - 3*b*r**3 - 12*b*r**3*v,
           30*b**2*r**3*v,
           40*b**2*r**2 + 8*b**2*r**4 - 8*b**2*r**4*v)

  stiffnessmatrixc <- alpha1 * matrix(c(r1, r2, r3, r4,
                              r5, r6, r7, r8, r9,
                              r10, r11, r12), byrow = T, nrow = ele_DOF)

  r11 <- c(30 - 15*r**2 + 30*r**4,-15*b*r**2 + 30*b*r**4,-30*b*r + 15*b*r**3,-30
           + 15*r**2 + 15*r**4,15*b*r**2 + 15*b*r**4,-30*b*r,-15 - 15*r**2 -
           15*r**4,15*b*r**4,-15*b*r,15 + 15*r**2 - 30*r**4,30*b*r**4,-15*b*r -
           15*b*r**3)
  r12 <- c(-15*b*r**2 + 30*b*r**4,40*b**2*r**4,30*b**2*r**3,15*b*r**2 +
            15*b*r**4,20*b**2*r**4,0,-15*b*r**4,10*b**2*r**4,0,-30*b*r**4,20*b**2*
            r**4,0)
  r13 <- c(-30*b*r +
            15*b*r**3,30*b**2*r**3,40*b**2*r**2,30*b*r,0,20*b**2*r**2,15*b*r,0,10*
            b**2*r**2,-15*b*r - 15*b*r**3,0,20*b**2*r**2)
  r14 <- c(-30 + 15*r**2 + 15*r**4,15*b*r**2 + 15*b*r**4,30*b*r,30 - 15*r**2 +
           30*r**4,-15*b*r**2 + 30*b*r**4,30*b*r - 15*b*r**3,15 + 15*r**2 -
           30*r**4,30*b*r**4,15*b*r + 15*b*r**3,-15 - 15*r**2 -
           15*r**4,15*b*r**4,15*b*r)
  r15 <- c(15*b*r**2 + 15*b*r**4,20*b**2*r**4,0,-15*b*r**2 +
           30*b*r**4,40*b**2*r**4,-30*b**2*r**3,-30*b*r**4,20*b**2*r**4,0,-15*b*
           r**4,10*b**2*r**4,0)
  r16 <- c(-30*b*r,0,20*b**2*r**2,30*b*r -
           15*b*r**3,-30*b**2*r**3,40*b**2*r**2,15*b*r +
           15*b*r**3,0,20*b**2*r**2,-15*b*r,0,10*b**2*r**2)
  r17 <- c(-15 - 15*r**2 - 15*r**4,-15*b*r**4,15*b*r,15 + 15*r**2 -
           30*r**4,-30*b*r**4,15*b*r + 15*b*r**3,30 - 15*r**2 +
           30*r**4,15*b*r**2 - 30*b*r**4,30*b*r - 15*b*r**3,-30 + 15*r**2 +
           15*r**4,-15*b*r**2 - 15*b*r**4,30*b*r)
  r18 <- c(15*b*r**4,10*b**2*r**4,0,30*b*r**4,20*b**2*r**4,0,15*b*r**2 -
           30*b*r**4,40*b**2*r**4,30*b**2*r**3,-15*b*r**2 -
           15*b*r**4,20*b**2*r**4,0)
  r19 <- c(-15*b*r,0,10*b**2*r**2,15*b*r + 15*b*r**3,0,20*b**2*r**2,30*b*r -
           15*b*r**3,30*b**2*r**3,40*b**2*r**2,-30*b*r,0,20*b**2*r**2)
  r20 <- c(15 + 15*r**2 - 30*r**4,-30*b*r**4,-15*b*r - 15*b*r**3,-15 - 15*r**2 -
           15*r**4,-15*b*r**4,-15*b*r,-30 + 15*r**2 + 15*r**4,-15*b*r**2 -
           15*b*r**4,-30*b*r,30 - 15*r**2 + 30*r**4,15*b*r**2 -
           30*b*r**4,-30*b*r + 15*b*r**3)
  r21 <- c(30*b*r**4,20*b**2*r**4,0,15*b*r**4,10*b**2*r**4,0,-15*b*r**2 -
           15*b*r**4,20*b**2*r**4,0,15*b*r**2 -
           30*b*r**4,40*b**2*r**4,-30*b**2*r**3)
  r22 <- c(-15*b*r -
           15*b*r**3,0,20*b**2*r**2,15*b*r,0,10*b**2*r**2,30*b*r,0,20*b**2*r**2,-
           30*b*r + 15*b*r**3,-30*b**2*r**3,40*b**2*r**2)
  stiffnessmatrixmcst <- alpha2 * matrix(c(r11, r12, r13, r14,
                                        r15, r16, r17, r18, r19,
                                        r20, r21, r22), byrow = T, nrow = ele_DOF)

  micro_stiffnessmatrix <-stiffnessmatrixc + stiffnessmatrixmcst
  return(micro_stiffnessmatrix)

  }


#' Mass matrix of a modified couple stress Kirchhoff plate element.
#'
#' Forms the mass matrix of a microscale Kirchhoff plate element (MKP)
#'
#' @param youngmod                Young's modulus
#' @param edge_a                  Length along x-axis
#' @param edge_b                  Length along y-axis
#' @param thickness               Plate thinckness
#' @param rho                     Mass density of the plate
#'
#' @return                       Mass matrix (combination of the
#'                               classical stiffness matrix of a Timoshenko beam with
#'                               the effect arising from the modified
#'                               couple stress theory).
#' @export

FormMassMKplate<- function(youngmod, edge_a, edge_b, thickness, rho){

  ele_DOF <- 12
  a <- edge_a
  b <- edge_b
  r <- a/b
  h <- thickness
  beta <- (rho * h * a * b)/6300

  r1 <- c(3454,922*b,-922*a,1226,398*b,548*a,394,-232*b,232*a,1226,-548*b,
          -398*a)
  r2 <- c(922*b,320*b**2,-252*a*b,398*b,160*b**2,168*a*b,232*b,-120*b**2,
          112*a*b,548*b,-240*b**2,-168*a*b)
  r3 <- c(-922*a,-252*a*b,320*a**2,-548*a,-168*a*b,-240*a**2,-232*a,112*a*b,
          -120*a**2,-398*a,168*a*b,160*a**2)
  r4 <- c(1226,398*b,-548*a,3454,922*b,922*a,1226,-548*b,398*a,394,-232*b,-232*a)
  r5 <- c(398*b,160*b**2,-168*a*b,922*b,320*b**2,252*a*b,548*b,-240*b**2,168*a*b,
          232*b,-120*b**2,-112*a*b)
  r6 <- c(548*a,168*a*b,-240*a**2,922*a,252*a*b,320*a**2,398*a,-168*a*b,
          160*a**2,232*a,-112*a*b,-120*a**2)
  r7 <- c(394,232*b,-232*a,1226,548*b,398*a,3454,-922*b,922*a,1226,-398*b,-548*a)
  r8 <- c(-232*b,-120*b**2,112*a*b,-548*b,-240*b**2,-168*a*b,-922*b,
          320*b**2,-252*a*b,-398*b,160*b**2,168*a*b)
  r9 <- c(232*a,112*a*b,-120*a**2,398*a,168*a*b,160*a**2,922*a,-252*a*b,
          320*a**2,548*a,-168*a*b,-240*a**2)
  r10 <- c(1226,548*b,-398*a,394,232*b,232*a,1226,-398*b,548*a,3454,-922*b,-922*a)
  r11 <- c(-548*b,-240*b**2,168*a*b,-232*b,-120*b**2,-112*a*b,-398*b,160*b**2,
           -168*a*b,-922*b,320*b**2,252*a*b)
  r12 <- c(-398*a,-168*a*b,160*a**2,-232*a,-112*a*b,-120*a**2,-548*a,168*a*b,
           -240*a**2,-922*a,252*a*b,320*a**2)

  massmatrixp <- beta * matrix(c(r1, r2, r3, r4,
                                 r5, r6, r7, r8, r9,
                                 r10, r11, r12), byrow = T, nrow = ele_DOF)
  return(massmatrixp)

}




#' Expanded stiffness matrix of a microscale Kirchhoff plate
#'
#' Accepts the element's stiffness matrix and returns the expanded form
#' of the matrix to conform with the global degree of freedom.
#'
#' @param tdof          Total degree of freedom in a
#'                      connected system of plates.
#'
#' @param eMatrix       Element's stiffness matrix.
#' @param i             Index of the first node.
#' @param j             Index of the second node.
#' @param k             Index of the third node.
#' @param m             Index of the fourth node.
#'
#' @return              The expanded stiffness matrix
#'                      a microscale Kirchhoff plate element.
#' @export

ExpandStiffnessMKPlate <-  function(tdof, eMatrix, i, j, k, m) {

  r1 = 3 * i - 2; r2 = 3 * i - 1; r3 = 3 * i;
  r4 = 3 * j - 2; r5 = 3 * j - 1; r6 = 3 * j;
  r7 = 3 * k - 2; r8 = 3 * k - 1; r9 = 3 * k;
  r10 = 3 * m - 2; r11 = 3 * m - 1; r12 = 3 * m;

  vec <- c(r1, r2, r3, r4, r5, r6,
           r7, r8, r9, r10, r11, r12)

  mtbbigMatrix <- matrix(vector(l = tdof * tdof), nrow = tdof, byrow = T)
  mtbbigMatrix[vec, vec] <- eMatrix

  return (mtbbigMatrix)

}


#' Equivalent load vector (non-conforming rectangular Kirchhoff plate element)
#'
#' Generates the column matrix of equivalent loads for a
#'           non-conforming rectangular Kirchhoff plate.
#'
#' @param pressuremag             Pressure magnitude, e.g. q in N/m^2.
#' @param edge_a                  Length along x-axis
#' @param edge_b                  Length along y-axis
#'
#' @return                        A column matrix of the equivalent nodal loads.
#' @export
#'
FormLoadMKPlate <-  function(pressuremag, edge_a, edge_b){
  ele_DOF <- 12
  a <- edge_a
  b <- edge_b
  beta <- (pressuremag * a * b)/3
  contents <- c(3, b, -a, 3, b, a, 3, -b, a, 3, -b, -a)
  equivalentload <- beta * matrix(contents, nrow = ele_DOF, byrow = T)
  return (equivalentload)
}




#' Expanded vector of the equivalent load
#'
#' This function generates the expanded vector of the equivalent load
#'                                             (Euler-Bernoulli beam).
#'
#' @param tdof                       Total degree of freedom.
#' @param elementloadMatrix          The unexpanded vector of equivalent loads.
#' @param i                          Index of the first node.
#' @param j                          Index of the second node.
#' @param k                          Index of the third node.
#' @param m                          Index of the fourth node.
#'
#' @return                           Expanded vector (a column matrix) of equivalent loads.
#' @export

ExpandLoadMKPlate <- function(tdof, elementloadMatrix, i, j, k, m){
  r1 <-  3 * i - 2; r2 <-  3 * i - 1; r3 <-  3 * i;
  r4 <-  3 * j - 2; r5 <-  3 * j - 1; r6 <- 3 * j;
  r7 <-  3 * k - 2; r8 <-  3 * k - 1; r9 <-  3 * k;
  r10 <-  3 * m - 2; r11 <-  3 * m - 1; r12 <-  3 * m;

  contents <- c(r1, r2, r3, r4, r5, r6,r7, r8, r9, r10, r11, r12)

  bigLoadMKPlate <- matrix(0, nrow = tdof, byrow = T);
  bigLoadMKPlate[contents] <- elementloadMatrix;

  return (bigLoadMKPlate)
}


#' Reduced load vector.
#'
#' Forms the reduced load matrix after applying the boundary condition(s).
#'
#' @param bigColumnMatrix          A vector of applied nodal loads (force,moments etc)
#'                                 obtained from using FormLoadMKPlate() for a single element
#'                                 or the assembled load matrix
#' @param unrestrainednodes        A set of freenodes retrieved after applying the boundary condition.
#'
#' @return                        A column matrix of the reduced load vector
#' @export
#'
FormReducedLoad2D <- function(bigColumnMatrix, unrestrainednodes){

  reducedloadvector <- bigColumnMatrix[unrestrainednodes]
  return(cbind(reducedloadvector))
}

