
#' Stiffness matrix of a modified couple stress Timoshenko beam.
#'
#'Forms the stiffness matrix of a microscale Timoshenko beam (MTB)
#'
#' @param youngmod            Young's modulus.
#' @param shearmod            Shear modulus.
#' @param momentinertia       Moment of inertia (e.g. bh^3/12 for a rectangle).
#' @param area                Cross-section area.
#' @param shearfactor         Shear correction factor (constant k).
#' @param poissonratio        Poisson's ratio.
#' @param length              Element's length.
#' @param lengthscale         Material length scale parameter
#'
#' @return                    Stiffness matrix (combination of the
#'                            classical stiffness matrix of a Timoshenko beam with
#'                            the effect arising from the modified
#'                             couple stress theory).
#' @export

FormStiffnessMTB<- function(youngmod, shearmod,
                            momentinertia, area,
                            shearfactor, poissonratio,
                            length, lengthscale){

  L <- length
  ele_DOF <- 4
  phi <- (12 * youngmod * momentinertia) / (shearfactor * area* shearmod * L^2)
  beta1 <- (area * shearmod * lengthscale^2) / ((1 + phi)^2 * L^3)
  beta2 <- (youngmod * momentinertia)/((1 + phi) * L^3)
  bendshear <- beta2 * matrix(c(12, 6*L, -12, 6*L,
                                    6*L, (4+phi)*L^2, -6*L, (2-phi)*L^2,
                                    -12, -6*L, 12, -6*L,
                                    6*L, (2-phi)*L^2, -6*L,(4+phi)*L^2),
                                    nrow = ele_DOF, byrow = T)

  couplestress <- beta1 * matrix(c(12, 6*L, -12, 6*L,
                                      6*L, (4 + 2*phi + phi^2)*L^2,
                                      -6*L, -(-2 + 2*phi + phi^2)*L^2,
                                      -12, -6*L, 12, -6*L,
                                      6*L, -(-2 + 2*phi + phi^2)*L^2,
                                      -6*L, (4 + 2*phi + phi^2)*L^2),
                                      nrow = ele_DOF, byrow = T)
  totalstiffness <- bendshear + couplestress
  return(totalstiffness)
  }


#' Expanded stiffness matrix of a microscale Timoshenko beam
#'
#' Accepts the element's stiffness matrix and returns the expanded form
#' of the matrix to conform with the global degree of freedom.
#'
#' @param total_DOF    Total degree of freedom in a
#'                      connected system of beams.
#'
#' @param eMatrix       Element's stiffness matrix.
#' @param i             First node of an element.
#' @param j             Second node of an element.
#'
#' @return              The expanded stiffness matrix
#'                      a microscale Timoshenko beam element.
#' @export

ExpandStiffnessMTB <-  function(total_DOF,eMatrix,i,j) {

  r1 <- (i-1)+i
  r2 <- (i-1)+(i+1)
  r3 <- (j-2)+(j+1)
  r4 <- (j-2)+(j+2)

  mtbbigMatrix <- matrix(vector(l=total_DOF*total_DOF),nrow=total_DOF,
                         byrow=T)
  mtbbigMatrix[c(r1,r2,r3,r4),c(r1,r2,r3,r4)] <- eMatrix

  return (mtbbigMatrix)

  }


#' Reduced stiffness matrix
#'
#' @param globalK          Global stiffness matrix
#'                         (assembled matrices of all elements)
#'
#' @param knownforcenodes  The set of nodes for the applied loads
#'                          obtained by excluding the rows
#'                          of the fixed/hinged nodes in the global dof.
#'
#' @return                 Reduced stiffness matrix
#' @export
#'
FormReducedStiffness <- function(globalK,knownforcenodes){

  reducedstiff <- globalK[c(knownforcenodes),(knownforcenodes)]
  return(reducedstiff)
}


#' Reduced stiffness or mass matrix.
#'
#' Given the global stiffness or mass matrix, and
#' the vector of free nodes, this function returns the
#' reduced stiffness or reduced mass matrix.
#'
#' @param globalmat          Global mass matrix (assembled mass matrices)
#'
#' @param vec_freenodes      The set of nodes of freed nodes obtained by excluding the rows
#'                           of the fixed/hinged nodes in the global dof.
#'
#' @return                   Reduced stiffness or reduced mass matrix
#'
#' @export
#'
FormReducedMatrix <- function(globalmat,vec_freenodes){

  reducedmatrix <- globalmat[c(vec_freenodes),(vec_freenodes)]
  return(reducedmatrix)
}


#' Reduced force vector.
#'
#' Forms the reduced stiffness matrix after
#'          applying the boundary condition(s).
#'
#' @param forcevector  A vector of applied nodal loads (force,moments etc)
#'                     supplied in the form c(F1,M1, etc)
#'
#' @return             A column matrix of the reduced force vector
#' @export
FormReducedForce <- function(forcevector){

  reducedforcevector <- matrix(forcevector,ncol=1)
  return(reducedforcevector)
}


#' Nodal degrees of freedom
#'
#' Determines the unknown nodal degrees of freedom.
#'
#' @param reducedk         Reduced stiffness matrix
#'                         obtained with FormReducedStiffnes().
#'
#' @param reducedforce     Reduced force vector obtained with
#'                         FormReducedForce().
#'
#' @return                 A vector of all unknown nodal
#'                         degrees of freedom.
#' @export

FindNodalDOFs <- function(reducedk,reducedforce){

  nodal_dof <- solve(reducedk,reducedforce)
  return(nodal_dof)
}

#' Local elements forces and moments of a microscale Timoshenko beam.
#'
#' Solves for the nodal forces and moments of an element.
#'
#' @param youngmod           Young's modulus.
#' @param shearmod           Shear modulus.
#' @param momentinertia      Moment of inertia (e.g. bh^3/12 for rectangle).
#' @param area               Cross-section area.
#' @param shearfactor        Shear correction factor (constant k).
#' @param poissonratio       Poisson's ratio.
#' @param length             Element's length.
#' @param lengthscale        Material length scale parameter
#' @param global_dof         Global degrees of freedom
#'                           (combination of the dof found
#'                           with FindNodalDOF() and the boundary conditions).
#'
#' @param i                  Index of the first node.
#' @param j                  Index of the second node.
#'
#' @return                   Local forces and moments (MTB).
#' @export

FindForcesMomentsMTB <- function(youngmod, shearmod,
                                 momentinertia, area,
                                 shearfactor, poissonratio,
                                 length,lengthscale, global_dof,
                                 i,j){

    ele_DOF <- 4
    r1  <-  (i - 1) + i;
    r2  <-  (i - 1) + (i + 1);
    r3  <-  (j - 2) + (j + 1);
    r4  <-  (j - 2) + (j + 2);
    nodaldisp <- global_dof[c(r1,r2,r3,r4)];
    L <- length
    phi <- (12 * youngmod * momentinertia) / (shearfactor * area* shearmod * L^2)
    beta1 <- (area * shearmod * lengthscale^2) / ((1 + phi)^2 * L^3)
    beta2 <- (youngmod * momentinertia)/((1 + phi) * L^3)

    bendshear <- beta2 * matrix(c(12, 6*L, -12, 6*L,
                                  6*L, (4+phi)*L^2, -6*L,
                                  (2-phi)*L^2,-12, -6*L,
                                  12, -6*L,6*L,
                                  (2-phi)*L^2, -6*L,(4+phi)*L^2),
                                  nrow <- ele_DOF,byrow=T)

    couplestress <- beta1 * matrix(c(12, 6*L, -12, 6*L,
                                     6*L, (4 + 2*phi + phi^2)*L^2,
                                     -6*L, -(-2 + 2*phi + phi^2)*L^2,
                                     -12, -6*L, 12, -6*L,
                                     6*L, -(-2 + 2*phi + phi^2)*L^2,
                                     -6*L, (4 + 2*phi + phi^2)*L^2),
                                          nrow=ele_DOF,byrow=T)
    totalstiffness <- bendshear + couplestress
    local_forcesmoments <- totalstiffness %*% matrix(nodaldisp,ncol=1)
    return(local_forcesmoments)
  }


# Boundary condition functions -----------------------------------------

#' Hinged or simply supported boundary condition(s)
#'
#' Retrieves the translational displacement/dof of a node or a set of node.
#'
#'
#' @param nodes        A vector of node(s) for which the
#'                     translational displacement is to be suppressed.
#'
#' @return             A vector of rows (hinged nodes) to be eliminated
#'                     in the global matrix equation.
#' @export

HingeNodes <- function(nodes){

  len=length(nodes)
  hinged_dof=vector(mode="numeric")
  for (k in 1:len){
    hinged_dof[k] <- (nodes[k]-1)+nodes[k]

  }

  return(hinged_dof)
}


#' Fixed or clamped boundary condition(s)
#'
#' Retrieves the translational and rotational degrees of freedom
#'                                  of a node or a set of node.
#'
#' @param nodes        A vector of node(s) for which the
#'                     displacements are to be suppressed.
#'
#' @return             A vector of rows (fixed nodes) to be eliminated
#'                     in the global matrix equation.
#' @export

FixNodes=function(nodes){
  len <- length(nodes)
  fixed_dof <- vector(mode="numeric")
  if (length(nodes)==1){
    fixed_dof[nodes] <- (nodes-1)+nodes
    fixed_dof[nodes+1] <- ((nodes-1)+(nodes+1))
    return(fixed_dof)
  } else {

    hinged_dof <- vector(mode = "numeric")
    fixed_dof <- vector(mode = "numeric")
    for (k in 1:len){
      hinged_dof[k] <- (nodes[k]-1)+nodes[k]
    }

    for (k in (1:len)){
      fixed_dof[k] <- (nodes[k]-1)+(nodes[k]+1)
    }
  }

  return(sort(c(hinged_dof,fixed_dof)))
}


#' Rows of unrestrained nodes.
#'
#' Extracts the rows of unrestrained nodes.
#' This is required for reducing the global stiffness matrix.
#'
#'
#' @param total_DOF          A single number indicating the
#'                           total degree of freedom in a
#'                           connected system of beams.
#'
#' @param restrainednodes    A vector of restrained nodes.
#'
#' @return                   A vector of rows of known loads.
#'                           (or a vector of rows of unrestrained nodes).
#' @export
#'

ExtractFreeRows <- function(total_DOF,restrainednodes){
  vec=1:(total_DOF)
  nodes_knownloads=vec[-restrainednodes]
  return(nodes_knownloads)

}


