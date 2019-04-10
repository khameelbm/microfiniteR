
#' Fixed or clamped boundary condition(s)
#'
#' Retrieves the translational and rotational degrees of freedom
#'                                  of a node or a set of node.
#'
#' @param nodes        A vector of node(s) for which the the
#'                     translational and rotational displacements are to be suppressed.
#'
#' @return             A vector of rows (fixed nodes) to be eliminated
#'                     in the global matrix equation.
#' @export

FixNodes2D=function(nodes){

  len <- length(nodes)
  fixed_dof <- vector(mode = "numeric")

  if (length(nodes) == 1){
    fixed_dof[nodes] <- 3 * nodes
    fixed_dof[nodes+1] <- 3 * nodes - 2
    fixed_dof[nodes+2] <- 3 * nodes - 1
    return(sort(fixed_dof))

  } else {

    hinged_dof <- vector(mode = "numeric") #displacement
    fixed1_dof <- vector(mode = "numeric") #first slope
    fixed2_dof <- vector(mode = "numeric") #second slope

    for (k in (1:len)){hinged_dof[k] <- 3 * nodes[k]}
    for (k in (1:len)){fixed1_dof[k] <- 3 * nodes[k] - 2}
    for (k in (1:len)){fixed2_dof[k] <- 3 * nodes[k] - 1}
  }

  return(sort(c(hinged_dof, fixed1_dof, fixed2_dof)))
}


#' Simply supported boundary condition(s) along x
#'
#' Retrieves the displacement and theta_y of a node or a set of node along x.
#'
#'
#' @param nodes        A vector of node(s) for which the
#'                     translational displacement is to be suppressed.
#'
#' @return             A vector of rows (hinged nodes) to be eliminated
#'                     in the global matrix equation.
#' @export

HingeNodes2Dx <- function(nodes){

  len <- length(nodes)
  disp_dof <- vector(mode = "numeric")

  if (length(nodes) == 1){
    disp_dof[nodes] <- 3 * nodes
    disp_dof[nodes+1] <- 3 * nodes - 2
    return(sort(disp_dof))

  } else {

    disp_dof <- vector(mode = "numeric")          #displacements
    rotationy_dof <- vector(mode = "numeric")     #second slope or thetay

    for (k in (1:len)){disp_dof[k] <- 3 * nodes[k] - 2}
    for (k in (1:len)){rotationy_dof[k] <- 3 * nodes[k]}
  }

  return(sort(c(disp_dof, rotationy_dof)))
}



#' Simply supported boundary condition(s) along y
#'
#' Retrieves the displacement and theta_x of a node or a set of node along y.
#'
#'
#' @param nodes        A vector of node(s) for which the
#'                     translational displacement is to be suppressed.
#'
#' @return             A vector of rows (hinged nodes) to be eliminated
#'                     in the global matrix equation.
#' @export
HingeNodes2Dy <- function(nodes){

  len <- length(nodes)
  disp_dof <- vector(mode = "numeric")

  if (length(nodes) == 1){
    disp_dof[nodes] <- 3 * nodes -1
    disp_dof[nodes+1] <- 3 * nodes - 2
    return(sort(disp_dof))

  } else {

    disp_dof <- vector(mode = "numeric")          #displacements
    rotationx_dof <- vector(mode = "numeric")     #second slope or thetax

    for (k in (1:len)){disp_dof[k] <- 3 * nodes[k] - 2}
    for (k in (1:len)){rotationx_dof[k] <- 3 * nodes[k] - 1}
  }

  return(sort(c(disp_dof, rotationx_dof)))
}


#' Symmetric boundary condition(s) along y
#'
#' Retrieves theta_y of a node or a set of node along y.
#'
#'
#' @param nodes        A vector of nodes for which the symmetric boundary condition(s) is to be applied
#'
#' @return             A vector of rows (hinged nodes) to be eliminated
#'                                        in the global matrix equation.
#' @export
Symmetry2Dy <- function(nodes){

  len <- length(nodes)

  if (length(nodes) == 1){

    rotationy_dof <- 3 * nodes
    return(rotationy_dof)

  } else {

    rotationy_dof <- vector(mode = "numeric")     #second slope or thetay
    for (k in (1:len)){rotationy_dof[k] <- 3 * nodes[k]}
    return(sort(rotationy_dof))
  }


}


#' Symmetric boundary condition(s) along x
#'
#' Retrieves theta_y of a node or a set of node along x.
#'
#'
#' @param nodes        A vector of nodes for which the symmetric boundary condition(s) is to be applied
#'
#' @return             A vector of rows (hinged nodes) to be eliminated
#'                                        in the global matrix equation.
#' @export
Symmetry2Dx <- function(nodes){

  len <- length(nodes)

  if (length(nodes) == 1){

    rotationx_dof <- 3 * nodes - 1
    return(rotationx_dof)

  } else {

    rotationx_dof <- vector(mode = "numeric")     #second slope or thetay
    for (k in (1:len)){rotationx_dof[k] <- 3 * nodes[k] - 1}
    return(sort(rotationx_dof))
  }


}


