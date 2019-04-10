
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
#' Retrieves the translational and rotational degree of freedom
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

    hinged_dof <- vector(mode="numeric")
    fixed_dof <- vector(mode="numeric")
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
