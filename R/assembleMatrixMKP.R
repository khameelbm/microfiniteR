
#' Form global matrix from element matrix (stiffness, or mass) for a plate element
#'
#' @param enumber           Number of discretized plate elements.
#' @param listnodes         A list of nodal indices for the elements
#' @param youngmod          Young's modulus.
#' @param poissonratio      Poisson's ratio.
#' @param edge_a            Length along x-axis
#' @param edge_b            Length along y-axis
#' @param thickness         Plate thinckness
#' @param lengthscale       Material length scale parameter
#' @param rho               Mass density (kg/cubic meter)
#' @param case              Used to specify if stiffness matrix or mass matrix
#'                          or stability matrix. Note that case = 1 (stiffness),
#'                          or 2 (mass), or 3(stability))
#'
#'
#' @return                 Depending on what is being evaluated, it returns
#'                         the global stiffness matrix (case 1) or mass matrix (case 2)
#' @export
#'

FormGlobalKM2D=function(enumber, listnodes, youngmod, poissonratio, edge_a,
                       edge_b, thickness,lengthscale, rho, case){
  dof <- 12;
  numNodes = (2*enumber) + 1;
  dofpernode <- 3
  tdof <- numNodes * dofpernode;
  klist <- list();
  mlist <- list();


  if(case == 1){

    for(j in seq(enumber)){

      klist[[j]]=FormStiffnessMKplate(youngmod, poissonratio, edge_a,
                                      edge_b, thickness,lengthscale);
    }


    globalmatrix <- matrix(0, tdof, tdof);
    Klist <- list();

    for(j in seq(enumber)){
      nodalindex <- unlist(listnodes[j])
      Klist[[j]] <- ExpandStiffnessMKPlate(tdof, klist[[j]], nodalindex[1], nodalindex[2],
                                           nodalindex[3], nodalindex[4])
      globalmatrix <- globalmatrix + Klist[[j]]
    }


  }

  # Mass matrix
  if(case == 2){

    for(j in seq(enumber)){

      mlist[[j]]=FormMassMKplate(youngmod, edge_a, edge_b, thickness, rho);
    }

    globalmatrix <- matrix(0,tdof,tdof);
    Mlist <- list();

    for(j in seq(enumber)){
      nodalindex <- unlist(listnodes[j])
      Mlist[[j]] <- ExpandStiffnessMKPlate(tdof, mlist[[j]], nodalindex[1], nodalindex[2],
                                           nodalindex[3], nodalindex[4])
      globalmatrix <- globalmatrix + Mlist[[j]]
    }
  }


  return(globalmatrix)


}



#' Global load vector for a plate with surface UDL
#'
#' Form global load vector (a column matrix) for a plate with surface UDL
#'
#' @param enumber           Number of discretized plate elements.
#' @param listnodes         A list of nodal indices for the elements
#' @param edge_a            Length along x-axis
#' @param edge_b            Length along y-axis

#' @return                 Returns the global vector of equivalent load
#' @export
#'

FormGlobalLoad2D=function(enumber, listnodes,loadq, edgea, edgeb){
  dof <- 12;
  numNodes = (2*enumber) + 1;
  dofpernode <- 3
  tdof <- numNodes * dofpernode;
  loadlist <- list();
  case = 1
  if(case == 1){

    for(j in seq(enumber)){

      loadlist[[j]]=FormLoadMKPlate(loadq, edgea, edgeb);
    }

    globalloadmatrix <- matrix(0, tdof, byrow = T);
    Loadlist <- list();

    for(j in seq(enumber)){
      nodalindex <- unlist(listnodes[j])
      Loadlist[[j]] <- ExpandLoadMKPlate(tdof, loadlist[[j]], nodalindex[1], nodalindex[2],
                                           nodalindex[3], nodalindex[4])
      globalloadmatrix <- globalloadmatrix + Loadlist[[j]]
    }


  }

   return(globalloadmatrix)


}
