
#' Form global matrix from element matrix (stiffness, mass and stability)
#'
#' @param enumber           Number of discretized elements.
#' @param youngmod          Young's modulus.
#' @param shearmod          Shear modulus.
#' @param momentinertia     Moment of inertia (e.g. bh^3/12 for a rectangle).
#' @param area              Cross-sectional area
#' @param shearfactor       Shear correction factor (constant k).
#' @param poissonratio      Poisson's ratio.
#' @param totallength       Total length of the undiscretized structure
#' @param lengthscale       Material length scale parameter
#' @param rho               Mass density (kg/cubic meter)
#' @param case              Used to specify if stiffness matrix or mass matrix
#'                          or stability matrix. Note that case = 1 (stiffness),
#'                          or 2 (mass), or 3(stability))
#'
#'
#' @return                 Depending on what is being evaluated, it returns
#'                         the global stiffness, or mass or stability matrix
#' @export
#'

FormGlobalKMS=function(enumber, youngmod, shearmod,
                      momentinertia, area, shearfactor,
                      poissonratio, totallength, lengthscale,
                      rho, case){
  dof <- 4;
  dofpernode <- 2
  tdof <- (enumber + 1)*dofpernode;
  klist <- list();
  mlist <- list();
  slist <- list();
  length <- totallength/enum

  if(case == 1){

    for(j in seq(enumber)){

      klist[[j]]=FormStiffnessMTB(youngmod, shearmod,
                                  momentinertia, area,
                                  shearfactor, poissonratio,
                                  length, lengthscale);
    }


    globalmatrix <- matrix(0, tdof, tdof);
    Klist <- list();

    for(j in seq(enumber)){

      Klist[[j]] <- ExpandStiffnessMTB(tdof, klist[[j]], j, j+1)
      globalmatrix <- globalmatrix + Klist[[j]]
    }


  }

  # Mass matrix
  if(case == 2){

    for(j in seq(enumber)){

      mlist[[j]]=FormMassMTB(youngmod, shearmod,
                             momentinertia, area,
                             shearfactor, poissonratio,
                             length, rho);
    }

    globalmatrix <- matrix(0,tdof,tdof);
    Mlist <- list();

    for(j in seq(enumber)){

      Mlist[[j]] <- ExpandStiffnessMTB(tdof, mlist[[j]], j, j+1)
      globalmatrix <- globalmatrix + Mlist[[j]]
    }
  }

  # Stability matrix
  if(case == 3){

    for(j in seq(enumber)){

      slist[[j]]=FormStabilityMTB(youngmod, shearmod,
                                  momentinertia, area,
                                  shearfactor, poissonratio,
                                  length);
    }

    globalmatrix <- matrix(0, tdof, tdof);
    Slist <- list();

    for(j in seq(enumber)){

      Slist[[j]] <- ExpandStiffnessMTB(tdof, slist[[j]], j, j+1)
      globalmatrix <- globalmatrix + Slist[[j]]
    }
  }

  return(globalmatrix)


}
