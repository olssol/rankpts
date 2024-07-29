## p     = P(tox);
## q     = P(response);
#' Phase I functions
#'
#' @param pqr p q and odds ratio (p00p11/p10/p01)
#'
#' @export
#'
stb_tl_p1_p11 <- function(pqr) {
    p   <- pqr[1]
    q   <- pqr[2]
    rho <- pqr[3]

    if (1 == rho) {
        p11 <- p * q
        p01 <- (1 - p) * q
        p10 <- p * (1 - q)
        p00 <- (1 - p) * (1 - q)
    } else {
        tmp1 <- (p + q - p * rho - q * rho - 1)^2
        tmp2 <- 4 * (rho - 1) * p * q * rho
        tmp3 <- p + q - p * rho - q * rho - 1
        p11 <- -(sqrt(tmp1 - tmp2) + tmp3) / 2 / (rho - 1)
        p01 <- q - p11
        p10 <- p - p11
        p00 <- p01 * p10 * rho / p11
    }

    ##check
    if (0) {
        cp   <- p11 + p10
        cq   <- p11 + p01
        crho <- p11 * p00 / p01 / p10
        print(c(cp, cq, crho))
    }

    rst <- c(p00,p01,p10,p11)
    rst <- rst / sum(rst)

    rst
}

#' Phase I function random dirichelet
#'
#'
#'
#' @export
#'
stb_tl_p1_rdir <- function(n, alpha) {
    l  <- length(alpha)
    y  <- rgamma(l * n, alpha)
    x  <- matrix(y, ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x / as.vector(sm))
}

#' Set simulation scenario
#'
#' Simulation function. Get true \eqn{\theta}'s using marginal probabilities and
#' odds ratio \eqn{\rho} for all dose levels.
#'
#' @rdname vtScenario
#'
#' @param tox Vector of marginal DLT risk rates for all levels
#' @param res Vector of marginal immune response rates for all levels
#' @param rho Vector of odds ratio for all levels. If length of \code{rho} is
#'     shorter than the length of \code{tox} or \code{res}, vector \code{rho} is
#'     repeated to have the same length as \code{tox} and \code{res}.
#'
#' @details
#'
#' The calculation is as following. If \eqn{\rho = 1}, then \eqn{\theta_{11} =
#' pq}, \eqn{\theta_{01} = (1-p)q}, \eqn{\theta_{10} = p(1-q)}, and
#' \eqn{\theta_{00} = (1-p)(1-q)}. Otherwise, \eqn{ \theta_{11} = -(\sqrt{A+B}},
#' \eqn{\theta_{01} = q-\theta_{11}}, \eqn{\theta_{10} = p-\theta_{11}}, and
#' \eqn{\theta_{00} = \theta_{01}\theta_{10}\rho/\theta_{11}}, where
#' \eqn{A=(p+q-p \rho-q\rho-1)^2-4(\rho-1)pq\rho)} and
#' \eqn{B=(p+q-p\rho-q\rho-1))/2/(\rho-1)}.
#'
#'
#' @return a \code{VTTRUEPS} object containing all \eqn{\theta}'s in a matrix
#'     with its number of rows equaling the number of dose levels and its number
#'     of columns being 4.
#'
#' @examples
#'  rst.sce <- vtScenario(tox=c(0.05, 0.05, 0.08), res=c(0.2, 0.3, 0.5), rho=1)
#'
#' @export
#'
stb_tl_p1_sce <- function(tox = c(0.05, 0.05, 0.08),
                          res = c(0.2, 0.3, 0.5),
                          rho = 1) {

    if (!is.vector(res)) {
        res <- rep(res, length(tox))
    }

    if (!is.vector(rho)) {
        rho <- rep(rho, length(tox))
    }

    true_theta <- apply(cbind(tox, res, rho),
                        1,
                        function(x) {
        tt      <- stb_tl_p1_p11(x)
        tt[4]   <- 1 - sum(tt[1:3])
        tt
    })

    rst <- t(true_theta)

    colnames(rst) <- c("TR00", "TR01", "TR10", "TR11")
    invisible(rst)
}
