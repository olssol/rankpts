#' Merge lists
#'
#'
#' @export
#'
tl_merge_lists <- function(lst_to, lst_from) {
    for (i in seq(length(lst_from))) {
        cur_name <- names(lst_from)[i]

        if (cur_name %in% names(lst_to))
            next

        lst_to[[cur_name]] <- lst_from[[cur_name]]
    }

    lst_to
}


#' Get X Beta by formula
#'
#'
#'
tl_xbeta <- function(data, fml, beta) {
    fml     <- as.formula(fml)
    des_mat <- model.matrix(fml, data = data)
    rst     <- apply(des_mat, 1, function(x) sum(beta * x))

    rst
}

#' Assign text to numeric vector
#'
#'
#' @export
#'
tkt_assign <- function(txt, prefix = "c(", suffix = ")") {

    rst <- NULL
    txt <- paste("rst <-", prefix, txt, suffix)

    tryCatch({
        eval(parse(text = txt))
    }, error = function(e) {
    })

    rst
}
