#' @title Función mas.razon
#' @description La función mas.razon es utilizada para estudiar el parámetro razón de los totales poblacionales de X respecto a Y, cuando se realiza un muestreo aleatorio simple.



#' @param numerador: vector numérico de las observaciones realizadas para la variable X.
#' @param denominador: vector numérico de las observaciones realizadas para la variable Y.
#' @param N: tamaño de la población muestreada.
#' @param alfa: nivel de error (1-alfa: nivel de confianza).
#' @return
#'
#' Estudio Inferencial para el parámetro razón poblacional a partir de un muestreo aleatorio simple (M.A.S.):
#'
#' -Estimación: estimación del parámetro razón (estimador: razón muestral).
#'
#' -Varianza: estimación de la varianza de la razón muestral.
#'
#' -Error de muestreo: estimación de la desviación típica de la razón muestral.
#'
#'           -Intervalos de confianza: intervalos de confianza al nivel 1-alfa para el parámetro razón poblacional,
#'           bajo los supuestos de normalidad y no normalidad (s.n.-sin normalidad) en la población.
#'
#

#' @export mas.razon
#' @examples
#' mas.razon(rnorm(1000,75,3), rnorm(1000,54,1), N=10000, alfa=0.01)




mas.razon<-function (numerador, denominador, N, alfa = 0.05)
{
  n <- length(numerador)
  if (n == length(denominador)) {
    res <- razon(numerador, denominador, N)
    raz <- res[1]
    varr <- res[2]
    u <- qnorm(1 - alfa/2)
    desr <- u * sqrt(varr)
    k <- 1/sqrt(alfa)
    desrch <- k * sqrt(varr)
    cat("\n", "M.A.S.(", N, ",", n, ")", "\n", "\n", "Estudio Inferencial",
        "\n")
    cat("RAZÓN", "\n", "Estimación:", raz, "\n", "Varianza:",
        varr, "\n", "Error de muestreo:", sqrt(varr),
        "\n", "Intervalo de confianza al nivel", (1 -
                                                    alfa) * 100, "%: (", raz - desr, ",",
        raz + desr, ")", "\n", "Intervalo de confianza al nivel",
        (1 - alfa) * 100, "% s.n.:(", raz - desrch,
        ",", raz + desrch, ")", "\n")
  }
  else {
    cat("\n", "Los vectores de datos no tienen la misma longitud.",
        "\n")
  }
}
