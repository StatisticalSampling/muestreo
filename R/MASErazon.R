#' @title Función mase.razon
#' @description    La función mase.razon es utilizada para estudiar el parámetro razón
#' de los totales poblacionales de X respecto a Y, cuando se realiza un
#' muestreo aleatorio simple estratificado. Así mismo puede ser utilizada
#' para realizar un estudio de la razón en cada uno de los estratos
#' considerados (mediante un m.a.s.).



#' @param numerador: vector numérico de las observaciones realizadas para la variable, X.
#' @param denominador: vector numérico de las observaciones realizadas para la variable, Y.
#' @param M: vector que contiene los tamaños poblacionales de los estratos.
#' @param m: vector que contiene los tamaños de las muestras extraídas de cada estrato.
#' @param alfa: nivel de error (1-alfa: nivel de confianza).
#' @param estrato: número del estrato para el que se quiere realizar el estudio. (estrato = 0 indica que se estudia toda la muestra
#' conjuntamente.)
#' @return
#'
#' Valor:
#'
#'   -estrato = 0
#'
#'   Estudio Inferencial para el parámetro razón poblacional a partir
#'   de un muestreo aleatorio simple estratificado (M.A.S.E.):
#'
#'   -Estimación: estimación del parámetro razón (estimador: razón
#'                          ponderada muestral).
#'
#'       -Varianza: estimación de la varianza de la razón ponderada muestral.
#'
#'       -Error de muestreo: estimación de la desviación típica de la razón ponderada muestral.
#'        -Intervalos de confianza: intervalos de confianza al nivel 1-alfa para el parámetro razón poblacional, bajo los supuestos
#'        de normalidad y no normalidad (s.n.-sin normalidad) en la población.
#'
#'        -estrato = i (>0)
#'
#'        Estudio Inferencial para el parámetro razón poblacional (del estrato
#'                                                                   i-ésimo) a partir de un muestreo aleatorio simple (M.A.S.):
#'
#'      -Estimación: estimación del parámetro razón (estimador: razón muestral).
#'
#'      -Varianza: estimación de la varianza de la razón muestral.
#'
#'      -Error de muestreo: estimación de la desviación típica de la razón muestral.
#'
#'      -Intervalos de confianza: intervalos de confianza al nivel 1-alfa para el parámetro razón poblacional, bajo los supuestos de normalidad y no normalidad (s.n.-sin normalidad) en
#'        la población.


#' @export mase.razon
#' @examples
#' mase.razon(x<-c(rnorm(1000,50,3),rnorm(2000,45,2),rnorm(500,51,4)),
#' y<-c(rnorm(1000,30,1.5),rnorm(2000,55,2),rnorm(500,70,5)),
#' M=c(10000,20000,5000), m=c(1000,2000,500))





mase.razon<-function (numerador, denominador, M, m, alfa = 0.05, estrato = 0)
{
  if ((sum(m) != length(numerador)) | (length(M) != length(m))) {
    cat("\n", "El tamaño muestral no coincide con el tamaño de la muestra",
        "\n")
  }
  else {
    if (length(m) < estrato) {
      cat("\n", "El número de estrato no es correcto", "\n")
    }
    else {
      if (estrato == 0) {
        estra.razon(numerador, denominador, M, m,
                    alfa)
      }
      else {
        if (estrato != 1) {
          a <- sum(m[1:estrato - 1]) + 1
          b <- a + m[estrato] - 1
          mas.razon(numerador[a:b], denominador[a:b],
                    M[estrato], alfa)
        }
        else {
          mas.razon(numerador[1:m[1]], denominador[1:m[1]],
                    M[1], alfa)
        }
      }
    }
  }
}
