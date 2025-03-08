#' @title Función ms.media
#' @description La función ms.media es utilizada para estudiar el parámetro media de una variable cuando se realiza un muestreo sistemático. También
#' es válida para proporciones utilizando variables dicotómicas.



#' @param y: vector numérico de las observaciones realizadas.
#' @param N: tamaño de la población muestreada.
#' @param k: paso utilizado para extraer la muestra (N=k*n, siendo n el tamaño muestral).
#' @param d: número de intervalos o clases utilizados para construir el histograma de los datos observados y la tabla de frecuencias, cuando corresponda.
#' @param alfa: nivel de error (1-alfa: nivel de confianza).
#' @return
#'
#' Estudio Inferencial para el parámetro media poblacional a partir de un muestreo sistemático:
#'
#' -Estimación: estimación del parámetro media (estimador: media muestral).
#'
#' -Varianza: estimación de la varianza de la media muestral. Método Jackknife.
#'
#' -Error de muestreo: estimación de la desviación típica de la media muestral.
#'
#' -Intervalos de confianza: intervalos de confianza al nivel 1-alfapara el parámetro media poblacional, bajo los supuestos
#' de normalidad y no normalidad (s.n.-sin normalidad) en la población.
#'
#' Estudio Descriptivo:
#'
#' -Gráfico: Histograma de frecuencias absolutas y estimación de la función de densidad.
#'
#' -Estadísticos muestrales: Tamaño, media, mediana, varianza, desviación típica, coeficiente de variación, rango
#' intercuartílico, coeficiente de asimetría y coeficiente de curtosis.
#'
#' -Tabulación: Tabla de frecuencias (muestras mayores de 50) ó diagrama tallo-hoja (en otro caso).



#' @export ms.media
#' @examples
#' ms.media(rnorm(1000,75,3), N=10000, k=10, d=20)
#'





ms.media<-function (y, N, k, d = 5, alfa = 0.05)
{
  n <- length(y)
  histograma(y)
  if (k != (N/n)) {
    cat("\n", "El paso", k, "no es múltiplo del tamaño poblacional:",
        N, "\n")
  }
  else {
    n <- length(y)
    total <- sum(y)
    me <- total/n
    varm <- jacknife(y, N)[1]
    vart <- jacknife(y, N)[2]
    u <- qnorm(1 - alfa/2)
    desm <- u * sqrt(varm)
    dest <- u * sqrt(vart)
    K <- 1/sqrt(alfa)
    desmch <- K * sqrt(varm)
    destch <- K * sqrt(vart)
    cat("\n", "M.S.(", N, ",", k, ")", "\n", "\n", "Estudio Inferencial",
        "\n")
    cat("MEDIA", "\n", "Estimación:", me, "\n", "Varianza:",
        varm, "\n", "Error de muestreo:", sqrt(varm),
        "\n", "Intervalo de confianza al nivel", (1 -
                                                    alfa) * 100, "%: (", me - desm, ",", me +
          desm, ")", "\n", "Intervalo de confianza al nivel",
        (1 - alfa) * 100, "% s.n.:(", me - desmch,
        ",", me + desmch, ")", "\n")
    estudio.des(y, d)
  }
}
