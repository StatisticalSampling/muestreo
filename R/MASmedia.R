#' @title Función mas.media
#' @description La función mas.media es utilizada para estudiar el parámetro media de una variable cuando se realiza un muestreo aleatorio simple.
#' También es válida para proporciones utilizando variables dicotómicas.


#' @param y: vector numérico de las observaciones realizadas.
#' @param N: tamaño de la población muestreada.
#' @param d: número de intervalos o clases utilizados para construir el histograma de los datos observados y la tabla de frecuencias, cuando corresponda.
#' @param alfa: nivel de error (1-alfa: nivel de confianza).
#' @param delta: error relativo.
#' @return
#'
#' Estudio Inferencial para el parámetro media poblacional a partir de un muestreo aleatorio simple (M.A.S.):
#'
#' -Estimación: estimación del parámetro media (estimador: media muestral).
#'
#' -Varianza: estimación de la varianza de la media muestral.
#'
#' -Error de muestreo: estimación de la desviación típica de la media muestral.
#'
#'           -Intervalos de confianza: intervalos de confianza al nivel 1-alfa para el parámetro media poblacional,
#'           bajo los supuestos de normalidad y no normalidad (s.n.-sin normalidad) en la población.
#'
#' -Tamaño muestral: tamaño muestral necesario para cometer un error relativo delta para un nivel de confianza 1-alfa,
#'  bajo los supuestos de normalidad y no normalidad (s.n.-sin normalidad) en la población.
#'
#'
#
#'
#'Estudio Descriptivo:
#'
#' -Gráfico: Histograma de frecuencias absolutas y estimación de la
#'función de densidad.
#'
#' -Estadísticos muestrales: Tamaño, media, mediana, varianza,
#' desviación típica, coeficiente de variación, rango
#' intercuartílico, coeficiente de asimetría y coeficiente
#' de curtosis.
#'
#'  -Tabulación: Tabla de frecuencias (muestras mayores de 50) ó diagrama tallo-hoja (en otro caso).
#'
#'
#' @export mas.media
#' @examples
#' mas.media(rnorm(1000,75,3), N=10000, d=20, delta=0.001)
#'




mas.media<-function (y, N, d = 5, alfa = 0.05, delta = 0.5)
{
  histograma(y)
  n <- length(y)
  fre <- n/N
  total <- sum(y)
  me <- total/n
  varm <- (1 - fre)/n * var(y)
  vart <- (N^2) * varm
  u <- qnorm(1 - alfa/2)
  desm <- u * sqrt(varm)
  dest <- u * sqrt(vart)
  k <- 1/sqrt(alfa)
  desmch <- k * sqrt(varm)
  destch <- k * sqrt(vart)
  na <- (u^2) * var(y)/mean(y)^2/(delta^2)
  nas <- (k^2) * var(y)/mean(y)^2/(delta^2)
  if(ceiling(na/(1 + (na/N))) < N) tamanom <- ceiling(na/(1 + (na/N))) else tamanom <- N
  if(ceiling(nas/(1 + (nas/N))) < N) tamanoms <- ceiling(nas/(1 + (nas/N))) else tamanoms <- N
  cat("\n", "M.A.S.(", N, ",", n, ")", "\n", "\n", "Estudio Inferencial",
      "\n\n")
  cat("MEDIA", "\n", "Estimación:", me, "\n", "Varianza:",
      varm, "\n", "Error de muestreo:", sqrt(varm), "\n", "Intervalo de confianza al nivel",
      (1 - alfa) * 100, "%: (", me - desm, ",", me + desm,
      ")", "\n", "Intervalo de confianza al nivel", (1 - alfa) *
        100, "% s.n.:(", me - desmch, ",", me + desmch,
      ")", "\n", "Tamaño muestral con error relativo", delta,
      "al nivel",(1-alfa)*100,"%:", tamanom, "\n", "Tamaño muestral s.n. con error relativo",
      delta, "al nivel",(1-alfa)*100,"%:", tamanoms, "\n")
  estudio.des(y, d)
}
