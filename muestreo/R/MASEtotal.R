#' @title Función mase.total
#' @description   La función mase.total es utilizada para estudiar el parámetro total
#' de una variable cuando se realiza un muestreo aleatorio simple
#' estratificado. Así mismo puede ser utilizada para realizar un estudio
#' del total en cada uno de los estratos considerados (mediante un m.a.s.).



#' @param y: vector numérico de las observaciones realizadas.
#' @param 	M: vector que contiene los tamaños poblacionales de los estratos.
#' @param m: vector que contiene los tamaños de las muestras extraídas de cada estrato.
#' @param alfa: nivel de error (1-alfa: nivel de confianza).
#' @param estrato: número del estrato para el que se quiere realizar el estudio (estrato = 0 indica que se estudia toda la muestra conjuntamente).
#' @param delta: error relativo (para estrato > 0).
#' @return
#'
#'
#'Valor:
#'
#'-estrato = 0
#'
#'Estudio Inferencial para el parámetro total poblacional a partir
#'de un muestreo aleatorio simple estratificado (M.A.S.E.):
#'
#'-Estimación: estimación del parámetro total (estimador: suma
#'                                              ponderada muestral).
#'
#'-Varianza: estimación de la varianza de la suma ponderada muestral.
#'
#'-Error de muestreo: estimación de la desviación típica de la suma ponderada muestral.
#'
#'-Intervalos de confianza: intervalos de confianza al nivel 1-alfa para el parámetro total poblacional, bajo los supuestos
#'de normalidad y no normalidad (s.n.-sin normalidad) en la población.
#'
#'Estudio Descriptivo:
#'
#'-Gráfico: Box-plot múltiple.
#'
#'-Estadísticos muestrales: Tamaño, media, mediana, varianza, desviación típica, coeficiente de variación, rango
#'intercuartílico, coeficiente de asimetría y coeficiente de curtosis.
#'
#'-Tabulación: Tabla de frecuencias (muestras mayores de 50) ó diagrama tallo-hoja (en otro caso).
#'
#'-estrato = i (>0)
#'
#'Estudio Inferencial para el parámetro total poblacional (del estrato
#'                                                         i-ésimo) a partir de un muestreo aleatorio simple (M.A.S.):
#'
#'-Estimación: estimación del parámetro total (estimador: suma
#'                                               muestral).
#'
#'-Varianza: estimación de la varianza de la suma muestral.
#'
#'-Error de muestreo: estimación de la desviación típica de la suma muestral.
#'
#'-Intervalos de confianza: intervalos de confianza al nivel 1-alfa
#'para el parámetro total poblacional, bajo los supuestos
#'de normalidad y no normalidad (s.n.-sin normalidad) en la población.
#'
#' -Tamaño muestral: tamaño muestral necesario para cometer un error
#'relativo delta para un nivel de confianza 1-alfa, bajo
#'los supuestos de normalidad y no normalidad (s.n.-sin
#'                                             normalidad) en la población.
#'
#' Estudio Descriptivo:
#'
#' -Gráfico: Histograma de frecuencias absolutas y estimación de la
#' función de densidad.
#'
#' -Estadísticos muestrales: Tamaño, media, mediana, varianza,
#' desviación típica, coeficiente de variación, rango
#' intercuartílico, coeficiente de asimetría y coeficiente de curtosis.
#'
#' -Tabulación: Tabla de frecuencias (muestras mayores de 50) ó diagrama
#' tallo-hoja (en otro caso).
#'
#'
#' @export mase.total
#' @examples
#' 	mase.total(x<-c(rnorm(1000,50,3),rnorm(2000,45,2),rnorm(500,51,4)),
#' 	M=c(10000,20000,5000), m=c(1000,2000,500),d=20)





mase.total<-function (y, M, m, d = 5, alfa = 0.05, delta = 0.5, estrato = 0)
{
  if ((sum(m) != length(y)) | (length(M) != length(m))) {
    cat("\n", "El tamaño muestral no coincide con el tamaño de la muestra",
        "\n")
  }
  else {
    if (length(m) < estrato) {
      cat("\n", "El n? de estrato no es correcto", "\n")
    }
    else {
      if (estrato == 0) {
        estra.total(y, M, m, d, alfa, delta)
      }
      else {
        if (estrato != 1) {
          mas.total(y[(sum(m[1:(estrato - 1)]) + 1):sum(m[1:estrato])],
                    M[estrato], d, alfa, delta)
        }
        else {
          mas.total(y[1:m[1]], M[1], d, alfa, delta)
        }
      }
    }
  }
}
