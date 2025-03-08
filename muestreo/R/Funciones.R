estra.media<-function (y, M, m, d = 5, alfa = 0.05, delta = 0.5)
{
  res <- estra(y, M, m)
  a <- 1
  L <- length(M)
  j <- 1
  u <- qnorm(1 - alfa/2)
  boxplot(y~estratificar(m),col=2,xlab="Estratos", main = "Box-plot múltiple")
  N <- sum(M)
  n <- length(y)
  k <- 1/sqrt(alfa)
  desm <- u * sqrt(res[2])
  desmch <- k * sqrt(res[2])
  cat("\n", "M.A.S.E.(", N, ",", n, ",", L)
  while (j <= L) {
    cat(",{", M[j], ",", m[j], "}")
    j <- j + 1
  }
  cat(")")
  cat("\n", "\n", "Estudio Inferencial", "\n")
  cat("MEDIA", "\n", "Estimación:", res[1], "\n", "Varianza:",
      res[2], "\n", "Error de muestreo:", sqrt(res[2]), "\n",
      "Intervalo de confianza al nivel", (1 - alfa) * 100,
      "%: (", res[1] - desm, ",", res[1] + desm, ")",
      "\n", "Intervalo de confianza al nivel", (1 - alfa) *
        100, "% s.n.: (", res[1] - desmch, ",", res[1] +
        desmch, ")", "\n")
  j <- 1
  while (j <= L) {
    num <- ((u^2) * M[j]/sum(M) * sqrt(var(y[a:sum(m[1:j])])) * res[5])
    den<- ((delta^2) + (u^2) * res[6]/sum(M))
    if (num/den < M[j]) b<-num/den else b<-M[j]
    cat(" Afijación de Neyman para el estrato", j, "con error absoluto",
        delta, "al nivel",(1-alfa)*100,"%:",ceiling(b), "\n")
    j <- j + 1
    a <- a + m[j-1]
  }
  estudio.des(y, d)
}
estra.total<-function (y, M, m, d = 5, alfa = 0.05, delta = 0.5)
{
  res <- estra(y, M, m)
  a <- 1
  L <- length(M)
  j <- 1
  u <- qnorm(1 - alfa/2)
  boxplot(y~estratificar(m),col=2,xlab="Estratos", main = "Box-plot múltiple")
  N <- sum(M)
  n <- length(y)
  k <- 1/sqrt(alfa)
  dest <- u * sqrt(res[4])
  destch <- k * sqrt(res[4])
  cat("\n", "M.A.S.E.(", N, ",", n)
  while (j <= L) {
    cat(",{", M[j], ",", m[j], "}")
    j <- j + 1
  }
  cat(")")
  cat("\n", "\n", "Estudio Inferencial", "\n")
  cat("TOTAL", "\n", "Estimación:", res[3], "\n", "Varianza:",
      res[4], "\n", "Error de muestreo:", sqrt(res[4]), "\n",
      "Intervalo de confianza al nivel", (1 - alfa) * 100,
      "%: (", res[3] - dest, ",", res[3] + dest, ")",
      "\n","Intervalo de confianza al nivel", (1 - alfa) *
        100, "% s.n.: (", res[3] - destch, ",", res[3] +
        destch, ")", "\n")
  j <- 1
  while (j <= L) {
    num <- (sum(M)^2)*((u^2) * M[j]/sum(M) * sqrt(var(y[a:sum(m[1:j])])) * res[5])
    den<- ((delta^2) + (u^2) * res[6]*sum(M))
    if (num/den < M[j]) b<-num/den else b<-M[j]
    cat(" Afijación de Neyman para el estrato", j, "con error absoluto",
        delta, "al nivel",(1-alfa)*100,"%:",ceiling(b), "\n")
    j <- j + 1
    a <- a + m[j-1]
  }
  estudio.des(y, d)
}
estra.razon<-function (numerador, denominador, M, m, alfa = 0.05)
{
  res <- razone(numerador, denominador, M, m)
  a <- 1
  L <- length(M)
  j <- 1
  u <- qnorm(1 - alfa/2)
  N <- sum(M)
  n <- length(numerador)
  if (n == length(denominador)) {
    desr <- u * sqrt(res[2])
    k <- 1/sqrt(alfa)
    desrch <- k * sqrt(res[2])
    cat("\n", "M.A.S.E.(", N, ",", n, ",", L)
    while (j <= L) {
      cat(",{", M[j], ",", m[j], "}")
      j <- j + 1
    }
    cat(")")
    cat("\n", "\n", "Estudio Inferencial", "\n")
    cat("RAZÓN", "\n", "Estimación:", res[1], "\n", "Varianza:",
        res[2], "\n", "Error de muestreo:", sqrt(res[2]),
        "\n", "Intervalo de confianza al nivel", (1 -
                                                    alfa) * 100, "%: (", res[1] - desr, ",",
        res[1] + desr, ")", "\n", "Intervalo de confianza al nivel",
        (1 - alfa) * 100, "% s.n.: (", res[1] - desrch,
        ",", res[1] + desrch, ")", "\n")
  }
  else {
    cat("\n", "Los vectores de datos no tienen la misma longitud.",
        "\n")
  }
}
estra<-function (y, M, m)
{
  L <- length(M)
  N <- sum(M)
  me <- 0
  varm <- 0
  a <- 1
  NE <- 0
  NEe <- 0
  for (i in 1:L) {
    su <- y[a:sum(m[1:i])]
    me <- M[i]/N * mean(su) + me
    varm <- (1 - m[i]/M[i])/m[i] * ((M[i]/N)^2) * var(su) +
      varm
    NE <- sqrt(var(su)) * M[i]/N + NE
    NEe <- var(su) * M[i]/N + NEe
    a <- a + m[i]
  }
  total <- N * me
  vart <- (N^2) * varm
  c(me, varm, total, vart, NE, NEe)
}
congl.media<-function (y, V, v, N, M, d = 5, alfa = 0.05)
{
  if (length(y) == sum(v)) {
    res <- conglomerado(y, V, v, N, M)
    u <- qnorm(1 - alfa/2)
    histograma(y)
    desm <- u * sqrt(res[2])
    k <- 1/sqrt(alfa)
    desmch <- k * sqrt(res[2])
    cat("\n", "Estudio Inferencial", "\n")
    cat("MEDIA", "\n", "Estimación:", res[1], "\n", "Varianza:",
        res[2], "\n", "Error de muestreo:", sqrt(res[2]),
        "\n", "Intervalo de confianza al nivel", (1 -
                                                    alfa) * 100, "%: (", res[1] - desm, ",",
        res[1] + desm, ")", "\n", "Intervalo de confianza al nivel",
        (1 - alfa) * 100, "% s.n.: (", res[1] - desmch,
        ",", res[1] + desmch, ")", "\n")
    estudio.des(y, d)
  }
  else {
    cat("\n", "El tamaño muestral no coincide con el tamaño de la muestra.",
        "\n")
  }
}
congl.total<-function (y, V, v, N, M, d = 5, alfa = 0.05)
{
  if (length(y) == sum(v)) {
    res <- conglomerado(y, V, v, N, M)
    u <- qnorm(1 - alfa/2)
    histograma(y)
    desm <- u * sqrt((N^2) * res[2])
    k <- 1/sqrt(alfa)
    desmch <- k * sqrt((N^2) * res[2])
    cat("\n", "Estudio Inferencial", "\n")
    cat("TOTAL", "\n", "Estimación:", N * res[1], "\n",
        "Varianza:", (N^2) * res[2], "\n", "Error de muestreo:",
        sqrt((N^2) * res[2]), "\n", "Intervalo de confianza al nivel",
        (1 - alfa) * 100, "%: (", N * res[1] - desm,
        ",", N * res[1] + desm, ")", "\n", "Intervalo de confianza al nivel",
        (1 - alfa) * 100, "% s.n.: (", N * res[1] -
          desmch, ",", N * res[1] + desmch, ")", "\n")
    estudio.des(y, d)
  }
  else {
    cat("\n", "El tamaño muestral no coincide con el tamaño de la muestra.",
        "\n")
  }
}
congl.razon<-function (numerador, denominador, V, v, N, M, alfa = 0.05)
{
  if (length(denominador) == length(numerador)) {
    if (length(numerador) == sum(v)) {
      res <- razoncon(numerador, denominador, V, v,
                      N, M)
      u <- qnorm(1 - alfa/2)
      desm <- u * sqrt(res[2])
      k <- 1/sqrt(alfa)
      desmch <- k * sqrt(res[2])
      cat("\n", "Estudio Inferencial", "\n")
      cat("RAZÓN", "\n", "Estimación:", res[1], "\n",
          "Varianza:", res[2], "\n", "Error de muestreo:",
          sqrt(res[2]), "\n", "Intervalo de confianza al nivel",
          (1 - alfa) * 100, "%: (", res[1] - desm,
          ",", res[1] + desm, ")", "\n", "Intervalo de confianza al nivel",
          (1 - alfa) * 100, "% s.n.: (", res[1] -
            desmch, ",", res[1] + desmch, ")", "\n")
    }
    else {
      cat("\n", "El tamaño muestral no coincide con el tamaño de la muestra.",
          "\n")
    }
  }
  else {
    cat("\n", "El tamaño muestral del numerador y del denominador NO coinciden.",
        "\n")
  }
}
conglomerado<-function (y, V, v, N, M)
{
  g <- length(v)
  a <- M * V[1]/N * mean(y[1:v[1]])
  b <- (V[1]/N)^2 * ((1/v[1]) - (1/V[1])) * var(y[1:v[1]])
  for (i in 2:g) {
    d <- y[(sum(v[1:(i - 1)]) + 1):(sum(v[1:i]))]
    a <- c(a, M * V[i]/N * mean(d))
    b <- c(b, (V[i]/N)^2 * ((1/v[1]) - (1/V[1])) * var(d))
    res <- rbind(a, b)
  }
  media <- mean(res[1, ])
  va <- ((1/g) - (1/M)) * var(res[1, ]) + M/g * sum(res[2,
  ])
  c(media, va)
}
razon<-function (numerador, denominador, N)
{
  n <- length(numerador)
  if (n == length(denominador)) {
    ra <- sum(numerador)/sum(denominador)
    varian <- (((1 - n/N)/n)/((sum(denominador)/n)^2)) *
      (var(numerador) + (ra^2) * var(denominador) -
         2 * ra * cov(numerador, denominador))
    res <- c(ra, varian)
    res
  }
  else {
    cat("\n", "Los vectores de datos no tienen la misma longitud.",
        "\n")
  }
}
razone<-function (numerador, denominador, M, m)
{
  L <- length(M)
  N <- sum(M)
  num <- numerador[1:m[1]]
  dem <- denominador[1:m[1]]
  ra <- razon(num, dem, M[1])[1] * M[1]/N
  va <- razon(num, dem, M[1])[2] * ((M[1]/N)^2)
  for (i in 2:L) {
    num <- numerador[(sum(m[1:(i - 1)]) + 1):sum(m[1:i])]
    dem <- denominador[(sum(m[1:(i - 1)]) + 1):sum(m[1:i])]
    ra <- ra + razon(num, dem, M[i])[1] * M[i]/N
    va <- va + razon(num, dem, M[i])[2] * ((M[i]/N)^2)
  }
  c(ra, va)
}
razoncon<-function (numerador, denominador, V, v, N, M)
{
  if (length(numerador) == length(denominador)) {
    num <- conglomerado(numerador, V, v, N, M)
    den <- conglomerado(denominador, V, v, N, M)
    ra <- num[1]/den[1]
    g <- length(v)
    a <- M * V[1]/N * (mean(numerador[1:v[1]]) - ra * mean(denominador[1:v[1]]))
    b <- (V[1]/N)^2 * ((1/v[1]) - (1/V[1])) * (var(numerador[1:v[1]]) +
                                                 (ra^2) * var(denominador[1:v[1]]) - 2 * ra * cov(numerador[1:v[1]],
                                                                                                  denominador[1:v[1]]))
    for (i in 2:g) {
      d <- M * V[i]/N * (mean(numerador[(sum(v[1:(i - 1)]) +
                                           1):(sum(v[1:i]))]) - ra * mean(denominador[(sum(v[1:(i -
                                                                                                  1)]) + 1):(sum(v[1:i]))]))
      a <- c(a, d)
      dd <- (V[i]/N)^2 * ((1/v[i]) - (1/V[i])) * (var(numerador[(sum(v[1:(i -
                                                                            1)]) + 1):(sum(v[1:i]))]) + (ra^2) * var(denominador[(sum(v[1:(i -
                                                                                                                                             1)]) + 1):(sum(v[1:i]))]) - 2 * ra * cov(numerador[(sum(v[1:(i -
                                                                                                                                                                                                            1)]) + 1):(sum(v[1:i]))], denominador[(sum(v[1:(i -
                                                                                                                                                                                                                                                              1)]) + 1):(sum(v[1:i]))]))
      b <- c(b, dd)
      res <- rbind(a, b)
    }
    va <- (1/(den[1]^2)) * (((1/g) - (1/M)) * var(res[1,
    ]) + M/g * sum(res[2, ]))
    c(ra, va)
  }
  else {
    cat("\n", "NO coincide el tamaño muestral de ambas muestra.",
        "\n")
  }
}
estudio.des<-function (y, d = 5)
{
  library(e1071)
  cat("\n", "Estudio Descriptivo de la muestra", "\n", "Parámetros de interés",
      "\n", "Tamaño muestral:", length(y), "\n", "Media muestral:",
      mean(y), "\n", "Mediana muestral:", median(y), "\n",
      "Varianza muestral:", var(y), "\n", "Desviación típica:",
      sqrt(var(y)), "\n", "Coeficiente de variación:", sqrt(var(y))/mean(y) *
        100, "\n", "Rango intercuartílico:", IQR(y), "\n",
      "Coeficiente de asimetría:", skewness(y), "\n", "Coeficiente de curtosis:",
      kurtosis(y), "\n", "\n", "La información se resume en la siguiente tabla:",
      "\n")
  if (length(y) <= 50) {
    stem(y)
  }
  else {
    tabla(y, d)
  }
}
jacknife<-function (y, N)
{
  n <- length(y)
  varm <- 0
  vart <- 0
  resa <- rep(0, n)
  resb <- rep(0, n)
  for (i in 1:n) {
    a <- sum(y) - (n - 1) * mean(y[-i])
    b <- n * sum(y) - (n - 1) * sum(y[-i])
    resa[i] <- a
    resb[i] <- b
  }
  for (i in 1:n) {
    jacha <- sum(resa)/n
    jachb <- sum(resb)/n
    varm <- varm + 1/(n * (n - 1)) * ((resa[i] - jacha)^2)
    vart <- vart + 1/(n * (n - 1)) * ((resb[i] - jachb)^2)
  }
  c(varm, vart)
}
tabla<-function (x, cl)
{
  f <- table(cut(x, seq(min(x) - 0.001, max(x) + 0.001, length = cl +
                          1)))
  h <- f/length(x)
  F <- cumsum(f)
  H <- cumsum(h)
  res <- cbind(f, h, F, H)
  colnames(res) <- c("  f", "  h", "  F", "  H")
  cat("\nTabla de Frecuencias\n\n")
  res
}
histograma<-function(y){
  #library(Mass)
  truehist(y)
  lines(density(y))
  title(main="Histograma de frecuencias y densidad")
}
estratificar<-function (n)
{
  a <- length(n)
  res <- rep(1, n[1])
  for (i in 2:a) {
    res <- c(res, rep(i, n[i]))
  }
  res
}
truehist<-function (data, nbins = nclass.scott(data), h, x0 = -h/1000,
                    breaks, prob = TRUE, xlim = range(breaks), ymax = max(est),
                    col = 5, xlab = deparse(substitute(data)), bty = "n", ...)
{
  eval(xlab)
  data <- data[!is.na(data)]
  if (missing(breaks)) {
    if (missing(h))
      h <- diff(pretty(data, nbins))[1]
    first <- floor((min(data) - x0)/h)
    last <- ceiling((max(data) - x0)/h)
    breaks <- x0 + h * c(first:last)
  }
  if (any(diff(breaks) <= 0))
    stop("breaks must be strictly increasing")
  if (min(data) < min(breaks) || max(data) > max(breaks))
    stop("breaks do not cover the data")
  db <- diff(breaks)
  if (!prob && sqrt(var(db)) > mean(db)/1000)
    warning("Uneven breaks with prob = FALSE will give a misleading plot")
  bin <- cut(data, breaks, include.lowest = TRUE)
  est <- tabulate(bin, length(levels(bin)))
  if (prob)
    est <- est/(diff(breaks) * length(data))
  n <- length(breaks)
  plot(xlim, c(0, ymax), type = "n", xlab = xlab, ylab = "",
       bty = bty)
  rect(breaks[-n], 0, breaks[-1], est, col = col, ...)
  invisible()
}
nclass.scott<-function (x)
{
  h <- 3.5 * sqrt(var(x)) * length(x)^(-1/3)
  ceiling(diff(range(x))/h)
}
