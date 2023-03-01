hMVK <- function(t, A, B, delta){
  E <- exp((B - A) * t)
  out <- log(1 - E) - log(B * E - A)
  A * B * delta * exp(out)
}
pMVK <- function(t, A, B, delta){
  out <- log(B - A) + B * t - log(B * exp((B - A) * t) - A)
  1 - exp(delta * out)
}
dMVK <- function(t, A, B, delta){
  E <- exp((B - A) * t)
  out <- delta * (log(B - A) + B * t) + log(1 - E) - log(B * E - A) * (delta + 1)
  A * B * delta * exp(out)
}
