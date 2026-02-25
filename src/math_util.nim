import std/complex

proc `/=`*[T](x: var Complex[T]; y: T) =
  x.re /= y
  x.im /= y
proc `*=`*[T](x: var Complex[T]; y: T) =
  x.re *= y
  x.im *= y
