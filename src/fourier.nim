import field
import math
import std/assertions
import std/complex

proc `/=`[T](x: var Complex[T]; y: T) =
  x.re /= y
  x.im /= y

proc raw_fft(data: var open_array[Complex64], offset: int, stride: int, length: int) =
  if length == 1: return
  assert is_power_of_two length 

  var chunk = length
  while chunk >= 2:
    let half = chunk shr 1
    let angle = -2.0 * PI / float64(chunk)
    let factor_step = rect(1.0, angle)

    var i = 0
    while i < length:
      var factor = complex 1.0 
      var j = 0
      while j < half:
        let even = data[offset + stride * (i + j)]
        let odd = data[offset + stride * (i + j + half)]

        data[offset + stride * (i + j)] = even + odd
        data[offset + stride * (i + j + half)] = (even - odd) * factor

        factor *= factor_step
        
        inc j
      i += chunk
    chunk = chunk shr 1

  var i = 1
  var j = 0
  while i < length:
    var bit = length shr 1
    while (j and bit) != 0:
      j = j xor bit
      bit = bit shr 1
    j = j xor bit

    if i < j:
      let tmp = data[offset + stride * i]
      data[offset + stride * i] = data[offset + stride * j]
      data[offset + stride * j] = tmp

    inc i

## Performs FFT in-place with a negative exponent and 1/(4π^2) normalization
proc fft*[Nx, Ny: static int](field: var Field[Nx, Ny, Complex64]) =
  assert is_power_of_two Nx
  assert is_power_of_two Ny

  for row in 0..(Ny - 1):
    raw_fft(field.data, row * Nx, 1, Nx)
  for column in 0..(Nx - 1):
    raw_fft(field.data, column, Nx, Ny)

## Performs iFFT in-place with a positive exponent and no normalization
proc ifft*[Nx, Ny: static int](field: var Field[Nx, Ny, Complex64]) =
  assert is_power_of_two Nx
  assert is_power_of_two Ny

  for entry in field.data.mitems: entry = entry.conjugate

  for row in 0..(Ny - 1):
    raw_fft(field.data, row * Nx, 1, Nx)
  for column in 0..(Nx - 1):
    raw_fft(field.data, column, Nx, Ny)

  for entry in field.data.mitems: entry = entry.conjugate
  let norm = 1.0 / float64(Nx * Ny)
  for entry in field.data.mitems: entry = entry * norm
