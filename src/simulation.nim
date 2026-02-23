import std/complex
from math import cos, sin, PI, exp
import pixie
import linalg
import field
import fourier

const tolerance = 1e-7

proc k(i, N: int): float64 =
  if i < N div 2: float64(i)
  elif i == N div 2: 0.0
  else: float64(i - N)

proc dealias[Nx, Ny: static int](field: var Field[Nx, Ny, Complex64]) =
  const Kx = Nx div 3
  const Ky = Ny div 3
  for x, y in field.indices:
    let kx = (if x < Nx div 2: x else: x - Nx)
    let ky = (if y < Ny div 2: y else: y - Ny)
    if kx.abs > Kx or ky.abs > Ky:
      field[x, y] = complex(0.0)
type
  Simulation[Nx, Ny: static int] = ref object
    fourier_vorticity*: Field[Ny, Nx, Complex64]
    Δt: float64
    kinematic_viscosity: float64

proc create_simulation*[Nx, Ny: static int](
  vorticity: Field[Ny, Nx, float64],
  kinematic_viscosity: float64,
  Δt: float64,
): Simulation[Nx, Ny] =
  var fourier_vorticity = createField[Nx, Ny, Complex64](complex(0.0))
  for i in low(FieldIndex[Nx, Ny])..high(FieldIndex[Nx, Ny]):
    let idx: FieldIndex[Nx, Ny] = i
    fourier_vorticity[i] = complex(vorticity[idx])
  fft fourier_vorticity 
  dealias fourier_vorticity

  result = Simulation[Nx, Ny](
    fourier_vorticity: fourier_vorticity,
    kinematic_viscosity: kinematic_viscosity,
    Δt: Δt,
  )

proc calculate_velocity[Nx, Ny: static int](
  fourier_ω: Field[Nx, Ny, Complex64]
): Field[Nx, Ny, Vec2f] =
  var velocity_x = createField[Nx, Ny, Complex64](complex(0.0))
  var velocity_y = createField[Nx, Ny, Complex64](complex(0.0))
  
  for x, y in fourier_ω.indices:
    let k = (x: k(x, Nx), y: k(y, Ny))
    if k.x == 0.0 and k.y == 0.0: continue
    let ω = fourier_ω[x, y]
    velocity_x[x, y] = im(1.0) * k.y / k.abs2 * ω
    velocity_y[x, y] = im(1.0) * -k.x / k.abs2 * ω

  ifft velocity_x
  ifft velocity_y

  var velocity = createField[Nx, Ny, Vec2f]((x: 0.0, y: 0.0))
  for x, y in velocity.indices:
    assert(velocity_x[x, y].im.abs < tolerance, msg = $velocity_x[x, y].im.abs)
    assert(velocity_y[x, y].im.abs < tolerance, msg = $velocity_y[x, y].im.abs)
    velocity[x, y] = (
      x: velocity_x[x, y].re,
      y: velocity_y[x, y].re,
    )

  result = velocity

proc calculate_nonlinear[Nx, Ny: static int](
  velocity: Field[Nx, Ny, Vec2f],
  fourier_vorticity: Field[Nx, Ny, Complex64],
): Field[Nx, Ny, Complex64] =
  var test_ω = create_field[Nx, Ny, Complex64](complex(0.0))
  for x, y in fourier_vorticity.indices:
    test_ω[x, y] = fourier_vorticity[x, y]

  ifft test_ω
  
  for c in test_ω.data:
    assert(c.im.abs < tolerance, msg = $c.im.abs)

  var nonlinear = create_field[Nx, Ny, Complex64](complex(0.0))

  var ∇ω_x = create_field[Nx, Ny, Complex64](complex(0.0))
  var ∇ω_y = create_field[Nx, Ny, Complex64](complex(0.0))

  for x, y in ∇ω_x.indices:
    ∇ω_x[x, y] = fourier_vorticity[x, y] * im(1.0) * k(x, Nx)

  for x, y in ∇ω_y.indices:
    ∇ω_y[x, y] = fourier_vorticity[x, y] * im(1.0) * k(y, Ny)

  ifft ∇ω_x
  ifft ∇ω_y

  for x, y in nonlinear.indices:
    let v = velocity[x, y]
    assert(∇ω_x[x, y].im.abs < tolerance, msg = $∇ω_x[x, y].im.abs)
    assert(∇ω_y[x, y].im.abs < tolerance, msg = $∇ω_y[x, y].im.abs)
    nonlinear[x, y] = complex(
      -(v.x * ∇ω_x[x, y].re + v.y * ∇ω_y[x, y].re)
    )

  fft nonlinear
  dealias nonlinear
  result = nonlinear

proc step*[Nx, Ny: static int](simulation: Simulation[Nx, Ny]) =
  let Δt = simulation.Δt
  let kinematic_viscosity = simulation.kinematic_viscosity
  var fourier_ω = simulation.fourier_vorticity
  let velocity = calculate_velocity simulation.fourier_vorticity
  let nonlinear = calculate_nonlinear(velocity, simulation.fourier_vorticity)
  
  for x, y in fourier_ω.indices:
    let k = (
      x: k(x, Nx),
      y: k(y, Ny),
    )

    if k.x == 0 and k.y == 0:
      fourier_ω[x, y] = complex(0.0)
      continue
    
    let E = exp(-kinematic_viscosity * k.abs2 * Δt)
    let ω = fourier_ω[x, y]
    let B = (1.0 - E) / (kinematic_viscosity * k.abs2)
    let N = nonlinear[x, y]

    fourier_ω[x, y] = E * ω + B * N
  dealias fourier_ω

  var test_ω = create_field[Nx, Ny, Complex64](complex(0.0))
  for x, y in fourier_ω.indices:
    test_ω[x, y] = fourier_ω[x, y]

  ifft test_ω
  
  for c in test_ω.data:
    assert(c.im.abs < tolerance, msg = $c.im.abs)

  simulation.fourier_vorticity = fourier_ω

proc visualize*[Nx, Ny: static int](simulation: Simulation[Nx, Ny]): Image =
  var ω = create_field[Nx, Ny, Complex64](complex(0.0))
  let fourier_ω = simulation.fourier_vorticity
  for x, y in fourier_ω.indices:
    ω[x, y] = fourier_ω[x, y]
  ifft ω

  for x, y in ω.indices:
    assert ω[x, y].im.abs < tolerance

  var max_abs = 0.0
  for c in ω.data: max_abs = if c.re.abs > max_abs: c.re.abs else: max_abs

  echo "max_abs: ", $max_abs

  var image = newImage(Nx, Ny)
  for x, y in ω.indices:
    let re = ω[x, y].re
    let a = re / max_abs
    let c = 
      if a >= 0: rgba(uint8(a * 255.0), 0, 0, 255)
      else: rgba(0, 0, uint8(-a * 255.0), 255)
    image.data[y * Nx + x] = c

  result = image
