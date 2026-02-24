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

proc smooth_gaussian*[Nx, Ny: static int](
  sharp_mask: Field[Nx, Ny, float64],
  smoothing: float64,
): Field[Nx, Ny, float64] =
  var smooth_mask = create_field[Nx, Ny, float64](0.0)
  var fourier_mask = create_field[Nx, Ny, Complex64](complex(0.0))
  for x, y in sharp_mask.indices:
    fourier_mask[x, y] = complex(sharp_mask[x, y])
  fft fourier_mask
  for x, y in fourier_mask.indices:
    let kx = k(x, Nx)
    let ky = k(y, Ny)
    fourier_mask[x, y] = fourier_mask[x, y] * exp(-smoothing * (kx * kx / Nx / Nx + ky * ky / Ny / Ny))
  ifft fourier_mask
  for x, y in smooth_mask.indices:
    smooth_mask[x, y] = fourier_mask[x, y].re
  result = smooth_mask

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
    fourier_vorticity*: Field[Nx, Ny, Complex64]
    nonlinear: Field[Ny, Nx, Complex64]
    object_velocity*: Field[Ny, Nx, Vec2f]
    mask*: Field[Ny, Nx, float64]
    kinematic_viscosity: float64
    permeability: float64
    cfl: float64
    time*: float64
    Δt: float64

proc create_simulation*[Nx, Ny: static int](
  vorticity: Field[Ny, Nx, float64],
  object_velocity: Field[Ny, Nx, Vec2f],
  sharp_mask: var Field[Ny, Nx, float64],
  kinematic_viscosity: float64,
  permeability: float64,
  cfl: float64,
  smoothing: float64,
): Simulation[Nx, Ny] =
  var fourier_vorticity = createField[Nx, Ny, Complex64](complex(0.0))
  for i in low(FieldIndex[Nx, Ny])..high(FieldIndex[Nx, Ny]):
    let idx: FieldIndex[Nx, Ny] = i
    fourier_vorticity[i] = complex(vorticity[idx])
  fft fourier_vorticity 
  dealias fourier_vorticity

  result = Simulation[Nx, Ny](
    fourier_vorticity: fourier_vorticity,
    object_velocity: object_velocity,
    mask: smooth_gaussian(sharp_mask, smoothing),
    nonlinear: create_field[Ny, Nx, Complex64](complex(0.0)),
    kinematic_viscosity: kinematic_viscosity,
    cfl: cfl,
    permeability: permeability,
    Δt: 0.0,
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
  object_velocity: Field[Nx, Ny, Vec2f],
  mask: Field[Nx, Ny, float64],
  permeability: float64,
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

  var ∂gx∂y = create_field[Nx, Ny, Complex64](complex(0.0))
  var ∂gy∂x = create_field[Nx, Ny, Complex64](complex(0.0))
  for x, y in mask.indices:
    let m = mask[x, y]
    let v = velocity[x, y]
    let ov = object_velocity[x, y]
    ∂gx∂y[x, y] = complex(1.0 * m / permeability * (v.x - ov.x))
    ∂gy∂x[x, y] = complex(1.0 * m / permeability * (v.y - ov.y))

  fft ∂gx∂y
  fft ∂gy∂x

  dealias ∂gx∂y
  dealias ∂gy∂x

  for x, y in mask.indices:
    ∂gx∂y[x, y] = ∂gx∂y[x, y] * im(1.0) * k(y, Ny)
    ∂gy∂x[x, y] = ∂gy∂x[x, y] * im(1.0) * k(x, Nx)

  for x, y in nonlinear.indices:
    let v = velocity[x, y]
    assert(∇ω_x[x, y].im.abs < tolerance, msg = $∇ω_x[x, y].im.abs)
    assert(∇ω_y[x, y].im.abs < tolerance, msg = $∇ω_y[x, y].im.abs)
    nonlinear[x, y] = complex(-(v.x * ∇ω_x[x, y].re + v.y * ∇ω_y[x, y].re))
  fft nonlinear
  for x, y in nonlinear.indices:
    nonlinear[x,y] = nonlinear[x,y] + (∂gx∂y[x,y] - ∂gy∂x[x,y])
  dealias nonlinear

  result = nonlinear

proc step*[Nx, Ny: static int](simulation: Simulation[Nx, Ny]) =
  let time = simulation.time
  let Δx = 2 * PI / Nx
  let Δt_old = simulation.Δt
  let cfl = simulation.cfl
  let permeability = simulation.permeability
  let kinematic_viscosity = simulation.kinematic_viscosity
  var fourier_ω = simulation.fourier_vorticity
  let velocity = calculate_velocity simulation.fourier_vorticity
  let nonlinear_old = simulation.nonlinear
  let nonlinear_new = calculate_nonlinear(
    velocity, 
    simulation.fourier_vorticity,
    simulation.object_velocity,
    simulation.mask,
    simulation.permeability,
  )

  var max_velocity = 0.0
  for v in velocity.data: 
    if v.abs > max_velocity: 
      max_velocity = v.abs
  let Δt_new = 
    if cfl * Δx / max_velocity > permeability: permeability
    else: cfl * Δx / max_velocity

  if time == 0.0:
    for x, y in fourier_ω.indices:
      let k = (
        x: k(x, Nx),
        y: k(y, Ny),
      )

      if k.x == 0 and k.y == 0:
        continue

      let E = exp(-kinematic_viscosity * k.abs2 * Δt_new)
      let ω = fourier_ω[x, y]
      let B = (1.0 - E) / (kinematic_viscosity * k.abs2)
      let N = nonlinear_new[x, y]

      fourier_ω[x, y] = E * ω + B * N
  else:
    let β_10 = 0.5 * Δt_new / Δt_old * (Δt_new + 2 * Δt_old)
    let β_11 = -0.5 * Δt_new * Δt_new / Δt_old

    for x, y in fourier_ω.indices:
      let k = (
        x: k(x, Nx),
        y: k(y, Ny),
      )

      if k.x == 0 and k.y == 0:
        continue

      let E_new = exp(-kinematic_viscosity * Δt_new * k.abs2)
      let E_old = exp(-kinematic_viscosity * Δt_old * k.abs2)
      let ω = fourier_ω[x, y]
      let N_new = nonlinear_new[x, y]
      let N_old = nonlinear_old[x, y]

      fourier_ω[x, y] = E_new * (ω + β_10 * N_new + β_11 * E_old * N_old)

  dealias fourier_ω

  var test_ω = create_field[Nx, Ny, Complex64](complex(0.0))
  for x, y in fourier_ω.indices:
    test_ω[x, y] = fourier_ω[x, y]

  ifft test_ω

  for c in test_ω.data:
    assert(c.im.abs < tolerance, msg = $c.im.abs)

  simulation.nonlinear = nonlinear_new
  simulation.Δt = Δt_new
  simulation.time += Δt_new
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

  var image = newImage(Nx, Ny)
  for x, y in ω.indices:
    let re = ω[x, y].re
    let a = re / max_abs
    let c = 
      if a >= 0: rgba(uint8(a * 255.0), 0, 0, 255)
      else: rgba(0, 0, uint8(-a * 255.0), 255)
    image.data[y * Nx + x] = c

  result = image
