import std/complex
from math import exp, PI

import pixie

import fourier
import field
import linalg
import math_util

const tolerance = 1e-9

type
  SimulationConstants* = object
    cfl*: float64
    kinematic_viscosity*: float64
    permeability*: float64
    smoothing*: float64
  Simulation*[Nx, Ny: static int] = ref object # little reason to ever copy Simulation
    Fω: Field[Nx, Ny, Complex64] # fourier vorticity
    N: Field[Nx, Ny, Complex64] # nonlinear term
    χ*: Field[Nx, Ny, float64] # object mask
    χ_v: Field[Nx, Ny, Vec2D] # object velocity
    t*: float64
    Δt: float64
    extents: Vec2D
    constants: SimulationConstants

    v: Field[Nx, Ny, Vec2D]
    t_v_x: Field[Nx, Ny, Complex64]
    t_v_y: Field[Nx, Ny, Complex64]
    t_N: Field[Nx, Ny, Complex64]
    t_∇ω_x: Field[Nx, Ny, Complex64]
    t_∇ω_y: Field[Nx, Ny, Complex64]
    t_∂gx∂y: Field[Nx, Ny, Complex64]
    t_∂gy∂x: Field[Nx, Ny, Complex64]
    N_old: Field[Nx, Ny, Complex64] 

proc mode(i, N: int): float64 =
  if i < N div 2: float64(i)
  elif i == N div 2: 0.0
  else: float64(i - N)

proc κ(i, N: int; L: float64): float64 = (2.0 * PI / L) * mode(i, N)

proc create_simulation*[Nx, Ny: static int](
  ω: Field[Nx, Ny, float64],
  χ: Field[Nx, Ny, float64],
  χ_v: Field[Nx, Ny, Vec2D],
  extents: Vec2D,
  constants: SimulationConstants,
): Simulation[Nx, Ny] =
  var Fω = create_field[Nx, Ny, Complex64](complex(0.0))
  for x, y in ω.indices:
    Fω[x, y] = complex(ω[x, y])
  fft Fω
  dealias Fω

  result = Simulation[Nx, Ny](
    Fω: Fω,
    N: create_field[Nx, Ny, Complex64](complex(0.0)),
    χ: χ,
    χ_v: χ_v,
    t: 0.0,
    Δt: 0.0,
    extents: extents,
    constants: constants,

    v: create_field[Nx, Ny, Vec2D]((0.0, 0.0)),
    t_v_x: create_field[Nx, Ny, Complex64](complex(0.0)),
    t_v_y: create_field[Nx, Ny, Complex64](complex(0.0)),
    t_N: create_field[Nx, Ny, Complex64](complex(0.0)),
    t_∇ω_x: create_field[Nx, Ny, Complex64](complex(0.0)),
    t_∇ω_y: create_field[Nx, Ny, Complex64](complex(0.0)),
    t_∂gx∂y: create_field[Nx, Ny, Complex64](complex(0.0)),
    t_∂gy∂x: create_field[Nx, Ny, Complex64](complex(0.0)),
    N_old: create_field[Nx, Ny, Complex64](complex(0.0)),
  )

proc smooth_gaussian*[Nx, Ny: static int](
  χ: var Field[Nx, Ny, float64],
  smoothing: float64,
  extents: Vec2D,
) =
  var Fχ = create_field[Nx, Ny, Complex64](complex(0.0))
  for x, y in χ.indices:
    Fχ[x, y] = complex(χ[x, y])
  fft Fχ
  for x, y in Fχ.indices:
    let kx = κ(x, Nx, extents.x)
    let ky = κ(y, Ny, extents.y)
    Fχ[x, y] *= exp(-smoothing * (kx * kx / Nx / Nx + ky * ky / Ny / Ny))
  ifft Fχ
  for x, y in χ.indices:
    χ[x, y] = Fχ[x, y].re

proc step*[Nx, Ny: static int](
  sim: var Simulation[Nx, Ny],
) =
  let μ = sim.constants.permeability
  let L = (x: sim.extents.x, y: sim.extents.y)
  
  # calculate velocity for this step
  for x, y in sim.Fω.indices:
    let κ = (x: κ(x, Nx, sim.extents.x), y: κ(y, Ny, sim.extents.y))
    if κ.x == 0.0 and κ.y == 0.0: continue
    let Fω_c = sim.Fω[x, y]
    sim.t_v_x[x, y] = im(1.0) * κ.y / κ.abs2 * Fω_c
    sim.t_v_y[x, y] = im(1.0) * -κ.x / κ.abs2 * Fω_c
  ifft sim.t_v_x
  ifft sim.t_v_y
  for x, y in sim.v.indices:
    assert sim.t_v_x[x, y].im.abs < tolerance
    assert sim.t_v_y[x, y].im.abs < tolerance
    sim.v[x, y] = (
      sim.t_v_x[x, y].re,
      sim.t_v_y[x, y].re,
    )

  # calculate nonlinear term
  for x, y in sim.Fω.indices:
    let Fω_c = sim.Fω[x, y]
    sim.t_∇ω_x[x, y] = Fω_c * im(1.0) * κ(x, Nx, L.x)
    sim.t_∇ω_y[x, y] = Fω_c * im(1.0) * κ(y, Ny, L.y)
  ifft sim.t_∇ω_x
  ifft sim.t_∇ω_y

  for x, y in sim.χ.indices:
    let v_c = sim.v[x, y]
    sim.N[x, y] = complex(-(v_c.x * sim.t_∇ω_x[x, y].re + v_c.y * sim.t_∇ω_y[x, y].re))
  fft sim.N

  for x, y in sim.χ.indices:
    let χ_c = sim.χ[x, y]
    let v_c = sim.v[x, y]
    let χ_v_c = sim.χ_v[x, y]
    sim.t_∂gx∂y[x, y] = complex(1.0 * χ_c / μ * (v_c.x - χ_v_c.x))
    sim.t_∂gy∂x[x, y] = complex(1.0 * χ_c / μ * (v_c.y - χ_v_c.y))
  fft sim.t_∂gx∂y
  fft sim.t_∂gy∂x

  dealias sim.t_∂gx∂y
  dealias sim.t_∂gy∂x
    
  for x, y in sim.N.indices:
    sim.N[x, y] += im(1.0) * (sim.t_∂gx∂y[x,y] * κ(y, Ny, L.y) - sim.t_∂gy∂x[x,y] * κ(x, Nx, L.x))

  dealias sim.N

  # do step
  let cfl = sim.constants.cfl
  let ν = sim.constants.kinematic_viscosity
  let Δx = min(L.x / Nx, L.y / Ny)
  var max_v = 1e-12
  for v in sim.v.data:
    if v.abs > max_v: max_v = v.abs
  let Δt_new =
    if cfl * Δx / max_v > μ: μ
    else: cfl * Δx / max_v 
  
  if sim.t == 0.0:
    for x, y in sim.Fω.indices:
      let κ = (
        x: κ(x, Nx, L.x),
        y: κ(y, Ny, L.y),
      )

      if κ.x == 0.0 and κ.y == 0.0: continue
      
      let E = exp(-ν * κ.abs2 * Δt_new)
      let Fω_c = sim.Fω[x, y]
      let B = (1.0 - E) / (ν * κ.abs2)
      let N = sim.N[x, y]

      sim.Fω[x, y] = E * Fω_c + B * N
  else:
    let Δt_old = sim.Δt
    let β_10 = 0.5 * Δt_new / Δt_old * (Δt_new + 2 * Δt_old)
    let β_11 = -0.5 * Δt_new * Δt_new / Δt_old

    for x, y in sim.Fω.indices:
      let κ = (
        x: κ(x, Nx, L.x),
        y: κ(y, Ny, L.y),
      )

      if κ.x == 0.0 and κ.y == 0.0: continue

      let E_new = exp(-ν * Δt_new * κ.abs2)
      let E_old = exp(-ν * Δt_old * κ.abs2)
      let Fω_c = sim.Fω[x, y]
      let N_new = sim.N[x, y]
      let N_old = sim.N_old[x, y]

      sim.Fω[x, y] = E_new * (Fω_c + β_10 * N_new + β_11 * E_old * N_old)

  swap(sim.N_old, sim.N)
  sim.Δt = Δt_new
  sim.t += sim.Δt

type
  VorticityView* = enum
    vvLinear,
    vvSymLog

proc visualize_vorticity*[Nx, Ny: static int](
  sim: Simulation[Nx, Ny],
  view: VorticityView = vvLinear,
): Image =
  var ω = create_field[Nx, Ny, Complex64](complex(0.0))
  for x, y in sim.Fω.indices:
    ω[x, y] = sim.Fω[x, y]
  ifft ω

  for c in ω.data:
    assert c.im.abs < tolerance

  var maxAbs = 0.0
  for c in ω.data:
    let a = abs(c.re)
    if a > maxAbs: maxAbs = a

  var image = newImage(Nx, Ny)

  if maxAbs <= 0.0:
    for i in 0 ..< image.data.len:
      image.data[i] = rgba(245'u8, 245'u8, 245'u8, 255'u8)
    return image

  let linThresh = 0.03
  let symlogDen = 1.0 + ln(1.0 / linThresh)

  for x, y in ω.indices:
    var a = ω[x, y].re / maxAbs

    if a > 1.0: a = 1.0
    elif a < -1.0: a = -1.0

    var s = a
    if view == vvSymLog:
      let aa = abs(a)
      if aa <= linThresh:
        s = (a / linThresh) * (1.0 / symlogDen)
      else:
        let signA = (if a < 0.0: -1.0 else: 1.0)
        s = signA * ((1.0 + ln(aa / linThresh)) / symlogDen)

    let t = 0.5 * (s + 1.0)

    var r, g, b: float64
    if t <= 0.5:
      let u = t / 0.5
      let w = u * u * (3.0 - 2.0 * u)
      r = 49.0  + (247.0 - 49.0)  * w
      g = 54.0  + (247.0 - 54.0)  * w
      b = 149.0 + (247.0 - 149.0) * w
    else:
      let u = (t - 0.5) / 0.5
      let w = u * u * (3.0 - 2.0 * u)
      r = 247.0 + (165.0 - 247.0) * w
      g = 247.0 + (0.0   - 247.0) * w
      b = 247.0 + (38.0  - 247.0) * w

    if r < 0.0: r = 0.0 elif r > 255.0: r = 255.0
    if g < 0.0: g = 0.0 elif g > 255.0: g = 255.0
    if b < 0.0: b = 0.0 elif b > 255.0: b = 255.0

    image.data[y * Nx + x] = rgba(uint8(r), uint8(g), uint8(b), 255'u8)

  result = image
