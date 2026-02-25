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
    vv_linear,
    vv_sym_log

proc visualize_vorticity*[Nx, Ny: static int](
  sim: Simulation[Nx, Ny],
  view: VorticityView = vv_linear,
): Image =
  # internal helpers / constants so colors are not hardcoded inline
  proc clamp01_255(x: float64): float64 =
    if x < 0.0:
      0.0
    elif x > 255.0:
      255.0
    else:
      x

  proc smoothstep01(u: float64): float64 =
    u * u * (3.0 - 2.0 * u)

  type
    RgbColor = tuple[r, g, b: float64]

  let zero_field_color: RgbColor = (245.0, 245.0, 245.0)
  let neg_color: RgbColor = (49.0, 54.0, 149.0)
  let mid_color: RgbColor = (247.0, 247.0, 247.0)
  let pos_color: RgbColor = (165.0, 0.0, 38.0)

  var omega = create_field[Nx, Ny, Complex64](complex(0.0))
  for x, y in sim.Fω.indices:
    omega[x, y] = sim.Fω[x, y]
  ifft omega

  for c in omega.data:
    assert c.im.abs < tolerance

  var max_abs = 0.0
  for c in omega.data:
    let a = abs(c.re)
    if a > max_abs:
      max_abs = a

  var image = new_image(Nx, Ny)

  if max_abs <= 0.0:
    for i in 0 ..< image.data.len:
      image.data[i] = rgba(
        uint8(zero_field_color.r),
        uint8(zero_field_color.g),
        uint8(zero_field_color.b),
        255'u8
      )
    return image

  let lin_thresh = 0.03
  let symlog_den = 1.0 + ln(1.0 / lin_thresh)

  for x, y in omega.indices:
    var a = omega[x, y].re / max_abs

    if a > 1.0:
      a = 1.0
    elif a < -1.0:
      a = -1.0

    var s = a
    if view == vv_sym_log:
      let abs_a = abs(a)
      if abs_a <= lin_thresh:
        s = (a / lin_thresh) * (1.0 / symlog_den)
      else:
        let sign_a = (if a < 0.0: -1.0 else: 1.0)
        s = sign_a * ((1.0 + ln(abs_a / lin_thresh)) / symlog_den)

    let t = 0.5 * (s + 1.0)

    var r, g, b: float64
    if t <= 0.5:
      let u = t / 0.5
      let w = smoothstep01(u)
      r = neg_color.r + (mid_color.r - neg_color.r) * w
      g = neg_color.g + (mid_color.g - neg_color.g) * w
      b = neg_color.b + (mid_color.b - neg_color.b) * w
    else:
      let u = (t - 0.5) / 0.5
      let w = smoothstep01(u)
      r = mid_color.r + (pos_color.r - mid_color.r) * w
      g = mid_color.g + (pos_color.g - mid_color.g) * w
      b = mid_color.b + (pos_color.b - mid_color.b) * w

    r = clamp01_255(r)
    g = clamp01_255(g)
    b = clamp01_255(b)

    image.data[y * Nx + x] = rgba(uint8(r), uint8(g), uint8(b), 255'u8)

  result = image

type
  quiver_color_mode* = enum
    qcm_white,
    qcm_black

proc draw_quiver_overlay*[Nx, Ny: static int](
  image: var Image,
  sim: Simulation[Nx, Ny],
  stride: int,
  color_mode: quiver_color_mode = qcm_black,
) =
  ## Draws a simple quiver overlay (directed line segments) on top of an existing image.
  ## One line is drawn every `stride` cells.
  ##
  ## Direction is shown by drawing from the sample point in the direction of velocity.
  ## (No arrowhead; just a signed line segment.)
  ##
  ## Assumes:
  ## - sim.v is a physical-space velocity field (not Fourier-space)
  ## - image dimensions match Nx x Ny

  let min_stride = 1
  let sample_stride =
    if stride < min_stride: min_stride
    else: stride

  let line_length_fraction_of_stride = 0.9
  let min_visible_speed = 1e-12
  let min_drawn_length = 0.0
  let line_alpha = 255'u8

  let line_color =
    case color_mode
    of qcm_white: rgba(255'u8, 255'u8, 255'u8, line_alpha)
    of qcm_black: rgba(0'u8, 0'u8, 0'u8, line_alpha)

  # Compute global max speed for normalization.
  var max_speed = 0.0
  for v in sim.v.data:
    let speed = v.abs
    if speed > max_speed:
      max_speed = speed

  if max_speed <= min_visible_speed:
    return

  let max_line_length = line_length_fraction_of_stride * float64(sample_stride)

  var y = 0
  while y < Ny:
    var x = 0
    while x < Nx:
      let v = sim.v[x, y]
      let speed = v.abs

      if speed > min_visible_speed:
        let x_dir = v.x / speed
        let y_dir = v.y / speed

        var line_length = max_line_length * (speed / max_speed)
        if line_length < min_drawn_length:
          line_length = min_drawn_length

        # Tail is the sample point; tip points downstream (signed direction).
        let x_tail_f = float64(x)
        let y_tail_f = float64(y)
        let x_tip_f = x_tail_f + x_dir * line_length
        let y_tip_f = y_tail_f + y_dir * line_length

        var x0 = int(round(x_tail_f))
        var y0 = int(round(y_tail_f))
        let x1 = int(round(x_tip_f))
        let y1 = int(round(y_tip_f))

        # Bresenham line rasterization (single-pixel line).
        let dx_abs = abs(x1 - x0)
        let dy_abs = abs(y1 - y0)
        let x_step = (if x0 < x1: 1 else: -1)
        let y_step = (if y0 < y1: 1 else: -1)
        var error = dx_abs - dy_abs

        while true:
          if x0 >= 0 and x0 < Nx and y0 >= 0 and y0 < Ny:
            image.data[y0 * Nx + x0] = line_color

          if x0 == x1 and y0 == y1:
            break

          let double_error = 2 * error
          if double_error > -dy_abs:
            error -= dy_abs
            x0 += x_step
          if double_error < dx_abs:
            error += dx_abs
            y0 += y_step

      x += sample_stride
    y += sample_stride

proc draw_streamline_overlay*[Nx, Ny: static int](
  image: var Image,
  sim: Simulation[Nx, Ny],
  stride: int,
  color_mode: quiver_color_mode = qcm_black,
) =
  ## Draws simple streamline overlays on top of an existing image.
  ## Streamlines are seeded every `stride` cells and traced using nearest-neighbor
  ## velocity sampling with fixed pixel step length.
  ##
  ## Assumes:
  ## - sim.v is a physical-space velocity field (not Fourier-space)
  ## - image dimensions match Nx x Ny

  let min_stride = 1
  let sample_stride =
    if stride < min_stride: min_stride
    else: stride

  let min_visible_speed = 1e-12
  let step_size_pixels = 1.0              # integration step in pixels
  let max_stream_steps = 24               # per direction
  let line_alpha = 255'u8

  let line_color =
    case color_mode
    of qcm_white: rgba(255'u8, 255'u8, 255'u8, line_alpha)
    of qcm_black: rgba(0'u8, 0'u8, 0'u8, line_alpha)

  let x_min = 0.0
  let y_min = 0.0
  let x_max = float64(Nx - 1)
  let y_max = float64(Ny - 1)

  var seed_y = 0
  while seed_y < Ny:
    var seed_x = 0
    while seed_x < Nx:
      # Trace both directions from the seed to form a streamline through the point.
      for direction_index in 0 .. 1:
        let direction_sign =
          if direction_index == 0: -1.0
          else: 1.0

        var x_pos = float64(seed_x)
        var y_pos = float64(seed_y)

        var step_count = 0
        while step_count < max_stream_steps:
          # Stop if outside image bounds.
          if x_pos < x_min or x_pos > x_max or y_pos < y_min or y_pos > y_max:
            break

          let x_sample = int(round(x_pos))
          let y_sample = int(round(y_pos))

          if x_sample < 0 or x_sample >= Nx or y_sample < 0 or y_sample >= Ny:
            break

          let v = sim.v[x_sample, y_sample]
          let speed = v.abs
          if speed <= min_visible_speed:
            break

          let x_dir = v.x / speed
          let y_dir = v.y / speed

          let x_next = x_pos + direction_sign * step_size_pixels * x_dir
          let y_next = y_pos + direction_sign * step_size_pixels * y_dir

          var x0 = int(round(x_pos))
          var y0 = int(round(y_pos))
          let x1 = int(round(x_next))
          let y1 = int(round(y_next))

          # Skip zero-length segments caused by rounding.
          if x0 == x1 and y0 == y1:
            x_pos = x_next
            y_pos = y_next
            inc step_count
            continue

          # Pixel-perfect 1px Bresenham segment.
          # Important: do NOT draw the final endpoint of each segment, so adjacent
          # segments do not double-plot joints (which creates "extra pixels").
          let dx_abs = abs(x1 - x0)
          let dy_abs = abs(y1 - y0)
          let x_step = (if x0 < x1: 1 else: -1)
          let y_step = (if y0 < y1: 1 else: -1)
          var error = dx_abs - dy_abs

          while true:
            # Stop before drawing the final endpoint; next segment will start there.
            if x0 == x1 and y0 == y1:
              break

            if x0 >= 0 and x0 < Nx and y0 >= 0 and y0 < Ny:
              image.data[y0 * Nx + x0] = line_color

            let double_error = 2 * error
            if double_error > -dy_abs:
              error -= dy_abs
              x0 += x_step
            if double_error < dx_abs:
              error += dx_abs
              y0 += y_step

          # Advance integrator state.
          x_pos = x_next
          y_pos = y_next

          # Early out if we've left the image.
          if x_pos < x_min or x_pos > x_max or y_pos < y_min or y_pos > y_max:
            break

          inc step_count

      seed_x += sample_stride
    seed_y += sample_stride
