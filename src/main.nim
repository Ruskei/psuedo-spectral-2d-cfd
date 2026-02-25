import std/times
import std/monotimes
import math

import pixie

import sim
import field
import linalg

proc set_circle_χ[Nx, Ny: static int](
  χ: var Field[Nx, Ny, float64],
  extents: Vec2D,
  t: float64,
  o: Vec2D,
  v: Vec2D,
  r: float64,
  smoothing: float64,
) =
  # convert from real-space to node-space
  let f = (x: Nx / extents.x, y: Ny / extents.y)
  let p = (o + v * t) ⊙ f
  let rr = (r * Nx / extents.x) * (r * Nx / extents.x)

  for x, y in χ.indices:
    var d = (float64(x), float64(y)) - p
    d.x -= round(d.x / Nx) * Nx
    d.y -= round(d.y / Ny) * Ny

    χ[x, y] =
      if d.abs2 <= rr: 1.0
      else: 0.0

  smooth_gaussian(χ, smoothing, extents)

proc main() =
  const max_time = 0.8
  const Nx = 128
  const Ny = 64
  const extents = (8, 4)
  const constants = SimulationConstants(
    cfl: 0.1,
    kinematic_viscosity: 1e-2,
    permeability: 1e-3,
    smoothing: 16.0,
  )

  const o = (1, 2)
  const v = (8, 0)
  const r = 0.8
  const stride = 8

  let startup = get_mono_time()

  var ω = create_field[Nx, Ny, float64](0.0)

  var χ = create_field[Nx, Ny, float64](0.0)
  set_circle_χ(χ, extents, 0.0, o, v, r, constants.smoothing)

  var χ_v = create_field[Nx, Ny, Vec2D](v)

  var sim = create_simulation(ω, χ, χ_v, extents, constants)

  sim.visualize_vorticity.write_file("images/initial.png")

  let finished_startup = get_mono_time()

  var snapshots_taken = 0
  var snapshot_interval = 0.2

  while sim.t < max_time:
    set_circle_χ(sim.χ, extents, sim.t, o, v, r, constants.smoothing)
    sim.step()

    if sim.t > float64(snapshots_taken + 1) * snapshot_interval:
      inc snapshots_taken
      var progress_image = sim.visualize_vorticity
      draw_streamline_overlay(progress_image, sim, 11, qcm_black)
      progress_image.write_file("images/progress" & $sim.t & ".png")

  let finished_sim = get_mono_time()

  var final_image = sim.visualize_vorticity
  draw_quiver_overlay(final_image, sim, stride, qcm_black)
  final_image.write_file("images/final.png")

  echo "startup took ", in_milliseconds(finished_startup - startup), "ms"
  echo "stepping took ", in_milliseconds(finished_sim - finished_startup), "ms"
main()
