import std/times
import std/monotimes
import std/complex
import std/math
import pixie
import field
import simulation
import linalg

proc circle_test[Nx, Ny: static int](
  mask: var Field[Nx, Ny, float64],
  time: float64,
  origin: Vec2f,
  v: Vec2f,
  radius: float64,
  smoothing: float64
): Field[Nx, Ny, float64] =
  let d = (x: 2.0 * PI / Nx, y: 2.0 * PI / Ny)
  let current_origin = (v ⊙/ d * time) + origin
  for x, y in mask.indices:
    var f = (x: float64(x), y: float64(y)) - current_origin
    f.x -= float64(Nx) * round(f.x / float64(Nx))
    f.y -= float64(Ny) * round(f.y / float64(Ny))
    
    mask[x, y] =
      if f.abs2 <= radius * radius: 1.0
      else: 0.0

  result = smooth_gaussian(mask, smoothing)

proc main() =
  const permeability = 1e-3
  const kinematic_velocity = 1e-2
  const cfl = 0.1
  const smoothing = 16.0
  const max_time = 1

  const nx = 128
  const ny = 128

  const origin = (x: 64.0, y: 64.0)
  const v = (x: 1.0, y: 0.0)
  const r = 10

  var vorticity = create_field[nx, ny, float64](0.0)
  for x, y in vorticity.indices:
    vorticity[x, y] = -2 * sin(x / nx * 2 * PI) * sin(y / ny * 2 * PI)
  var mask = create_field[nx, ny, float64](0.0)
  var real_mask = circle_test(
    mask, 0.0, origin, v,
    r, 0.0
  )
  let object_velocity = create_field[nx, ny, Vec2f](v)
  let simulation = create_simulation[nx, ny](
    vorticity,
    object_velocity,
    real_mask,
    kinematic_velocity,
    permeability,
    cfl,
    smoothing,
  )

  simulation.visualize.write_file("images/initial.png")

  let start = get_mono_time()

  var steps = 0
  while simulation.time < max_time:
    simulation.mask = circle_test(simulation.mask, simulation.time, origin, v, r, smoothing)
    simulation.step()
    inc steps

  let finish = get_mono_time()

  echo in_milliseconds(finish - start), "ms"
  echo steps, " steps"

  simulation.visualize.write_file("images/final.png")

main()
