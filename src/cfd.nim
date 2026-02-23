import std/times
import std/monotimes
import std/complex
import std/math
import pixie
import field
import simulation

proc main() =
  const kinematic_velocity = 1e-2
  const cfl = 0.1
  const max_time = 5.0

  const nx = 128
  const ny = 128

  var vorticity = create_field[nx, ny, float64](0.0)

  for x, y in vorticity.indices:
    vorticity[x, y] = -2 * sin(x / nx * 2 * PI) * sin(y / ny * 2 * PI)

  let simulation = create_simulation(vorticity, kinematic_velocity, cfl)

  simulation.visualize.write_file("initial.png")

  let start = get_mono_time()

  while simulation.time < max_time:
    simulation.step()

  let finish = get_mono_time()

  echo in_milliseconds(finish - start), "ms"

  simulation.visualize.write_file("final.png")

main()
