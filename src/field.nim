import std/sequtils

type
  Field*[Nx, Ny: static int; T] = object
    data*: seq[T]

proc `[]`*[Nx, Ny: static int; T](
  f: Field[Nx, Ny, T];
  x, y: int;
): T =
  result = f.data[y * Nx + x]
proc `[]`*[Nx, Ny: static int, T](
  f: var Field[Nx, Ny, T];
  x, y: int;
): var T =
  result = f.data[y * Nx + x]

proc `[]=`*[Nx, Ny: static int, T](
  f: var Field[Nx, Ny, T];
  x, y: int;
  value: sink T,
) =
  f.data[y * Nx + x] = value

proc create_field*[Nx, Ny: static int; T](initial: T): Field[Nx, Ny, T] =
  result = Field[Nx, Ny, T](data: new_seq_with(Nx * Ny, initial))

iterator indices*[Nx, Ny: static int; T](field: Field[Nx, Ny, T]): (int, int) {.inline.} =
  for y in 0..(Ny - 1):
    for x in 0..(Nx - 1):
      yield (x, y)
