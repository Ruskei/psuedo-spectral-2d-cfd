import std/sequtils
type
  Field*[Nx, Ny: static int, T] = object
    data*: seq[T]
  FieldIndex*[Nx, Ny: static int] = range[0..(Nx * Ny - 1)]

proc `[]`*[Nx, Ny: static int, T](
  field: Field[Nx, Ny, T],
  i: FieldIndex[Nx, Ny],
): T =
  result = field.data[i]

proc `[]=`*[Nx, Ny: static int, T](
  field: var Field[Nx, Ny, T],
  i: FieldIndex[Nx, Ny],
  value: T,
) =
  field.data[i] = value

proc `[]`*[Nx, Ny: static int, T](
  field: Field[Nx, Ny, T],
  x: int,
  y: int,
): T =
  result = field.data[y * Nx + x]

proc `[]=`*[Nx, Ny: static int, T](
  field: var Field[Nx, Ny, T],
  x: int,
  y: int,
  value: T
) =
  field.data[y * Nx + x] = value

proc create_field*[Nx, Ny: static int, T](initial: T): Field[Nx, Ny, T] =
  result = Field[Nx, Ny, T](data: new_seq_with(Nx * Ny, initial))

iterator indices*[Nx, Ny: static int; T](field: Field[Nx, Ny, T]): (int, int) {.inline.} =
  for y in 0..(Ny - 1):
    for x in 0..(Nx - 1):
      yield (x, y)
