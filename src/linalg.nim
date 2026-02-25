import std/complex
import math

type
  Vec2*[T] = tuple[x: T, y: T]
  Vec2C* = Vec2[Complex64]
  Vec2D* = Vec2[float64]

#vec vec ops
#immutable
proc `∙`*[T](a, b: Vec2[T]): T =
  result = a.x * b.x + a.y * b.y
proc `+`*[T](a, b: Vec2[T]): Vec2[T] =
  result.x = a.x + b.x
  result.y = a.y + b.y
proc `-`*[T](a, b: Vec2[T]): Vec2[T] =
  result.x = a.x - b.x
  result.y = a.y - b.y
#mutators
proc `+=`*[T](a: var Vec2[T], b: Vec2[T]) =
  a.x += b.x
  a.y += b.y
proc `⊙`*[T](a: Vec2[T], b: Vec2[T]): Vec2[T] =
  result.x = a.x * b.x
  result.y = a.y * b.y

#vec scalar ops
proc `*`*[T](a: Vec2[T], b: T): Vec2[T] =
  result.x = a.x * b
  result.y = a.y * b
proc `/`*[T](a: Vec2[T], b: T): Vec2[T] =
  result.x = a.x / b
  result.y = a.y / b

#unitary ops
proc abs*[T](a: Vec2[T]): T =
  result = sqrt(a.x * a.x + a.y * a.y)
proc abs2*[T](a: Vec2[T]): T =
  result = a.x * a.x + a.y * a.y

converter to_vec2d*(v: Vec2[int]): Vec2D =
  result = (float64(v.x), float64(v.y))
