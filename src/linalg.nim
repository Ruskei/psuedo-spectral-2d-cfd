import std/complex
import math

type
  Vec2*[T] = tuple[x: T, y: T]
  Vec2f* = Vec2[float64]
  Vec2C* = Vec2[Complex64]

proc `∙`*[T](a, b: Vec2[T]): T =
  result = a.x * b.x + a.y * b.y

proc `+=`*[T](a: var Vec2[T], b: Vec2[T]) =
  a.x += b.x
  a.y += b.y

proc `+`*[T](a, b: Vec2[T]): Vec2[T] =
  result.x = a.x + b.x
  result.y = a.y + b.y

proc `-`*[T](a, b: Vec2[T]): Vec2[T] =
  result.x = a.x - b.x
  result.y = a.y - b.y

proc `*`*[T](a: Vec2[T], b: T): Vec2[T] =
  result.x = a.x * b
  result.y = a.y * b

proc `/`*[T](a: Vec2[T], b: T): Vec2[T] =
  result.x = a.x / b
  result.y = a.y / b

proc `⊙/`*[T](a: Vec2[T], b: Vec2[T]): Vec2[T] =
  result.x = a.x / b.x
  result.y = a.y / b.y

proc abs*[T](a: Vec2[T]): T =
  result = sqrt(a.x * a.x + a.y * a.y)

proc abs2*[T](a: Vec2[T]): T =
  result = a.x * a.x + a.y * a.y
