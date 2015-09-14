##########################################################################

# En este módulo se busca construir la diferenciación automática para
# derivadas de primer orden.
# Este programa se ha realizado en base al libro:
# Valitade Numerics: A Short Introduction to Rigorous Computation,
# Warwick Tucker, Princenton University Press.

##########################################################################

module AutDiff # Empezamos por definir el módulo, esto se escribirá cuando se haga referencia a este módulo en un notebook.

export Dual, dual_var, dual_cons# El módulo debe reconocer tres diferentes entradas: Duales, variables y constantes. 

type Dual{T<:Real}# Lo primero es definir los duales, los cuales serán "arreglos" con dos variables, se nos piden que sean de un subgrupo de reales, de ahí el uso de <: Real.
    x :: T
    y :: T
end

Dual(x,y) = Dual(promote(x,y)...)# Usamos promote para que ambos campos tengan el mismo tipo.

# Se nos pide definir una función dual_var que devuelva la variable independiente en un dual, ésto es sencillo.

function dual_var(x0)
    return Dual(x0,1.0)
end

function dual_cons(x0)
    return Dual(x0,0.0)
end

# Se nos pide definir métodos para que el dual de un sólo número sea lo que esperamos, ésto es sencillo usando la función dual_var:

Dual(a) = dual_cons(a)

# Paso siguiente, el programa debe realizar suma, resta, producto y división de duales, ésto se presenta a continuación.
# Sin embargo, se debe ser honesto en este sentido y decir que esta parte se copió del trabajo de César Bertoni, pues su notebook corrige los errores que tenía el mío (no estaba bien definido: -Dual(a,b), marcaba error).

import Base.+
+(a::Dual) = Dual(+a.x,+a.y)
+(a::Dual, b::Dual) = Dual( a.x + b.x, a.y + b.y )
+(a::Dual, b::Real) = a + Dual(b)
+(b::Real, a::Dual) = a + Dual(b)

import Base.-
-(a::Dual) = Dual(-a.x,-a.y)
-(a::Dual, b::Dual) = Dual( a.x - b.x, a.y - b.y )
-(a::Dual, b::Real) = a - Dual(b)
-(b::Real, a::Dual) = a - Dual(b)

import Base.*
*(a::Dual, b::Dual) = Dual( a.x * b.x, a.y * b.x + b.y * a.x)
*(a::Dual, b::Real) = a*Dual(b)
*(b::Real, a::Dual) = a*Dual(b)

import Base./
/(a::Dual, b::Dual) = Dual( a.x / b.x, ( a.y - (a.x/b.x)b.y)/(b.x))
/(a::Dual, b::Real) = a/Dual(b)
/(b::Real, a::Dual) = Dual(b)/a

import Base.^
^(a::Dual, b::Dual) = Dual( a.x ^b.x, b.x * a.x^(b.x-1) * a.y)
^(a::Dual, b::Integer) = Dual( a.x ^ b, b * a.x^(b-1) * a.y)
^(a::Dual, b::Real) = Dual( a.x ^ b, b * a.x^(b-1) * a.y)

# Lo siguiente es definir las operaciones básicas para duales.

import Base.abs
abs(a::Dual)=Dual(abs(a.x),a.y * sign(a.x)/sign(-a.x))

import Base.sin
sin(a::Dual)=Dual(sin(a.x),cos(a.x)*a.y)

import Base.cos
cos(a::Dual)=Dual(cos(a.x),-sin(a.x)*a.y)

import Base.tan
tan(a::Dual)=Dual(tan(a.x),sec(a.x)^2 * a.y)

import Base.cot
cot(a::Dual)=Dual(cot(a.x),-csc(a.x)^2 * a.y)

import Base.sec
sec(a::Dual)=Dual(sec(a.x),sec(a.x)*tan(a.x)*a.y)

import Base.csc
csc(a::Dual)=Dual(csc(a.x),-csc(a.x)*cot(a.x)*a.y)

import Base.asin
asin(a::Dual)=Dual(asin(a.x),a.y/sqrt(1-a.x^2))

import Base.acos
acos(a::Dual)=Dual(acos(a.x),-a.y/sqrt(1-a.x^2))

import Base.atan
atan(a::Dual)=Dual(atan(a.x),a.y/(a.x^2+1))

import Base.acot
acot(a::Dual)=Dual(acot(a.x),-a.y/(1+a.x^2))

import Base.asec
asec(a::Dual)=Dual(asec(a.x),a.y/(abs(a.x)*sqrt(a.x^2-1)))

import Base.acsc
acsc(a::Dual)=Dual(acsc(a.x),-a.y/(abs(a.x)*sqrt(a.x^2-1)))

import Base.exp
exp(a::Dual)=Dual(exp(a.x),exp(a.x)*a.y)

import Base.sqrt
sqrt(a::Dual)=Dual(sqrt(a.x),a.y / (2 * sqrt(a.x)))

import Base.log
log(a::Dual)=Dual(log(a.x),a.y/a.x)

end