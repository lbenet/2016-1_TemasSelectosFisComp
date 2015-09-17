##########################################################################

# En este módulo se busca construir la diferenciación automática para
# derivadas de primer orden.
# Este programa se ha realizado en base al libro:
# Valitade Numerics: A Short Introduction to Rigorous Computation,
# Warwick Tucker, Princenton University Press.

##########################################################################

module AutDiffC # Empezamos por definir el módulo, esto se escribirá cuando se haga referencia a este módulo en un notebook.

export DualC, dual_varC, dual_consC# El módulo debe reconocer tres diferentes entradas: Duales, variables y constantes. 

type DualC{T<:Number}# Lo primero es definir los duales, los cuales serán "arreglos" con dos variables, en este caso se requiere emplear duales con números complejos, es decir usaremos un número cualquiera (real o complejo).
    x :: T
    y :: T
end

DualC(x,y) = DualC(promote(x,y)...)# Usamos promote para que ambos campos tengan el mismo tipo.

# Se nos pide definir una función dual_var que devuelva la variable independiente en un dual, ésto es sencillo.

function dual_varC(x0)
    return DualC(x0,1.0)
end

function dual_consC(x0)
    return DualC(x0,0.0)
end

# Se nos pide definir métodos para que el dual de un sólo número sea lo que esperamos, ésto es sencillo usando a dunción dual_var:

DualC(x) = dual_consC(x)

# Paso siguiente, el programa debe realizar suma, resta, producto y división de duales, ésto se presenta a continuación.
# Sin embargo, se debe ser honesto en este sentido y decir que esta parte se copió del trabajo de César Bertoni, pues su notebook corrige los errores que tenía el mío (no estaba bien definido: -Dual(a,b), marcaba error).

import Base.+
+(a::DualC) = DualC(+a.x,+a.y)
+(a::DualC, b::DualC) = DualC( a.x + b.x, a.y + b.y )
+(a::DualC, b::Real) = DualC(b) + a
+(a::DualC, b::Complex) = DualC(b) + a
+(b::Real, a::DualC) = DualC(b) + a
+(b::Complex, a::DualC) = DualC(b) + a

import Base.-
-(a::DualC) = DualC(-a.x,-a.y)
-(a::DualC, b::DualC) = DualC( a.x - b.x, a.y - b.y )
-(a::DualC, b::Real) = a - DualC(b)
-(a::DualC, b::Complex) = a - DualC(b)
-(b::Real, a::DualC) = DualC(b) - a
-(b::Complex, a::DualC) = DualC(b) - a

import Base.*
*(a::DualC, b::DualC) = DualC( a.x * b.x, a.y * b.x + b.y * a.x)
*(a::DualC, b::Real) = DualC(b)*a
*(a::DualC, b::Complex) = DualC(b)*a
*(b::Real, a::DualC) = DualC(b)*a
*(b::Complex, a::DualC) = DualC(b)*a

import Base./
/(a::DualC, b::DualC) = DualC( a.x / b.x, ( a.y - (a.x/b.x)b.y)/(b.x))
/(a::DualC, b::Real) = a/DualC(b)
/(a::DualC, b::Complex) = a/DualC(b)
/(b::Real, a::DualC) = DualC(b)/a
/(b::Complex, a::DualC) = DualC(b)/a

import Base.^
^(a::DualC, b::DualC) = DualC( a.x ^b.x, b.x * a.x^(b.x-1) * a.y)
^(a::DualC, b::Integer) = DualC( a.x ^ b, b * a.x^(b-1) * a.y)
^(a::DualC, b::Real) = DualC( a.x ^ b, b * a.x^(b-1) * a.y)

# Lo siguiente es definir las operaciones básicas para duales.

import Base.abs
abs(a::DualC)=DualC(abs(a.x),a.y * sign(a.x)/sign(-a.x))

import Base.sin
sin(a::DualC)=DualC(sin(a.x),cos(a.x)*a.y)

import Base.cos
cos(a::DualC)=DualC(cos(a.x),-sin(a.x)*a.y)

import Base.tan
tan(a::DualC)=DualC(tan(a.x),sec(a.x)^2 * a.y)

import Base.cot
cot(a::DualC)=DualC(cot(a.x),-csc(a.x)^2 * a.y)

import Base.sec
sec(a::DualC)=DualC(sec(a.x),sec(a.x)*tan(a.x)*a.y)

import Base.csc
csc(a::DualC)=DualC(csc(a.x),-csc(a.x)*cot(a.x)*a.y)

import Base.asin
asin(a::DualC)=DualC(asin(a.x),a.y/sqrt(1-a.x^2))

import Base.acos
acos(a::DualC)=DualC(acos(a.x),-a.y/sqrt(1-a.x^2))

import Base.atan
atan(a::DualC)=DualC(atan(a.x),a.y/(a.x^2+1))

import Base.acot
acot(a::DualC)=DualC(acot(a.x),-a.y/(1+a.x^2))

import Base.asec
asec(a::DualC)=DualC(asec(a.x),a.y/(abs(a.x)*sqrt(a.x^2-1)))

import Base.acsc
acsc(a::DualC)=DualC(acsc(a.x),-a.y/(abs(a.x)*sqrt(a.x^2-1)))

import Base.exp
exp(a::DualC)=DualC(exp(a.x),exp(a.x)*a.y)

import Base.sqrt
sqrt(a::DualC)=DualC(sqrt(a.x),a.y / (2 * sqrt(a.x)))

import Base.log
log(a::DualC)=DualC(log(a.x),a.y/a.x)

end