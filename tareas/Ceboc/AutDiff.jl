module Autdiff.jl

export Dual, dual_var, dual_cte

type Dual{T<:Number}
    x :: T
    y :: T
end

Dual(x,y) = Dual(promote(x,y)...)

dual_var(x0::Number) = Dual(x0,1.0)
dual_cte(x0::Number) = Dual(x0,0.0)

import Base.+
+(a::Dual, b::Dual) = Dual( a.x + b.x, a.y + b.y )
+(a::Dual, b::Real) = Dual( a.x + b, a.y)
+(b::Real, a::Dual) = Dual( a.x - b, a.y)
+(a::Dual, b::Complex) = Dual( a.x + b, a.y)
+(b::Complex, a::Dual) = Dual( a.x + b, a.y)
+(a::Dual) = Dual(+a.x,+a.y)

import Base.+
-(a::Dual, b::Dual) = Dual( a.x - b.x, a.y - b.y )
-(a::Dual, b::Real) = Dual( a.x - b, a.y)
-(b::Real, a::Dual) = Dual( a.x - b, a.y)
-(a::Dual, b::Complex) = Dual( a.x - b, a.y)
-(b::Complex, a::Dual) = Dual( a.x - b, a.y)
-(a::Dual) = Dual(-a.x,-a.y)

import Base.*
*(a::Dual, b::Dual) = Dual( a.x * b.x, a.y * b.x + b.y * a.x)
*(a::Dual, b::Real) = Dual( a.x * b, a.y * b)
*(b::Real, a::Dual) = Dual( a.x * b, a.y * b)
*(a::Dual, b::Complex) = Dual( a.x * b, a.y * b)
*(b::Complex, a::Dual) = Dual( a.x * b, a.y * b)

import Base./
/(a::Dual, b::Dual) = Dual( a.x / b.x, ( a.y - (a.x/b.x)b.y)/(b.x))
/(a::Dual, b::Real) = Dual( a.x / b, a.y/b )
/(b::Real, a::Dual) = Dual( a.x / b, a.y/b )
/(a::Dual, b::Complex) = Dual( a.x / b, a.y/b )
/(b::Complex, a::Dual) = Dual( a.x / b, a.y/b )

import Base.^
^(a::Dual, b::Dual) = Dual( a.x ^b.x, b.x * a.x^(b.x-1) * a.y)
^(a::Dual, b::Integer) = Dual( a.x ^ b, b * a.x^(b-1) * a.y)
^(a::Dual, b::Real) = Dual( a.x ^ b, b * a.x^(b-1) * a.y)

import Base.abs
abs(a::Dual)=Dual(abs(a.x),a.y * sign(a.x)/sign(-a.x)) #con el fin que esté bien definida en 0.

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
exp(a::Dual)=Dual(exp(a.x),a.x*exp(a.x)*a.y)

import Base.sqrt
sqrt(a::Dual)=Dual(sqrt(a.x),a.y / (2 * sqrt(a.x)))

import Base.log
log(a::Dual)=Dual(log(a.x),a.y/a.x)

import Base.sinh
sinh(a::Dual)=Dual(sinh(a.x),a.y * cosh(a.x))

import Base.cosh
cosh(a::Dual)=Dual(cosh(a.x),a.y * sinh(a.x))

import Base.tanh
tanh(a::Dual)=sinh(a)/cosh(a)

import Base.coth
coth(a::Dual)=cosh(a)/sinh(a)

import Base.sech
sech(a::Dual)=1/cosh(a)

import Base.csch
csch(a::Dual)=1/sinh(a)

import Base.asinh
asinh(a::Dual)=Dual(asinh(a.x),a.y/sqrt(a.x^2+1))

import Base.acosh
acosh(a::Dual)=Dual(acosh(a.x),a.y/sqrt(a.x^2-1))

import Base.atanh
atanh(a::Dual)=Dual(atanh(a.x),a.y/(1-a.x^2))

import Base.acoth
acoth(a::Dual)=Dual(acoth(a.x),a.y/(1-a.x^2))

import Base.asech
asech(a::Dual)=Dual(asech(a.x),-a.y/(a.x*sqrt(1-a.x^2)))

import Base.acsch
acsch(a::Dual)=Dual(acsch(a.x),-a.y/(abs(a.x)*sqrt(1+a.x^2)))

end
