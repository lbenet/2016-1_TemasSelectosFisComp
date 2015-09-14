module Autdiff.jl

export Dual, dual_var

type Dual{T<:Number}
    x :: T
    y :: T
end

Dual(x,y) = Dual(promote(x,y)...)

dual_var(x0::Number) = Dual(x0,1.0)

import Base: +, -, * , / , ^
+(a::Dual, b::Dual) = Dual( a.x + b.x, a.y + b.y )
+(a::Dual, b::Real) = Dual( a.x + b, a.y)
+(b::Real, a::Dual) = Dual( a.x + b, a.y)
+(a::Dual, b::Complex) = Dual( a.x + b, a.y)
+(b::Complex, a::Dual) = Dual( a.x + b, a.y)
+(a::Dual) = Dual(+a.x,+a.y)

-(a::Dual, b::Dual) = Dual( a.x - b.x, a.y - b.y )
-(a::Dual, b::Real) = Dual( a.x - b, a.y)
-(b::Real, a::Dual) = Dual( a.x - b, a.y)
-(a::Dual, b::Complex) = Dual( a.x - b, a.y)
-(b::Complex, a::Dual) = Dual( a.x - b, a.y)
-(a::Dual) = Dual(-a.x,-a.y)


*(a::Dual, b::Dual) = Dual( a.x * b.x, a.y * b.x + b.y * a.x)
*(a::Dual, b::Real) = Dual( a.x * b, a.y * b)
*(b::Real, a::Dual) = Dual( a.x * b, a.y * b)
*(a::Dual, b::Complex) = Dual( a.x * b, a.y * b)
*(b::Complex, a::Dual) = Dual( a.x * b, a.y * b)

/(a::Dual, b::Dual) = Dual( a.x / b.x, ( a.y - (a.x/b.x)b.y)/(b.x))
/(a::Dual, b::Real) = Dual( a.x / b, a.y/b )
/(b::Real, a::Dual) = Dual( a.x / b, a.y/b )
/(a::Dual, b::Complex) = Dual( a.x / b, a.y/b )
/(b::Complex, a::Dual) = Dual( a.x / b, a.y/b )

^(a::Dual, b::Dual) = Dual( a.x ^b.x, b.x * a.x^(b.x-1) * a.y)
^(a::Dual, b::Integer) = Dual( a.x ^ b, b * a.x^(b-1) * a.y)
^(a::Dual, b::Real) = Dual( a.x ^ b, b * a.x^(b-1) * a.y)

import Base: abs, sin, cos, tan, cot, sec, sec, asin, acos, atan, acot,
asec, acsc, exp, sqrt, log, sinh, cosh, tanh, coth, sech, csch, asinh,
acosh, atanh, acoth, asech, acsch

abs(a::Dual)=Dual(abs(a.x),a.y * sign(a.x)/sign(-a.x)) #con el fin que esté bien definida en 0.
sin(a::Dual)=Dual(sin(a.x),cos(a.x)*a.y)
cos(a::Dual)=Dual(cos(a.x),-sin(a.x)*a.y)
tan(a::Dual)=Dual(tan(a.x),sec(a.x)^2 * a.y)
cot(a::Dual)=Dual(cot(a.x),-csc(a.x)^2 * a.y)
sec(a::Dual)=Dual(sec(a.x),sec(a.x)*tan(a.x)*a.y)
csc(a::Dual)=Dual(csc(a.x),-csc(a.x)*cot(a.x)*a.y)
asin(a::Dual)=Dual(asin(a.x),a.y/sqrt(1-a.x^2))
acos(a::Dual)=Dual(acos(a.x),-a.y/sqrt(1-a.x^2))
atan(a::Dual)=Dual(atan(a.x),a.y/(a.x^2+1))
acot(a::Dual)=Dual(acot(a.x),-a.y/(1+a.x^2))
asec(a::Dual)=Dual(asec(a.x),a.y/(abs(a.x)*sqrt(a.x^2-1)))
acsc(a::Dual)=Dual(acsc(a.x),-a.y/(abs(a.x)*sqrt(a.x^2-1)))
exp(a::Dual)=Dual(exp(a.x),a.x*exp(a.x)*a.y)
sqrt(a::Dual)=Dual(sqrt(a.x),a.y / (2 * sqrt(a.x)))
log(a::Dual)=Dual(log(a.x),a.y/a.x)
sinh(a::Dual)=Dual(sinh(a.x),a.y * cosh(a.x))
cosh(a::Dual)=Dual(cosh(a.x),a.y * sinh(a.x))
tanh(a::Dual)=sinh(a)/cosh(a)
coth(a::Dual)=cosh(a)/sinh(a)
sech(a::Dual)=1/cosh(a)
csch(a::Dual)=1/sinh(a)
asinh(a::Dual)=Dual(asinh(a.x),a.y/sqrt(a.x^2+1))
acosh(a::Dual)=Dual(acosh(a.x),a.y/sqrt(a.x^2-1))
atanh(a::Dual)=Dual(atanh(a.x),a.y/(1-a.x^2))
acoth(a::Dual)=Dual(acoth(a.x),a.y/(1-a.x^2))
asech(a::Dual)=Dual(asech(a.x),-a.y/(a.x*sqrt(1-a.x^2)))
acsch(a::Dual)=Dual(acsch(a.x),-a.y/(abs(a.x)*sqrt(1+a.x^2)))

end
