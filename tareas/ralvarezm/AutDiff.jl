###  Módulo que implementa números Duales      ###
###  para obtener la primera derivada de una   ###
###  función usando diferenciación automática. ###

module AutDiff
export Dual,dual_var
immutable Dual{T<:Number} <:Number
    f::T
    fp::T
end

Dual(f,fp)=Dual(promote(f,fp)...)
dual_var(x0::Number)=Dual(x0,1.0)


import Base:+,-,*,/,^

+(x::Dual,y::Complex) = Dual(x.f+y,x.fp)
+(x::Complex,y::Dual) = Dual(x+y.f,x.fp)
+(x::Dual,y::Dual)=Dual(x.f+y.f,x.fp+y.fp)
+(x::Dual,y::Real)=Dual(x.f+y,x.fp)
+(x::Real,y::Dual)=Dual(x+y.f,y.fp)
+(x::Dual) = Dual(+x.f,+x.fp)  

-(x::Dual,y::Complex) = Dual(x.f-y,x.fp)
-(x::Complex,y::Dual) = Dual(x-y.f,-y.fp)
-(x::Dual,y::Dual)=Dual(x.f-y.f,x.fp-y.fp)
-(x::Dual,y::Real)=Dual(x.f-y,x.fp)
-(x::Real,y::Dual)=Dual(x-y.f,-y.fp)
-(x::Dual)=Dual(-x.f,-x.fp)

*(x::Dual, y::Dual) = Dual(x.f*y.f,x.fp*y.f+y.fp*x.f)
*(x::Dual, y::Real) = Dual(x.f*y,x.fp*y)
*(x::Real, y::Dual) = Dual(x*y.f,x*y.fp)
*(x::Dual, y::Complex) = Dual(x.f*y,x.fp*y)
*(x::Complex, y::Dual) = Dual(x*y.f,x*y.fp)

/(x::Dual, y::Dual) = Dual(x.f/y.f,(x.fp-(x.f/y.f)y.fp)/(y.f))
/(x::Dual, y::Real) = Dual(x.f/y,x.fp/y)
/(x::Real, y::Dual) = Dual(x/y.f,(-x*y.fp)/(y.f)^2)
/(x::Dual, y::Complex) = Dual(x.f/y,x.fp/y)
/(x::Complex, y::Dual) = Dual(x/y.f,(-x*y.fp)/((y.f)^2))

^(x::Dual,y::Dual) = Dual(x.f^y.f,y.f*x.f^(y.f-1)*x.fp)
^(x::Dual,y::Integer) = Dual(x.f^y,y*x.f^(y-1)*x.fp)
^(x::Dual,y::Real) = Dual(x.f^y,y*x.f^(y-1)*x.fp)

import Base: abs,sin, cos, tan, cot, sec, sec, asin, acos, atan, acot, 
asec, acsc, exp, sqrt, log, sinh, cosh, tanh, coth, sech, csch, asinh, 
acosh, atanh, acoth, asech, acsch

abs(x::Dual)=Dual(abs(x.f),x.f*(sign(x.f)/sign(x.f)))  #Tomo prestada de la tarea de Ceboc esta original forma de definir el dual de abs.
sin(x::Dual)=Dual(sin(x.f),cos(x.f)*x.fp)
cos(x::Dual)=Dual(cos(x.f),-sin(x.f)*x.fp)
tan(x::Dual)=Dual(tan(x.f),sec(x.f)^2 * x.fp)
cot(x::Dual)=Dual(cot(x.f),-csc(x.f)^2 * x.fp)
sec(x::Dual)=Dual(sec(x.f),sec(x.f)*tan(x.f)*x.fp)
csc(x::Dual)=Dual(csc(x.f),-csc(x.f)*cot(x.f)*x.fp)
asin(x::Dual)=Dual(asin(x.f),x.fp/sqrt(1-x.f^2))
acos(x::Dual)=Dual(acos(x.f),-x.fp/sqrt(1-x.f^2))
atan(x::Dual)=Dual(atan(x.f),x.fp/(x.f^2+1))
acot(x::Dual)=Dual(acot(x.f),-x.fp/(x.f^2+1))
asec(x::Dual)=Dual(asec(x.f),x.fp/(abs(x.f)*sqrt(x.f^2-1)))
acsc(x::Dual)=Dual(acsc(x.f),-x.fp/(abs(x.f)*sqrt(x.f^2-1)))
exp(x::Dual)=Dual(exp(x.f),exp(x.f)*x.fp)
sqrt(x::Dual)=Dual(sqrt(x.f),x.fp/(2*sqrt(x.f)))
log(x::Dual)=Dual(log(x.f),x.fp/x.f)
sinh(x::Dual)=Dual(sinh(x.f),cosh(x.f)*x.fp)
cosh(x::Dual)=Dual(cosh(x.f),sinh(x.f)*x.fp)
tanh(x::Dual)=sinh(x)/cosh(x)
coth(x::Dual)=cosh(x)/sinh(x)
sech(x::Dual)=1/cosh(x)
csch(x::Dual)=1/sinh(x)
asinh(x::Dual)=Dual(asinh(x.f),x.fp/sqrt(x.f^2+1))
acosh(x::Dual)=Dual(acosh(x.f),x.fp/sqrt(x.f^2-1))
atanh(x::Dual)=Dual(atanh(x.f),x.fp/(1-x.f^2))
acoth(x::Dual)=Dual(acoth(x.f),x.fp/(1-x.f^2))
asech(x::Dual)=Dual(asech(x.f),-x.fp/(x.f*sqrt(1-x.f^2)))
acsch(x::Dual)=Dual(acsch(x.f),-x.fp/(abs(x.f)*sqrt(1+x.f^2)))

end
