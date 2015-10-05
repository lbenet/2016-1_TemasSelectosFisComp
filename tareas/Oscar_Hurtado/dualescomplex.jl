module dualescomplex
export Dual,variable
type Dual {T<:Number}
    u0:: T #función  evaluada en x_0
    up:: T #derivada
end
Dual(a,b)=Dual(promote(a,b)...)
Dual(a)=Dual(a,0)

import Base: +,-,*,/,^
+(a::Dual, b::Dual) = Dual( a.u0 + b.u0, a.up+ b.up )
+(a::Dual,b::Number)=a+Dual(b)
+(a::Number,b::Dual)=b+a

-(a::Dual, b::Dual) = Dual( a.u0 - b.u0, a.up- b.up )
-(a::Dual,b::Number)=a-Dual(b)
-(a::Number,b::Dual)=b-a

*(a::Dual,b::Dual)=Dual(a.u0*b.u0,a.up*b.u0+a.u0*b.up)
*(a::Dual,b::Number)=a*Dual(b)
*(a::Number,b::Dual)=b*a

h(a::Dual,b::Dual)=a.u0/b.u0

/(a::Dual,b::Dual)=Dual(h(a,b),(a.up-h(a,b)*b.up)/b.u0)
/(a::Dual,b::Number)=a/Dual(b)
/(a::Number,b::Dual)=Dual(a)/b

^(a::Dual,b::Int64)=Dual(a.u0^b,b*a.u0^(b-1)*a.up)

variable(x)=Dual(x)+Dual(0,1)

import Base: exp,log,sin,cos,tan
exp(a::Dual)=Dual(exp(a.u0),a.up*exp(a.u0))
log(a::Dual)=Dual(log(a.u0),a.up*1/a.u0)
sin(a::Dual)=Dual(sin(a.u0),a.up*cos(a.u0))
cos(a::Dual)=Dual(cos(a.u0),-a.up*sin(a.u0))
tan(a::Dual)=Dual(tan(a.u0),(-sec(a.u0)^2)*a.up)
sqrt(a::Dual)=Dual(sqrt(a.u0),-a.up/sqrt(a.u0))

end
