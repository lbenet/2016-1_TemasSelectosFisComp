
module AutDiff

export Dual, dual_num, dual_var



#Se define el vector dual.
type Dual{T<:Number} 
    f::T
    df::T
end

#Se promueve el par ordenado (a,b) para que tengan el mismo tipo.
Dual(a,b) = Dual(promote(a,b)...)

#Función que convierte de real a Dual.
dual_num(x::Real) = Dual(x,0)

#Función que devuelve el Dual de la identidad evaluado en x0.
dual_var(x0) = Dual(x0,1)




#LA siguiente línea es para no editar arbitrariamente las funciones.
import Base: +,-,*,/,^, sin, cos, exp, log
#Ahora se definen las operaciones aritméticas básicas entre los duales.

+(u::Dual,v::Dual) = Dual(u.f + v.f, u.df + v.df)
+(u::Dual,a::Real) = u + dual_num(a)
+(a::Real,u::Dual) = +(u::Dual,a::Real)


-(u::Dual,v::Dual) = Dual(u.f - v.f, u.df - v.df)
-(u::Dual,a::Real) = u-dual_num(a)
-(a::Real,u::Dual) = dual_num(a) - u


*(u::Dual,v::Dual) = Dual(u.f * v.f, u.df * v.f + v.df * u.f)
*(u::Dual,a::Real) = u*dual_num(a)
*(a::Real,u::Dual) = *(u::Dual,a::Real)

#La división es especial, se difine en más de una línea.
function /(u::Dual,v::Dual)
    w = u.f/v.f
    Dual(w,(u.df - w*v.df)/v.f)
    
end

/(u::Dual,a::Real) = u/dual_num(a)
/(a::Real,u::Dual) = dual_num(a)/u


^(u::Dual,n::Float64) = Dual(u.f^n, n*u.f^(n-1) * u.df)

sin(u::Dual) = Dual(sin(u.f), u.df*cos(u.f))

cos(u::Dual) = Dual(cos(u.f), -u.df*sin(u.f))

exp(u::Dual) = Dual(exp(u.f), u.df*exp(u.f))

log(u::Dual) = Dual(log(u.f), u.df/u.f)

end

