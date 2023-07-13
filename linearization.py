import sympy
from sympy import sin, cos, Matrix
from sympy import symbols

x, y, vx, vy, dt, u1, u2, m, theta, g, omega, i, r = symbols('x y vx vy dt u1 u2 m theta g omega i r')

f = Matrix([x+dt*vx, 
            vx+(dt*((-(u1+u2)*sin(theta))/m)), 
            y+dt*vy,
            vy+(dt*((((u1+u2)*cos(theta))/m)-g)),
            theta+dt*omega,
            omega+(dt*(r*(u1-u2))/i)])

z = Matrix([x,
            vx,
            y,
            vy,
            theta,
            omega])

u = Matrix([u1,
            u2])

# sympy.pprint(f)
# sympy.pprint(z)
# sympy.pprint(u)

df_dz = f.jacobian(z)
df_du = f.jacobian(u)

sympy.pprint(df_dz)
sympy.pprint(df_du)