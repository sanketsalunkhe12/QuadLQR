import sympy
from sympy import sin, cos, Matrix, hessian
from sympy import symbols

x, y, vx, vy, dt, u1, u2, m, theta, g, omega, i, r, x_goal, y_goal, theta_goal = symbols('x y vx vy dt u1 u2 m theta g omega i r x_goal y_goal theta_goal')

z = Matrix([x,
            vx,
            y,
            vy,
            theta,
            omega])

u = Matrix([u1,
            u2])

z_des = Matrix([x_goal,
            0,
            y_goal,
            0,
            theta_goal,
            0])

Q = Matrix([[1000, 0, 0, 0, 0, 0],
           [0, 1000, 0, 0, 0, 0],
           [0, 0, 1000, 0, 0, 0],
           [0, 0, 0, 1000, 0, 0],
           [0, 0, 0, 0, 1000, 0],
           [0, 0, 0, 0, 0, 1000]])

R  = Matrix([[0.1, 0],
            [0, 0.1]])

u_des = Matrix([m*g/2,
                 m*g/2])

c = (z - z_des).T * Q * (z - z_des) + (u - u_des).T * R * (u - u_des)

q_ = c.jacobian(z) #q
r_ = c.jacobian(u) #r
Q_ = hessian(c, z) #Q
R_ = hessian(c, u) #R
S_ = hessian(c, (u,z))

sympy.pprint(r_)