function solve_lorenz(u0, tspan, p, time_interval=0.01)
    function lorenz(du, u, p, t)
            du[1] = p[1]*(u[2]-u[1])
            du[2] = u[1]*(p[2] - u[3]) - u[2]
            du[3] = u[1]*u[2] - p[3]*u[3]
    end
    problem = ODEProblem(lorenz, u0, tspan, p)
    sol = solve(problem, saveat=time_interval)
    return sol
end