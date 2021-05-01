function solve_ode(u0, tspan, p, time_interval=1)
        function SIR(du,u,p,t)
                du[1]=-p[1]*u[1]*u[2]/(u[1]+u[2]+u[3])
                du[2]=p[1]*u[1]*u[2]/(u[1]+u[2]+u[3]) - p[2]*u[2]
                du[3] = p[2]*u[2]
        end
        problem = ODEProblem(SIR, u0, tspan,p)
        sol=solve(problem, saveat=time_interval)
        return sol
end

function solve_ode_seasonal_forcing(u0, tspan, p,time_interval=1)
        function SIR_forcing(du,u,p,t)
                L, beta0, beta1, theta, gamma = p
                # 17 weeks (~ 4 months) seasonal?
                beta_t = beta0*(1 + beta1*cos( (2*pi*t - theta)/26 ) ) 
                # du[1]=(Ntot-u[1]-u[2])/L - beta_t*u[2]*u[1]/Ntot
                # du[2]=beta_t*u[1]*u[2]/Ntot - u[2]/D
                du[1] =  -beta_t*u[2]*u[1]/(u[1]+u[2]+u[3]) + u[3]/L
                du[2] = beta_t*u[1]*u[2]/(u[1]+u[2]+u[3]) - gamma*u[2]
                du[3] = gamma*u[2]
               
        end
        problem = ODEProblem(SIR_forcing, u0, tspan,p)
        sol=solve(problem, saveat=time_interval)
        return sol
end

function solve_seir(u0, tspan, p, time_interval=1)
        function SEIR(du,u,p,t)
                du[1]=p[1]*(u[1]+u[2]+u[3]) - p[1]*u[1] - p[2]*u[1]*u[3]/(u[1]+u[2]+u[3]+u[4])
                du[2]=p[2]*u[1]*u[3]/(u[1]+u[2]+u[3]+u[4]) - (p[1] + p[3])*u[2]
                du[3] = p[3]*u[2] - (p[4]+p[1])*u[3]
                du[4] = p[4]*u[3] - p[1]*u[4]
        end
        problem = ODEProblem(SEIR, u0, tspan,p)
        sol=solve(problem, saveat=time_interval)
        return sol
end