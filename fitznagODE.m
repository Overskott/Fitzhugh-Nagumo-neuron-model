
function [T,V] = fitznagODE(R, I_app, a, b, tau, IC, t)

% Function that solves the Fizthugh-Nagumo ODE with parameters.

[T,V] = ode45(@fn, t, IC);

    function func = fn(~, IC)
        
    v = IC(1);
    w = IC(2);
    
    func = [v-((v.^3)./3)-w+R.*I_app; (v+a-b.*w).*(1/tau)];

    end
end