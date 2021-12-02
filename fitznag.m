function func = fitznag(~, IC)
    
    % Fizthugh-Nagumo with hard coded parameters

    % Parameters from https://en.wikipedia.org/w/index.php?title=FitzHugh%E2%80%93Nagumo_model&oldid=1037464531
    R = 1;
    I_app = 0.5;
    a = 0.8;
    b = 0.7;
    tau = 12.5;
    
    v = IC(1);
    w = IC(2);
    
    func = [v-((v.^3)./3)-w+R.*I_app; (v+a-b.*w).*(1./tau)];
end