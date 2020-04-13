%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function param = set_params6_samek(X,param)
    
    %%%! Assimilation efficiency lambda (constant across everything)
    param.Lambda = X(1);

    %%%! Metabolism constants (activity and basal)
    param.amet = X(4);    % coeff on met
    param.gam  = X(5);    % coeff on search area
    param.kc   = X(6);    % coeff on cmax T-dep fn 
    param.ke   = X(6);    % coeff on enc T-dep fn
    param.kt   = X(6);    % coeff on met T-dep fn (orig 0.063)
    param.bpow = X(2);    % power on metab fn (orig 0.25)
    param.benc = X(3);    % power on enc fn (orig 0.20)
    
    
%-----
end
