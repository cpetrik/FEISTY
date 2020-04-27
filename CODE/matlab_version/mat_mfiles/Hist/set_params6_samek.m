%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function set_params6_samek(X)
    global Lambda gam bpow benc amet kt kc ke
    
    %%%! Assimilation efficiency lambda (constant across everything)
    Lambda = X(1);

    %%%! Metabolism constants (activity and basal)
    amet = X(4);    % coeff on met
    gam  = X(5);    % coeff on search area
    kc   = X(6);    % coeff on cmax T-dep fn 
    ke   = X(6);    % coeff on enc T-dep fn
    kt   = X(6);    % coeff on met T-dep fn (orig 0.063)
    bpow = X(2);    % power on metab fn (orig 0.25)
    benc = X(3);    % power on enc fn (orig 0.20)
    
    
%-----
end
