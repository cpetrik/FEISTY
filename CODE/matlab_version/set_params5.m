%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function set_params5(X)
    global Lambda gam bpow benc amet 
    
    %%%! Assimilation efficiency lambda (constant across everything)
    Lambda = X(1);

    %%%! Metabolism constants (activity and basal)
    amet = X(4);    % coeff on met
    gam  = X(5);    % coeff on search area
    bpow = X(2);    % power on metab fn (orig 0.25)
    benc = X(3);    % power on enc fn (orig 0.20)
    
    
%-----
end
