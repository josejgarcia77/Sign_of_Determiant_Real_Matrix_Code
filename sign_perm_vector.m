function sgn = sign_perm_vector(P)
% SIGN_PERM_VECTOR Calculates the sign of a permutation vector p.
% P is a row vector P(1,n), which represents the permutation.
% sign(P) = (-1)^(No. of even-length cycles)
% Complexity : O(n + ncyc) ~ O(n + Hn) ~~ O(n+log(n)) steps.
%
% Minor modification of the code by Derek O'Connor 20 March 2011.
%
n   = length(P);
visited(1:n) = false;                  % Logical vector which marks all p(k)
                                       % not visited.
sgn = 1;
for k = 1:n
    if ~visited(k)                     % k not visited, start of new cycle
        next = k;
        L = 0;
        while ~visited(next)           % Traverse the current cycle k
            L = L+1;                   % and find its length L
            visited(next) =  true;
            next    = P(next);
        end
        if rem(L,2) == 0               % If L is even, change sign.
            sgn = -sgn;
        end
    end % if ~visited(k)
end % for k 
