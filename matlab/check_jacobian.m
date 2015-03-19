function check_jacobian(Jac,M,lambda)

check_mat = M + lambda*Jac;

c = condest(check_mat);

fprintf('Condition number of the jacobian is about %g\n',c);
end