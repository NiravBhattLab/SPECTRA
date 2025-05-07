problem = MILPproblem;
x=zeros(17,1);
x([1,2,9])=1;
x([10,17])=1;
for i=1:numel(problem.b)
    if problem.csense(i)=='L'
        if problem.A(i,:)*x>problem.b(i,:)
            error(num2str(i))
        end
    elseif problem.csense(i)=='G'
        if problem.A(i,:)*x<problem.b(i,:)
            error(num2str(i))
        end
    elseif problem.csense(i)=='E'
        if problem.A(i,:)*x~=problem.b(i,:)
            error(num2str(i))
        end
    end
end
if sum(problem.lb>x)~=0
    ids = find(problem.lb>x)
    error('lower bounds are infeasible')
elseif sum(problem.ub<x)~=0
    ids = find(problem.ub<x)
    error('upper bounds are infeasible')
end