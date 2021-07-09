# matlab variant of function to compute LU factorisation of Lower Hessenberg matrix (WLOG) in 
# optimal time complexity with principles of numerical and conventional linear algebra

function inputmat=LUHessenberg(inputmat)
    n =length(inputmat);
  for i=1:(n-1)
    for j = (i+1):n
      inputmat(j,i) <- inputmat(j,i)/inputmat(i,i);
      inputmat(j,i+1) <- inputmat(j,i+1) - inputmat(j,i)*inputmat(i,i+1);
    end
  end
  inputmat;
end
#output is both L and U matrices composed into one dense matrix, per space complexity of numerical linear algebra

v = [1 2 3 4 5]
X = tril(v + v’) + diag(-v(1:end-1),1) #where v is an arbitrary vector choice

#test LUHessenberg(X)
