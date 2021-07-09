# matlab variant of function to compute LU factorisation of Lower Hessenberg matrix (W.L.O.G.) in 
2
# optimal time complexity with principles of numerical and conventional linear algebra
3
​
4
function inputmat=LUHessenberg(inputmat)
5
    n =length(inputmat);
6
  for i=1:(n-1)
7
    for j = (i+1):n
8
      inputmat(j,i) <- inputmat(j,i)/inputmat(i,i);
9
      inputmat(j,i+1) <- inputmat(j,i+1) - inputmat(j,i)*inputmat(i,i+1);
10
    end
11
  end
12
  inputmat;
13
end
14
#output is both L and U matrices composed into one dense matrix, per space complexity of numerical linear algebra
15
​
16
v = [1 2 3 4 5]
17
X = tril(v + v’) + diag(-v(1:end-1),1) #where v is an arbitrary vector chocie
18
#test LUHessenberg(X)
19
​
