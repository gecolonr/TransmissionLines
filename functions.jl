function x(k,N)
    return cos(k*pi/N)
end

function build_c(N)
    
    c = ones(N+1)
    c[1] = 2
    c[N+1] = 2

    return c
end

function build_D(N)
    c = build_c(N)
    
    D = zeros(N+1, N+1)
    D[1,1] = (2*N^2+1)/6
    D[N+1, N+1] = -(2*N^2+1)/6
    
    for i = 1:N-1
        D[i+1,i+1] = -1/2*(x(i,N)/(1-x(i,N)^2))
    end
    
    for i = 0:N, j = 0:N
        if i != j
            D[i+1,j+1] = c[i+1]*(-1)^(i+j)/(c[j+1]*(x(i,N) - x(j,N)))
        end
    end

    return D
end

function build_M(N, tau, A0, A1)
    
    D = build_D(N)

    C = -2*D/tau
    
    C_hat = C[1:N, 1:N+1]
    
    n, null = size(A0)
    eye = Matrix(1.0I, n, n)

    M_temp = kronecker(C_hat, eye)
    
    sec_row = zeros(n, n*(N+1))
    sec_row[1:n, 1:n] = A1
    sec_row[1:n, N*n + 1:(N+1)*n] = A0

    M = [M_temp
             sec_row]
    return M, n
end

function benchmarking()
    N = 1

    D1 = [0.5 -0.5
            0.5 -0.5]

    D1_our = build_D(N)

    @assert D1 ≈ D1_our
    
    N = 2

    D2 = [1.5 -2 0.5
            0.5 0 -0.5
                -0.5 2 -1.5]

    D2_our = build_D(N)

    @assert D2 ≈ D2_our
end
;