using SparseArrays
using LinearAlgebra
using Combinatorics

thop = 1 # initial hopping between site 1 and 2, J=t0+delta0
delta = 1.
Delta = 2 # initial local potential 
T = 0.01 #parse(Float64,ARGS[1]) # temperature

v = 0.1 #parse(Float64, ARGS[2]) # pumping speed
N = 1e2 # number of intervals for time evolution
U = 7 #parse(Float64, ARGS[3])
L= 4 #parse(Int64, ARGS[4]) # number of unit cells



n1 = Diagonal([0,1,1,2]) 
id = Diagonal([1,1,1,1])
N1=id
for i=1:L-1
    global N1=kron(N1,id)
end
N1=kron(N1,n1)
for i=1:L-1
    global N1=kron(N1,id)
end
# N1 is the number of atoms at a single site for sublattice A

size=binomial(2L,L)^2  # dimensino of Hamiltonian for total S_z=0 sector
# some parameters used to convert between binary and decimal numbers
expn=collect(0:1:4*L-1)
base=2*ones(Int64,4*L)
conv=base.^expn
two=2*ones(Int64,2L)
index=zeros(Int64,size)

# We choose a basis by ordering the creation operators as {A,up}, {A,down}, {B,up}, {B,down} for ascending unit cell index I. 
# This is the equivalent to a Jordan-Wigner transformation.

M=10
P=zeros(M)
ch=zeros(M)
c=zeros(M)

for q=1:M
    global U=0.5q # we change the interactions 
    dt=1/(N*v)
    H=zeros(size,size)
    evenL=collect(1:2L)*2       # site with even indices for spin down
    oddL=evenL-ones(Int64,2L)   # site with odd indices for spin up
    evenB=collect(combinations(evenL,L))    # spin down half-filled
    oddB=collect(combinations(oddL,L))      # spin up half-filled
    basis=[[] for i=1:size]
    counter=1
    for i=1:length(oddB)
        for j=1:length(evenB)
            basis[counter]=vcat(oddB[i],evenB[j])
            counter=counter+1
        end
    end
    #println(basis)
    
    # Diagonal elements
    for I=1:size
        global index[I]=sum(two.^(basis[I]-ones(Int64,2L))) # index in base 10
        ket=digits(index[I],pad=4L,base=2) # convert combination basis into ket in {0,1} basis
        onset=-Delta*sum(repeat([1,1,-1,-1],outer=L).*ket) # different sublattices have different local potentials
        intcount=zeros(Int64,2L)
        for j=1:2L
            intcount[j]=ket[2j-1]*ket[2j]   # n_{j,up}*n_{j,down}
        end
        int=U*sum(intcount)     # Hubbard interactions
        H[I,I]=onset+int
    end
    # hopping
    for I=1:size
        ket=digits(index[I],pad=4L,base=2)
        for j=1:4L-2
            if ket[j]==1&&ket[j+2]==0
                It=Int(sum(ket.*conv))-2^(j-1)+2^(j+1)
                K=findall(x->x==It,index)[1]
                H[I,K]=H[K,I]=-thop+(-1)^Int(ceil(j/2))*delta
            end
        end
    end
    #= edge hopping (turn on for PBC; turn off for OBC)
    if L>1
        for I=1:size
            ket=digits(index[I],pad=4L,base=2)
            if ket[end-1]==1&&ket[1]==0
                It=Int(sum(ket.*conv))-2^(4L-2)+1
                K=findall(x->x==It,index)[1]
                H[I,K]=H[K,I]=-thop+delta
            end
            if ket[end]==1&&ket[2]==0
                It=Int(sum(ket.*conv))-2^(4L-1)+2
                K=findall(x->x==It,index)[1]
                H[I,K]=H[K,I]=-thop+delta
            end
        end
    end
    =#
    #display(H)
    println("H generated")
    println(size)
    flush(stdout)
    
    rho=zeros(ComplexF64,size,size)
    E0 = eigvals(H)
    vec0 = eigvecs(H)
    for i=1:size
        rho=rho+exp(-(E0[i]-E0[1])/T)*kron(vec0[:,i],vec0[:,i]') # unnormalized density matrix
    end

    rho=rho/tr(rho)     # initial density matrix in the total S_z=0 sector basis
    rhoE=spzeros(ComplexF64,16^L,16^L)
    for I=1:size
        for J=1:size
            rhoE[index[I]+1,index[J]+1]=rho[I,J] # initial density matrix in the tenor product bais
        end
    end
    nA=real(tr(rhoE*N1))  # particle number at a single site for sublattice A
    println(nA)
    Itot=0
    println("initialized")

    flush(stdout)

    sub=4^L
    Hrho=view(rhoE, 1:sub, 1:sub)
    for z=1:sub-1
        Hrho=Hrho+view(rhoE, 1+z*sub:(z+1)*sub, 1+z*sub:(z+1)*sub)
    end
    
    # step 1
    for i=1:N
        for I=1:size
            global index[I]=sum(two.^(basis[I]-ones(Int64,2L))) # index in base 10
            ket=digits(index[I],pad=4L,base=2) # convert combination basis into ket in {0,1} bais
            H[I,I]=-Delta*sum(repeat([1,1,-1,-1],outer=L).*ket)*(1-2*v*i*dt)
            intcount=zeros(Int64,2L)
            for j=1:2L
                intcount[j]=ket[2j-1]*ket[2j]
            end
            int=U*sum(intcount)
            H[I,I]=H[I,I]+int
        end
        # hopping
        for I=1:size
            ket=digits(index[I],pad=4L,base=2)
            for j=1:4L-2
                if ket[j]==1&&ket[j+2]==0
                    It=Int(sum(ket.*conv))-2^(j-1)+2^(j+1)
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=-thop+(-1)^Int(ceil(j/2))*delta
                end
            end
        end
        #= edge hopping (turn on for PBC; turn off for OBC)
        if L>1
            for I=1:size
                ket=digits(index[I],pad=4L,base=2)
                if ket[end-1]==1&&ket[1]==0
                    It=Int(sum(ket.*conv))-2^(4L-2)+1
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=-thop+delta
                end
                if ket[end]==1&&ket[2]==0
                    It=Int(sum(ket.*conv))-2^(4L-1)+2
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=-thop+delta
                end
            end
        end
        =#
        V=exp(-1im*H*dt)
        rho=V*rho*V'        # time evolution
    end
    for I=1:size
        for J=1:size
            rhoE[index[I]+1,index[J]+1]=rho[I,J]
        end
    end
    nAt=real(tr(rhoE*N1))
    dI1=nA-nAt      # cumulated charge for step 1
    Itot=Itot+dI1
    nA=nAt
    println(dI1)
    flush(stdout)

    # step 2
    for i=1:N
        for I=1:size
            global index[I]=sum(two.^(basis[I]-ones(Int64,2L))) # index in base 10
            ket=digits(index[I],pad=4L,base=2) # convert combination basis into ket in {0,1} bais
            H[I,I]=Delta*sum(repeat([1,1,-1,-1],outer=L).*ket)
            intcount=zeros(Int64,2L)
            for j=1:2L
                intcount[j]=ket[2j-1]*ket[2j]
            end
            int=U*sum(intcount)
            H[I,I]=H[I,I]+int
        end
        # hopping
        for I=1:size
            ket=digits(index[I],pad=4L,base=2)
            for j=1:4L-2
                if ket[j]==1&&ket[j+2]==0
                    It=Int(sum(ket.*conv))-2^(j-1)+2^(j+1)
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop+(-1)^Int(ceil(j/2))*delta)*(1-i*v*dt)
                end
            end
        end
        #= edge hopping (turn on for PBC; turn off for OBC)
        if L>1
            for I=1:size
                ket=digits(index[I],pad=4L,base=2)
                if ket[end-1]==1&&ket[1]==0
                    It=Int(sum(ket.*conv))-2^(4L-2)+1
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop+delta)*(1-i*v*dt)
                end
                if ket[end]==1&&ket[2]==0
                    It=Int(sum(ket.*conv))-2^(4L-1)+2
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop+delta)*(1-i*v*dt)
                end
            end
        end
        =#
        V=exp(-1im*H*dt)
        rho=V*rho*V'
    end
    for I=1:size
        for J=1:size
            rhoE[index[I]+1,index[J]+1]=rho[I,J]
        end
    end
    nAt=real(tr(rhoE*N1))
    dI2=nA-nAt
    Itot=Itot+dI2
    nA=nAt
    println(dI2)
    flush(stdout)

    # step 3
    for i=1:N
        for I=1:size
            global index[I]=sum(two.^(basis[I]-ones(Int64,2L))) # index in base 10
            ket=digits(index[I],pad=4L,base=2) # convert combination basis into ket in {0,1} bais
            H[I,I]=Delta*sum(repeat([1,1,-1,-1],outer=L).*ket)
            intcount=zeros(Int64,2L)
            for j=1:2L
                intcount[j]=ket[2j-1]*ket[2j]
            end
            int=U*sum(intcount)
            H[I,I]=H[I,I]+int
        end
        # hopping
        for I=1:size
            ket=digits(index[I],pad=4L,base=2)
            for j=1:4L-2
                if ket[j]==1&&ket[j+2]==0
                    It=Int(sum(ket.*conv))-2^(j-1)+2^(j+1)
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop-(-1)^Int(ceil(j/2))*delta)*i*v*dt
                end
            end
        end
        #= edge hopping (turn on for PBC; turn off for OBC)
        if L>1
            for I=1:size
                ket=digits(index[I],pad=4L,base=2)
                if ket[end-1]==1&&ket[1]==0
                    It=Int(sum(ket.*conv))-2^(4L-2)+1
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop-delta)*i*v*dt
                end
                if ket[end]==1&&ket[2]==0
                    It=Int(sum(ket.*conv))-2^(4L-1)+2
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop-delta)*i*v*dt
                end
            end
        end
        =#
        V=exp(-1im*H*dt)
        rho=V*rho*V'
    end
    for I=1:size
        for J=1:size
            rhoE[index[I]+1,index[J]+1]=rho[I,J]
        end
    end
    nAt=real(tr(rhoE*N1))
    dI3=nAt-nA
    Itot=Itot+dI3
    println(dI3)
    nA=nAt
    flush(stdout)

    # step 4
    for i=1:N
        for I=1:size
            global index[I]=sum(two.^(basis[I]-ones(Int64,2L))) # index in base 10
            ket=digits(index[I],pad=4L,base=2) # convert combination basis into ket in {0,1} bais
            H[I,I]=Delta*sum(repeat([1,1,-1,-1],outer=L).*ket)*(1-2*v*i*dt)
            intcount=zeros(Int64,2L)
            for j=1:2L
                intcount[j]=ket[2j-1]*ket[2j]
            end
            int=U*sum(intcount)
            H[I,I]=H[I,I]+int
        end
        # hopping
        for I=1:size
            ket=digits(index[I],pad=4L,base=2)
            for j=1:4L-2
                if ket[j]==1&&ket[j+2]==0
                    It=Int(sum(ket.*conv))-2^(j-1)+2^(j+1)
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop-(-1)^Int(ceil(j/2))*delta)
                end
            end
        end
        #= edge hopping (turn on for PBC; turn off for OBC)
        if L>1
            for I=1:size
                ket=digits(index[I],pad=4L,base=2)
                if ket[end-1]==1&&ket[1]==0
                    It=Int(sum(ket.*conv))-2^(4L-2)+1
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop-delta)
                end
                if ket[end]==1&&ket[2]==0
                    It=Int(sum(ket.*conv))-2^(4L-1)+2
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop-delta)
                end
            end
        end
        =#
        V=exp(-1im*H*dt)
        rho=V*rho*V'
    end
    for I=1:size
        for J=1:size
            rhoE[index[I]+1,index[J]+1]=rho[I,J]
        end
    end
    nAt=real(tr(rhoE*N1))
    dI4=nAt-nA
    Itot=Itot+dI4
    println(dI4)
    nA=nAt
    flush(stdout)

    # step 5
    for i=1:N
        for I=1:size
            global index[I]=sum(two.^(basis[I]-ones(Int64,2L))) # index in base 10
            ket=digits(index[I],pad=4L,base=2) # convert combination basis into ket in {0,1} bais
            H[I,I]=-Delta*sum(repeat([1,1,-1,-1],outer=L).*ket)
            intcount=zeros(Int64,2L)
            for j=1:2L
                intcount[j]=ket[2j-1]*ket[2j]
            end
            int=U*sum(intcount)
            H[I,I]=H[I,I]+int
        end
        # hopping
        for I=1:size
            ket=digits(index[I],pad=4L,base=2)
            for j=1:4L-2
                if ket[j]==1&&ket[j+2]==0
                    It=Int(sum(ket.*conv))-2^(j-1)+2^(j+1)
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop-(-1)^Int(ceil(j/2))*delta)*(1-i*v*dt)
                end
            end
        end
        #= edge hopping (turn on for PBC; turn off for OBC)
        if L>1
            for I=1:size
                ket=digits(index[I],pad=4L,base=2)
                if ket[end-1]==1&&ket[1]==0
                    It=Int(sum(ket.*conv))-2^(4L-2)+1
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop-delta)*(1-i*v*dt)
                end
                if ket[end]==1&&ket[2]==0
                    It=Int(sum(ket.*conv))-2^(4L-1)+2
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop-delta)*(1-i*v*dt)
                end
            end
        end
        =#
        V=exp(-1im*H*dt)
        rho=V*rho*V'
    end
    for I=1:size
        for J=1:size
            rhoE[index[I]+1,index[J]+1]=rho[I,J]
        end
    end
    nAt=real(tr(rhoE*N1))
    dI5=nAt-nA
    Itot=Itot+dI5
    println(dI5)
    nA=nAt
    flush(stdout)

    # step 6
    for i=1:N
        for I=1:size
            global index[I]=sum(two.^(basis[I]-ones(Int64,2L))) # index in base 10
            ket=digits(index[I],pad=4L,base=2) # convert combination basis into ket in {0,1} bais
            H[I,I]=-Delta*sum(repeat([1,1,-1,-1],outer=L).*ket)
            intcount=zeros(Int64,2L)
            for j=1:2L
                intcount[j]=ket[2j-1]*ket[2j]
            end
            int=U*sum(intcount)
            H[I,I]=H[I,I]+int
        end
        # hopping
        for I=1:size
            ket=digits(index[I],pad=4L,base=2)
            for j=1:4L-2
                if ket[j]==1&&ket[j+2]==0
                    It=Int(sum(ket.*conv))-2^(j-1)+2^(j+1)
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop+(-1)^Int(ceil(j/2))*delta)*i*v*dt
                end
            end
        end
        #= edge hopping (turn on for PBC; turn off for OBC)
        if L>1
            for I=1:size
                ket=digits(index[I],pad=4L,base=2)
                if ket[end-1]==1&&ket[1]==0
                    It=Int(sum(ket.*conv))-2^(4L-2)+1
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop+delta)*i*v*dt
                end
                if ket[end]==1&&ket[2]==0
                    It=Int(sum(ket.*conv))-2^(4L-1)+2
                    K=findall(x->x==It,index)[1]
                    H[I,K]=H[K,I]=(-thop+delta)*i*v*dt
                end
            end
        end
        =#
        V=exp(-1im*H*dt)
        rho=V*rho*V'
    end
    for I=1:size
        for J=1:size
            rhoE[index[I]+1,index[J]+1]=rho[I,J]
        end
    end
    nAt=real(tr(rhoE*N1))
    dI6=nA-nAt
    Itot=Itot+dI6
    println(dI6)
    nA=nAt
    flush(stdout)

    ch[q]=dI3+dI4+dI5   # cumulated current for steps 3-5
    c[q]=Itot           # total cumulated charge
    println(ch[q])
    println(c[q])
    flush(stdout)
    
    Hrho=view(rhoE, 1:sub, 1:sub)
    for z=1:sub-1
        Hrho=Hrho+view(rhoE, 1+z*sub:(z+1)*sub, 1+z*sub:(z+1)*sub)  # partial trace
    end
    #display(Hrho)
    P[q]=real(tr(Hrho*Hrho))
    #
end

print("ch=")
println(ch)
print("Purity=")
println(P)
print("Total Charge=")
println(c)