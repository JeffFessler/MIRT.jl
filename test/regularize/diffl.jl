# diffl.jl

using MIRT: diffl, diffl!, diffl_adj, diffl_adj!, diffl_map

using LinearMapsAA: LinearMapAM, LinearMapAO
using Test: @test, @testset, @test_throws, @inferred



@testset "1D" begin
	x = rand(4)
	@test diffl(x)[2:end] == diff(x)
end

@testset "2D" begin
	x = [2 4; 6 16]
	g = @inferred diffl(x, 1)
	@test all(g[1,:] .== 0)
	@test g[2:end,:] == diff(x, dims=1)
	g = @inferred diffl(x, 2)
	@test all(g[:,1] .== 0)
	@test g[:,2:end] == diff(x, dims=2)

	g = @inferred diffl(x, 1 ; edge=:none)
	@test g[2:end,:] == diff(x, dims=1)

	g = @inferred diffl(x, 1 ; add=true)
	@test g[2:end,:][:] == x[1,:] + x[2,:]

	x = rand(3,4)
	g1 = diffl(x)
	g2 = similar(g1)
	@inferred diffl!(g2, x) # test default 1
	@test g2[2:end,:] == g1[2:end,:]
end

@testset "stack" begin
	x = reshape((1:(2*3*4)).^2, 2, 3, 4)
	g = diffl(x, [3, 1])
	@test all(g[:,:,1,1] .== 0)
	@test all(g[1,:,:,2] .== 0)
	@test g[:,:,2:end,1] == diff(x, dims=3)
	@test g[2:end,:,:,2] == diff(x, dims=1)
end

@testset "adj" begin
	x = rand(3)
	@inferred diffl_adj(rand(3))
	@test_throws ArgumentError diffl_adj(rand(3) ; edge=:test)
#	@inferred diffl_adj(rand(3,4,2), 1:2) # fails
	@test size(diffl_adj(rand(3,4,2), 1:2)) == (3,4)
end




N = (2,3); d = 2

@testset "basics" begin
	O = diffl_map(N, d ; T=Int32, edge=:zero, add=false)
	@test O isa LinearMapAM
	@test Matrix(O)' == Matrix(O')
	@test O.name == "diffl_map"
	@test_throws String diffl_map(N ; edge=:none) # unsupported
end

@testset "adjoint" begin
	for N in [(3,), (3,4), (2,3,4)]
		dlist = [1, [1,]]
		length(N) > 1 && push!(dlist, 2:-1:1, 1:length(N))
		length(N) > 2 && push!(dlist, [length(N), 1])
		isinteractive() && @show dlist
		for d in dlist
			for edge in (:zero, :circ)
				for add in (false, true)
					for op in (false, true)
						O = diffl_map(N, d ;
							T=Int32, edge=edge, add=add, operator = op)
						@test Matrix(O)' == Matrix(O')
						@test O isa (op ? LinearMapAO : LinearMapAM)
					end
				end
			end
		end
	end
end
