# test/test_all_mirt.jl

export test_all_mirt

using Test: @testset, detect_ambiguities

srclist = (
"../data",
"algorithm",
"fbp",
"io",
"mri",
"nufft",
"plot",
"regularize",
"system",
"utility",
)


"""
`test_all_mirt()`

Run all MIRT tests.
To avoid plotting during tests, first do `default(show=false)`
"""
function test_all_mirt()

	basedir = dirname(pathof(MIRT))
#	@testset "all MIRT" begin
		for root in srclist
		#	@show root
			@testset "MIRT:$root" begin
				tmp = joinpath(basedir, root, "z-test.jl")
			#	printstyled(tmp * "\n", color=:blue, bold=true)
				include(tmp)
			end
		end
#	end

	@test length(detect_ambiguities(MIRT)) == 0
	true
end
