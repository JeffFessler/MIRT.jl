# io/z-test.jl

using Test: @testset

list = [
"caller_name.jl"
"fld-read.jl"
"ir_dump.jl"
"prompt.jl"
"shows.jl"
"fld-write.jl"
]

# todo: temporary work around for:
# https://github.com/JeffFessler/MIRT.jl/issues/58
if Base.Sys.iswindows()
	list = list[1:(end-1)]
end

for file in list
	@testset "$file" begin
		include(file)
	end
end
