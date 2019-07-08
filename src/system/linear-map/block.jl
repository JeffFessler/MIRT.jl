# lm-block.jl
# block operations for linear maps
# 2019-06-15, Jeff Fessler, University of Michigan

export block_lm

using LinearMaps
using LinearAlgebra: UniformScaling, I
using Test: @test

#const BlockVector{T} = Vector{Union{LinearMap,AbstractMatrix{T}}}
if !@isdefined(BlockVector)
	BlockVector{T} = Vector{<:Any} # include I etc.
end

"""
`ob = block_lm(blocks, ...)`

Construct `LinearMap` object composed of blocks.
such as `block_diag(A_1, A_2, ..., A_M)`.
Each block can be anything like a `LinearMap` such as an `AbstractMatrix`.

See `block_lm(:test)` for example usage.

in
* `blocks::BlockVector`	array of the blocks

* option
* `how::Symbol` options:
  + `:col` for `[A_1; A_2; ...; A_M]`
  + `:diag` for block diagonal (default)
  + `:kron` for `kron(eye(Mkron), blocks[1])`
  + `:row` for `[A1, A2, ..., A_M]`
  + `:sum` for `A1 + A2 + ... + A_M`
* `T::DataType` default `promote_type(eltype.(blocks)...)`
* `Mkron::Int' required for `:kron` type
* `tomo::Bool` special support for tomography-type objects; default `false`
	(not implemented)

out
* `ob` `LinearMap`
"""
function block_lm(
		blocks::BlockVector;
		T::DataType = promote_type(eltype.(blocks)...),
		how::Symbol = :diag,
		tomo::Bool = false,
		Mkron::Int = 0,
	)

	if how == :col
		 return block_lm_col(blocks, T, tomo)
	elseif how == :diag
		return block_lm_diag(blocks, T)
	elseif how == :kron
		return block_lm_kron(blocks, T, Mkron)
	elseif how == :row
		return block_lm_row(blocks, T)
	elseif how == :sum
		return block_lm_sum(blocks, T)
	else
		throw(ArgumentError("unknown block type $how"))
	end
end


"""
`block_lm_col()`
"""
function block_lm_col(blocks::BlockVector, T::DataType, tomo::Bool)

	MM = length(blocks)
	dims = zeros(Int, MM, 2)
	for mm=1:MM
		B = blocks[mm]
		if isa(B, UniformScaling)
			throw("I not yet supported")
		end
		dims[mm,:] .= size(B)
		dims[mm,2] != dims[1,2] && throw("all blocks must have same #cols for :col")
	end

	# start/end indices for selecting parts of x and y
	istart = cumsum([1; dims[1:end-1,1]])
	iend = istart + dims[:,1] .- 1

	dim = [sum(dims[:,1]), dims[1,2]]

#	if tomo # prep for mat2cell later
#		for mm=1:MM
#			odim = blocks[mm].odim
#			tomo_ndim = length(odim)
#			tomo_dim = tomo_ndim # insist on last dimension
#			if mm==1
#				for id = 1:tomo_ndim
#					if id != tomo_dim
#						mat2cell_arg[id] = odim[id]
#					end
#				end
#			end
#			mat2cell_arg[tomo_dim][mm] = odim(tomo_dim)
#		end
#	end

	return LinearMap{T}(
			x -> vcat([blocks[mm] * x for mm=1:MM]...),
			y -> sum([blocks[mm]' * (@view y[istart[mm]:iend[mm]]) for mm=1:MM]),
			dim[1], dim[2])
end



"""
#
# block_lm_col_gram()
#
function [T, reuse] = block_lm_col_gram(ob, W, reuse, varargin)

blocks = ob.blocks
T = cell(size(blocks))
for mm=1:length(blocks)
	A = blocks[mm]
	if isnumeric(A)
		if isempty(W)
			T[mm] = A' * A
		else
			warn 'todo: this may not work, need piece of W!'
			T[mm] = A' * W * A
		end
	else
		if isempty(W)
			T[mm] = build_gram(A, [], reuse, varargin{:})
		else
			if isvar('W.blocks[mm]')
				T[mm] = build_gram(A, W.blocks[mm], ...
					reuse, varargin{:})
			else
				fail('block_lm_col_gram needs block diag W')
			end
		end
	end
end
T = block_lm(T, how=:sum)
"""


"""
`block_lm_row()`
trick: just reuse col via transpose!
cannot use hcat(blocks...) because there might be a mix of Matrix/LinearMaps
"""
function block_lm_row(blocks::BlockVector, T::DataType)
	tblocks = [b' for b in blocks] # trick: transpose
	return block_lm(tblocks, how=:col)' # trick: transpose
end


"""
`block_lm_diag()`
"""
function block_lm_diag(blocks::BlockVector, T::DataType)

	MM = length(blocks)
#	@show dims = vcat(map(size, blocks))
	dims = zeros(Int, MM, 2)
	for mm=1:MM
		dims[mm,:] .= size(blocks[mm])
	end

	# start/end indices for selecting parts of x and y
	istart = cumsum([1; dims[1:end-1,1]])
	iend = istart + dims[:,1] .- 1
	jstart = cumsum([1; dims[1:end-1,2]])
	jend = jstart + dims[:,2] .- 1

	dim = sum(dims, dims=1)

	return LinearMap{T}(
		x -> vcat([blocks[mm] * x[jstart[mm]:jend[mm]] for mm=1:MM]...),
		y -> vcat([blocks[mm]' * y[istart[mm]:iend[mm]] for mm=1:MM]...),
		dim[1], dim[2])
end


"""
#
# block_lm_diag_mtimes_block()
# caution: it is quite unclear whether these are useful!
#
function y = block_lm_diag_mtimes_block(arg, is_transpose, x, istart, nblock)
if is_transpose
	y = block_lm_diag_block_back(arg, x, istart, nblock)
else
	y = block_lm_diag_block_forw(arg, x, istart, nblock)
end


#
# block_lm_diag_block_forw()
#
function y = block_lm_diag_block_forw(arg, x, istart, nblock)

if nrow(x) ~= dim(2)
	throw('x size=%d vs dim(2)=%d', nrow(x), dim(2))
end
y = []
for mm=istart:nblock:length(blocks)
	t = blocks[mm] * x([jstart(mm):jend(mm)], :)
	y = [y; t]
end


#
# block_lm_diag_block_back()
#
function x = block_lm_diag_block_back(arg, y, istart, nblock)

if nrow(y) ~= dim(1), error 'bad y size', end
x = []
for mm=istart:nblock:length(blocks)
	t = blocks[mm]' * y([istart(mm):iend(mm)], :)
	x = [x; t]
end


#
# block_lm_diag_gram()
#
function [T, reuse] = block_lm_diag_gram(ob, W, reuse, varargin)

blocks = ob.blocks
T = cell(size(blocks))
for mm=1:length(blocks)
	A = blocks[mm]
	if isnumeric(A)
		if isempty(W)
			T[mm] = A' * A
		else
			T[mm] = A' * W * A
		end
	else
		T[mm] = build_gram(A, W, reuse, varargin{:})
	end
end
T = block_lm(T, 'type', 'diag')
"""


#
# block_lm_kron()
#
function block_lm_kron(blocks::BlockVector, T::DataType, Mkron::Int)

	length(blocks) != 1 && throw(ArgumentError("kron expects one block"))

	dims = repeat(collect(size(blocks[1]))', Mkron, 1)
	# start/end indices for selecting parts of x and y
	istart = cumsum([1; dims[1:end-1,1]])
	iend = istart + dims[:,1] .- 1
	jstart = cumsum([1; dims[1:end-1,2]])
	jend = jstart + dims[:,2] .- 1
	dim = sum(dims,dims=1)

	return LinearMap{T}(
			x -> vcat([blocks[1] * x[jstart[mm]:jend[mm]] for mm=1:Mkron]...),
			y -> vcat([blocks[1]' * y[istart[mm]:iend[mm]] for mm=1:Mkron]...),
			dim[1], dim[2])
end


"""
#
# block_lm_kron_mtimes_block(): y = A{ib} * x
#
function y = block_lm_kron_mtimes_block(arg, is_transpose, x, iblock, nblock)

bl1 = blocks{1}; % base block, already put through Gblock
warn 'todo: size check not done'

if ~is_transpose % forw
#	if nrow(x) ~= dim(2)
#		fail('x size=%d vs dim(2)=%d', nrow(x), dim(2))
#	end
	y = []
	for mm=1:Mkron
		t = bl1{iblock} * x([jstart(mm):jend(mm)], :)
		y = [y; t]
	end

else % back
	y = x
#	if nrow(y) ~= dim(1), error 'bad y size', end
	x = []
	for mm=1:Mkron
		tmp = [istart(mm):iend(mm)]; % i list (if all data)
		% todo: we need a certain subset of that list
		fail('todo: kron subset backprojector not done')
		t = bl1{iblock}' * y(tmp, :)
		x = [x; t]
	end
	y = x
end
"""


#
# block_lm_sum()
#
function block_lm_sum(blocks::BlockVector, T::DataType)
	dim = size(blocks[1])
	MM = length(blocks)
	for mm=1:MM
		!isequal(size(blocks[mm]), dim) &&
			throw(DimensionMismatch("blocks must have same size for :sum"))
	end

	return LinearMap{T}(
			x -> sum([blocks[mm] * x for mm=1:MM]),
			y -> sum([blocks[mm]' * y for mm=1:MM]),
			dim[1], dim[2])
end


"""
#
# block_lm_free()
#
function block_lm_free(arg)
if chat
	printm 'freeing block_lm object static memory'
end
for mm=1:length(blocks)
	try
		free(blocks[mm])
	catch
	end
end
"""


"""
#
# block_lm_update()
#
function out = block_lm_update(ob, varargin)
# todo: figure out how to update, e.g., new_zmap(s), ...
"""


"""
`block_lm(:test)`
self test
"""
function block_lm(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))

	# test :col
	A = ones(4,3)
	B = LinearMap(x -> A*x, y -> A'*y, 4, 3)
#	C = I # todo later
	C = [2A; 3A]
	blocks = [A, B, C]
	Tc = block_lm(blocks, how=:col)
	Tc * ones(3)
	Tc' * ones(16)
	@test Matrix(Tc)' == Matrix(Tc')

	# test :diag
	A = ones(3,2)
	M = ones(4,3)
	B = LinearMap(x -> M*x, y -> M'*y, 4, 3)
#	C = I # todo later
	C = [2M; 3M]
	blocks = [A, B, C]
	Td = block_lm(blocks, how=:diag)
	Td * ones(8)
	Td' * ones(15)
	@test Matrix(Td)' == Matrix(Td')

	# test :kron
	Tk = block_lm([ones(3,2)], how=:kron, Mkron=4)
	Tk * ones(4*2)
	Tk' * ones(4*3)
	@test Matrix(Tk)' == Matrix(Tk')

	# test :row
	A = ones(4,3)
	B = LinearMap(x -> A*x, y -> A'*y, 4, 3)
#	C = I # todo later
	C = [2A 3A]
	blocks = [A, B, C]
	Tr = block_lm(blocks, how=:row)
	Tr * ones(12)
	Tr' * ones(4)
	@test Matrix(Tr)' == Matrix(Tr')

	# test :sum
	A = ones(4,3)
	B = LinearMap(x -> A*x, y -> A'*y, 4, 3)
#	C = I # todo later
	C = 2A
	blocks = [A, B, C]
	Ts = block_lm(blocks, how=:sum)
	Ts * ones(3)
	Ts' * ones(4)
	@test Matrix(Ts)' == Matrix(Ts')

	@test_throws ArgumentError block_lm(blocks, how=:bad)
	@test_throws String block_lm([I], how=:col)

	true
end

# block_lm(:test)
