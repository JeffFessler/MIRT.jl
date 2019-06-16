# fatrix-block.jl
# 2019-06-15, Jeff Fessler, University of Michigan

using LinearMaps
using LinearAlgebra: I
using Test: @test

#const FatrixVector{T} = Vector{Union{LinearMap,AbstractMatrix{T}}}
if !@isdefined(FatrixVector)
	FatrixVector{T} = Vector{<:Any} # include I etc.
end

"""
`ob = block_fatrix(blocks, ...)`

Construct block_fatrix object composed of fatrix blocks.
such as block_diag(A_1, A_2, ..., A_M)

See `block_fatrix(:test)` for example usage.

in
* `blocks::FatrixVector`	array of the blocks

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

out
* `ob`	LinearMap
"""
function block_fatrix(
		blocks::FatrixVector;
		T::DataType = promote_type(eltype.(blocks)...),
		how::Symbol = :diag,
		tomo::Bool = false,
		Mkron::Int = 0,
	)

	if how == :col
		 return block_fatrix_col(blocks, T, tomo)
	elseif how == :diag
		return block_fatrix_diag(blocks, T)
	elseif how == :kron
		return block_fatrix_kron(blocks, Mkron)
	elseif how == :row
#		return block_fatrix_row(blocks)
		return hcat_lm(blocks...)
	elseif how == :sum
		return block_fatrix_sum(blocks, T)
	else
		throw(ArgumentError("unknown block type $how"))
	end
end


"""
`block_fatrix_col()`
"""
function block_fatrix_col(blocks::FatrixVector, T::DataType, tomo::Bool)

	MM = length(blocks)
	dims = zeros(Int, MM, 2)
	for mm=1:MM
		B = blocks[mm]
		if B == I
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
#			odim = blocks[mm].arg.odim
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

	# build Fatrix object
	return LinearMap{T}(
			x -> vcat([blocks[mm] * x for mm=1:MM]...),
			y -> sum([blocks[mm]' * (@view y[istart[mm]:iend[mm]]) for mm=1:MM]),
			dim[1], dim[2])
end



"""
#
# block_fatrix_col_gram()
#
function [T, reuse] = block_fatrix_col_gram(ob, W, reuse, varargin)

blocks = ob.arg.blocks
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
			if isvar('W.arg.blocks[mm]')
				T[mm] = build_gram(A, W.arg.blocks[mm], ...
					reuse, varargin{:})
			else
				fail('block_fatrix_col_gram needs block diag W')
			end
		end
	end
end
T = fatrix_plus(T{:})
"""



"""
#
# block_fatrix_row()
# trick: just reuse col via transpose!
#
function block_fatrix_row(blocks, arg)

for mm=1:length(blocks)
	blocks[mm] = blocks[mm]'; % trick: transpose
end

ob = block_fatrix(blocks, 'type', 'col')'; % trick: transpose
"""


"""
`block_fatrix_diag()`
"""
function block_fatrix_diag(blocks::FatrixVector, T::DataType)

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

	# build Fatrix object
	return LinearMap{T}(
		x -> vcat([blocks[mm] * x[jstart[mm]:jend[mm]] for mm=1:MM]...),
		y -> vcat([blocks[mm]' * y[istart[mm]:iend[mm]] for mm=1:MM]...),
		dim[1], dim[2])
end


"""
#
# block_fatrix_diag_mtimes_block()
# caution: it is quite unclear whether these are useful!
#
function y = block_fatrix_diag_mtimes_block(arg, is_transpose, x, istart, nblock)
if is_transpose
	y = block_fatrix_diag_block_back(arg, x, istart, nblock)
else
	y = block_fatrix_diag_block_forw(arg, x, istart, nblock)
end


#
# block_fatrix_diag_block_forw()
#
function y = block_fatrix_diag_block_forw(arg, x, istart, nblock)

if nrow(x) ~= arg.dim(2)
	error('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
end
y = []
for mm=istart:nblock:length(arg.blocks)
	t = arg.blocks[mm] * x([arg.jstart(mm):arg.jend(mm)], :)
	y = [y; t]
end


#
# block_fatrix_diag_block_back()
#
function x = block_fatrix_diag_block_back(arg, y, istart, nblock)

if nrow(y) ~= arg.dim(1), error 'bad y size', end
x = []
for mm=istart:nblock:length(arg.blocks)
	t = arg.blocks[mm]' * y([arg.istart(mm):arg.iend(mm)], :)
	x = [x; t]
end


#
# block_fatrix_diag_gram()
#
function [T, reuse] = block_fatrix_diag_gram(ob, W, reuse, varargin)

blocks = ob.arg.blocks
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
T = block_fatrix(T, 'type', 'diag')
"""


#
# block_fatrix_kron()
#
function block_fatrix_kron(blocks::FatrixVector, T::DataType, Mkron::Int)

	length(blocks) != 1 && throw(ArgumentError("kron expects exactly one block"))

	dims = repeat(size(blocks[1]), (Mkron,1))
	# start/end indices for selecting parts of x and y
	istart = cumsum([1; dims[1:end-1,1]])
	iend = istart + dims(:,1) .- 1
	jstart = cumsum([1; dims[1:end-1,2]])
	jend = jstart + dims[:,2] .- 1
	dim = sum(dims,1)

	return LinearMap{T}(
			x -> vcat([blocks[1] * x[jstart[mm]:jend[mm]] for mm=1:Mkron]...),
			y -> vcat([blocks[1]' * y[istart[mm]:iend[mm]] for mm=1:Mkron]...),
			dim[1], dim[2])
end


"""
#
# block_fatrix_kron_mtimes_block(): y = A{ib} * x
#
function y = block_fatrix_kron_mtimes_block(arg, is_transpose, x, iblock, nblock)

bl1 = arg.blocks{1}; % base block, already put through Gblock
warn 'todo: size check not done'

if ~is_transpose % forw
#	if nrow(x) ~= arg.dim(2)
#		fail('x size=%d vs dim(2)=%d', nrow(x), arg.dim(2))
#	end
	y = []
	for mm=1:arg.Mkron
		t = bl1{iblock} * x([arg.jstart(mm):arg.jend(mm)], :)
		y = [y; t]
	end

else % back
	y = x
#	if nrow(y) ~= arg.dim(1), error 'bad y size', end
	x = []
	for mm=1:arg.Mkron
		tmp = [arg.istart(mm):arg.iend(mm)]; % i list (if all data)
		% todo: we need a certain subset of that list
		fail('todo: kron subset backprojector not done')
		t = bl1{iblock}' * y(tmp, :)
		x = [x; t]
	end
	y = x
end
"""


#
# block_fatrix_sum()
#
function block_fatrix_sum(blocks::FatrixVector, T::DataType)
	dim = size(blocks[1])
	MM = length(blocks)
	for mm=1:MM
		!isequal(size(blocks[mm]), dim) &&
			throw(DimensionMismatch("blocks must have same size for :sum"))
	end

	# build Fatrix object
	return LinearMap{T}(arg,
			x -> sum([blocks[mm] * x for mm=1:MM]...),
			y -> sum([blocks[mm]' * y for mm=1:MM]...),
			dim[1], dim[2])
end


"""
#
# block_fatrix_free()
#
function block_fatrix_free(arg)
if arg.chat
	printm 'freeing block_fatrix object static memory'
end
for mm=1:length(arg.blocks)
	try
		free(blocks[mm])
	catch
	end
end
"""


"""
#
# block_fatrix_update()
#
function out = block_fatrix_update(ob, varargin)
# todo: figure out how to update, e.g., new_zmap(s), ...
"""


"""
`block_fatrix(:test)`
self test
"""
function block_fatrix(test::Symbol)
	test != :test && throw(ArgumentError("test $test"))

	# test :col
	A = ones(4,3)
	B = LinearMap(x -> A*x, y -> A'*y, 4, 3)
#	C = I # todo later
	C = [2A; 3A]
	blocks = [A, B, C]
	Tc = block_fatrix(blocks, how=:col)
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
	Td = block_fatrix(blocks, how=:diag)
	Td * ones(8)
#	Td' * ones(15)
#	@test Matrix(Td)' == Matrix(Td')

#	todo

	true
end

block_fatrix(:test) # todo
