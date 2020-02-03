using Test: @test, @test_throws, @inferred
"""
    rmse = my_rmse(I,ref,ig)
    Generate 100 * RMSE (root mean squared error) of I compared to ref within domain ig.mask.

in
- I : Computed matrix
- ref : Reference matrix
- ig : image_geom object representing the domain of the RMSE

out
- rmse : rmse of I and ref within mask ig.
"""
function my_rmse(I::Array,ref::Array,ig::MIRT_image_geom)
    return 100 * norm(I[ig.mask] - ref[ig.mask]) / sqrt(sum(ig.mask[:]))
end
function my_rmse(I::Array,ref::Array,ig::Array{Bool})
    return 100 * norm(I[ig] - ref[ig]) / sqrt(sum(ig))
end

#TODO: ask fessler about what the hell image_geom is doing.
function my_rmse(x::Symbol)
    x != :test && throw("Invalid symbol; use :test for testing functions.")
    a = [1 2; 4 5;7 8]
    b = [3 2; 4 6;8 7]
    c = [false true;true false;false false] #avg. 0
    d = [false false;false false;false true] #avg. 1
    e = [true true;true true;true false] #average here should be 1,0,1,0,4 = sqrt(6/5) = 109...
    tests = cat(c,d,e,dims = 3)
    pred = [0 100 100*sqrt(6/5)]
    for i in 1:size(tests,1)
        @test my_rmse(a,b,tests[:,:,i]) == pred[i]
    end
    a = ones(300,400)
    b = ones(300,400)
    b[100:200,100:200] .= 0
    #so there's a blot of zeros from 100,100 to 200,200
    b[200:250,300:320] .= 3
    #and a blot of threes over here
    ig1 = image_geom(nx = 300,ny = 400,dx = 1,dy = 1)
    @test my_rmse(a,b,ig1.mask) > 20
    return true
end
