# rmsd100

using MIRT: rmsd100

using Test: @test, @test_throws, @inferred


a = [1 2; 4 5; 7 8]
b = [3 2; 4 6; 8 7]
c = [false true; true false; false false] #avg. 0
d = [false false; false false; false true] #avg. 1
e = [true true; true true; true false] # average 1,0,1,0,4 = sqrt(6/5)
tests = cat(c, d, e, dims = 3)
pred = [0, 100, 100*sqrt(6/5)]
for i in 1:size(tests,1)
	@test rmsd100(a, b ; mask=tests[:,:,i]) == pred[i]
end
@test rmsd100([1, 0], [0, 1]) == 100 # no mask test

#=  not important to support or test.  user can pass mask=ig.mask if needed
a = ones(300,400)
b = ones(300,400)
b[100:200,100:200] .= 0
#so there's a blot of zeros from 100,100 to 200,200
b[200:250,300:320] .= 3
#and a blot of threes over here
ig1 = image_geom(nx = 300,ny = 400,dx = 1,dy = 1)
@test my_rmse(a,b,ig1) > 20
=#        
