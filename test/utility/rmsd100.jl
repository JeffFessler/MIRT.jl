# test/utility/rmsd100

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
