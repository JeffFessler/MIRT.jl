# utility/z-test.jl

using Test

@test downsample(:test)
@test map_many(:test)
@test mask_test()
@test ndgrid(:test)
