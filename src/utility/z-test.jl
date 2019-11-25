# utility/z-test.jl

using Test

@test downsample(:test)
@test map_many(:test)
@test mask_test()
@test ndgrid(:test)
@test interp1(:test)
@test jinc(:test)
@test max_percent_diff(:test)
@test rect(:test)
