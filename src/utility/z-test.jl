# utility/z-test.jl

using Test

@test downsample(:test)
@test map_many(:test)

@test mask_or(:test)
@test mask_outline(:test)
@test embed(:test)
@test maskit(:test)
@test ndgrid(:test)
