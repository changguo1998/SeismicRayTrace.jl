# SeismicRayTrace.jl

This module is used to calculate travel time between two points in a layerd model.

## Usage

```julia
using SeismicRayTrace
depth1 = 0.0 # depth of first point
depth2 = 1.0 # depth of another point
distance = 10.0 # horizontal distance
model_depth = [0.0, 1.0, 2.0] # depth of each interface
model_velocity = [1.0, 2.0, 4.0] # velocity under the related interface
phase_type = ["all"] # can be refraction, reflection and guide
travel_time = raytrace(depth1, depth2, distance, model_depth, model_velocity, phase_type)
```

## Return Value

The return values are like:

```julia
(x = 4.0, t = 4.0, p = 1.0, l = 0, type = "refraction"),
(x = 3.99999,t = 4.47213,p = 0.8944267450394816,l = 2,type = "reflection"),
...
(x = 4.0, t = 3.7320508075688776, p = 0.5, l = 2, type = "guide")
...
```

The program first solve the equation to get `ray parameter p`,
and then get the `horitonzal distance x` and `travel time t`.
`l` means the layer on which the ray reflect or became guide wave.
