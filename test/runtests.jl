using SeismicRayTrace
using Test

@testset "SeismicRayTrace.jl" begin
    # Write your tests here.
    mdp = [0.0, 1.0, 2.0]
    mv = [1.0, 2.0, 4.0]
    t = raytrace(0.0, 0.0, 4.0, mdp, mv, ["all"])
    t == [
        (x = 4.0, t = 4.0, p = 1.0, l = 0, type = "refraction"),
        (
            x = 3.9999900280513696,
            t = 4.472127035819801,
            p = 0.8944267450394816,
            l = 2,
            type = "reflection",
        ),
        (
            x = 3.999990005738986,
            t = 4.037639086551873,
            p = 0.4192054730931121,
            l = 3,
            type = "reflection",
        ),
        (x = 4.0, t = 3.7320508075688776, p = 0.5, l = 2, type = "guide"),
        (x = 4.0, t = 3.802517076888147, p = 0.25, l = 3, type = "guide")
    ]
end
