using Test
using StaticArrays
using MarineHydro

all_elements = [
    "named tuple of vectors" => (
                        (center=[0.0, 0.0, -1.0],),
                        (center=[1.0, 1.0, -2.0],
                         vertices=[-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0] .+ [1.0, 1.0, -2.0]',
                         normal=[0.0, 0.0, 1.0],
                         area=1.0,
                         radius=sqrt(2)/2)
                         ),
    "static elements" => (
                        MarineHydro.StaticElement(
                            SVector(0.0, 0.0, -1.0),
                            @SMatrix([-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0]) .+ SVector(0.0, 0.0, -1.0)',
                            SVector(0.0, 0.0, 1.0),
                            1.0,
                            sqrt(2)/2,
                        ),
                        MarineHydro.StaticElement(
                            SVector(1.0, 1.0, -2.0),
                            @SMatrix([-0.5 -0.5 0.0; 0.5 -0.5 0.0; 0.5 0.5 0.0; -0.5 0.5 0.0]) .+ SVector(1.0, 1.0, -2.0)',
                            SVector(0.0, 0.0, 1.0),
                            1.0,
                            sqrt(2)/2,
                            )
    )
]


@testset "Greens function for $elements_kind" for (elements_kind, elements) in all_elements

    e1 = elements[1]
    e2 = elements[2]

    @testset "Rankine" begin

        @test greens(Rankine(), e1, e2) ≈ 1/sqrt(1^2 + 1^2 + 1^2)  atol=1e-8
        @test gradient_greens(Rankine(), e1, e2, with_respect_to_first_variable=true) ≈ (e2.center - e1.center)/(sqrt(1^2 + 1^2 + 1^2))^3  atol=1e-8
        @test gradient_greens(Rankine(), e1, e2, with_respect_to_first_variable=false) ≈ (e1.center - e2.center)/(sqrt(1^2 + 1^2 + 1^2))^3  atol=1e-8

        @test integral(Rankine(), e1, e2) ≈ greens(Rankine(), e1, e2) * e2.area  atol=1e-1
        @test integral_gradient(Rankine(), e1, e2, with_respect_to_first_variable=true) ≈ gradient_greens(Rankine(), e1, e2, with_respect_to_first_variable=true) * e2.area  atol=1e-1
        @test integral_gradient(Rankine(), e1, e2, with_respect_to_first_variable=false) ≈ gradient_greens(Rankine(), e1, e2, with_respect_to_first_variable=false) * e2.area  atol=1e-1
    end

    @testset "ReflectedRankine" begin

        function vertical_reflection(x)
            return [x[1], x[2], -x[3]]
        end

        @testset "Alternative definitions of the reflected Rankine" begin
            import MarineHydro.free_surface_symmetry

            @test greens(Rankine(), free_surface_symmetry(e1), e2) ≈ greens(Rankine(), e1, free_surface_symmetry(e2))
            @test gradient_greens(Rankine(), free_surface_symmetry(e1), e2, with_respect_to_first_variable=true) ≈ vertical_reflection(gradient_greens(Rankine(), e1, free_surface_symmetry(e2), with_respect_to_first_variable=true))
            @test gradient_greens(Rankine(), free_surface_symmetry(e1), e2, with_respect_to_first_variable=false) ≈ vertical_reflection(gradient_greens(Rankine(), e1, free_surface_symmetry(e2), with_respect_to_first_variable=false))

            @test integral(Rankine(), free_surface_symmetry(e1), e2) ≈ integral(Rankine(), e1, free_surface_symmetry(e2))
            @test integral_gradient(Rankine(), free_surface_symmetry(e1), e2, with_respect_to_first_variable=true) ≈ vertical_reflection(integral_gradient(Rankine(), e1, free_surface_symmetry(e2), with_respect_to_first_variable=true))
            @test integral_gradient(Rankine(), free_surface_symmetry(e1), e2, with_respect_to_first_variable=false) ≈ vertical_reflection(integral_gradient(Rankine(), e1, free_surface_symmetry(e2), with_respect_to_first_variable=false))
        end


        @test greens(RankineReflected(), e1, e2) ≈ 1/sqrt(1^2 + 1^2 + 3^2)  atol=1e-8
        @test gradient_greens(RankineReflected(), e1, e2, with_respect_to_first_variable=true) ≈ (vertical_reflection(e2.center) - e1.center)/(sqrt(1^2 + 1^2 + 3^2))^3  atol=1e-8
        @test gradient_greens(RankineReflected(), e1, e2, with_respect_to_first_variable=false) ≈ vertical_reflection((e1.center - vertical_reflection(e2.center)))/(sqrt(1^2 + 1^2 + 3^2))^3  atol=1e-8
        @test integral(RankineReflected(), e1, e2) ≈ 1/sqrt(1^2 + 1^2 + 3^2) * e2.area  atol=1e-1
        @test integral_gradient(RankineReflected(), e1, e2, with_respect_to_first_variable=true) ≈ (vertical_reflection(e2.center) - e1.center)/(sqrt(1^2 + 1^2 + 3^2))^3 * e2.area atol=1e-1
        @test integral_gradient(RankineReflected(), e1, e2, with_respect_to_first_variable=false) ≈ vertical_reflection((e1.center - vertical_reflection(e2.center)))/(sqrt(1^2 + 1^2 + 3^2))^3 * e2.area atol=1e-1

    end

    @testset "Guével-Delhommeau" begin

        import MarineHydro._dimless_wave_term, MarineHydro._d_dimless_wave_term_dr
        # From cpt.Delhommeau().fortran_core.delhommeau_integrals.numerical_integration(1.0, -1.0, 1000)
        capy_ref = [-1.84003542, 1.76871985, -0.59117747, -1.01715699]
        @test _dimless_wave_term(1.0, -1.0) ≈ capy_ref[1] + im*capy_ref[2] atol=1e-4
        @test _d_dimless_wave_term_dr(1.0, -1.0) ≈ capy_ref[3] + im*capy_ref[4] atol=1e-4

        # From cpt.Delhommeau().fortran_core.green_wave.wave_part_infinite_depth(
        #              np.array([0, 0, -1]), np.array([1, 1, -2]), 1.0,
        #              1000, 0, np.array([]), np.array([]), np.array([]), cpt.Delhommeau().fortran_core.constants.low_freq)
        capy_g = -0.9289959399878619+0.17490912138602516im
        capy_gradg = [-0.06819068+0.12043414im, -0.06819068+0.12043414im, -0.32597325+0.17490912im]

        @test greens(ExactGuevelDelhommeau(), e1, e2, 1.0) ≈ capy_g atol=1e-4
        @test gradient_greens(ExactGuevelDelhommeau(), e1, e2, 1.0, with_respect_to_first_variable=true) ≈ capy_gradg atol=1e-4

    end

    @testset "GFWu, k=$k" for k in [0.1, 1.0, 10.0]
        @test greens(ExactGuevelDelhommeau(), e1, e2, k) ≈ greens(GFWu(), e1, e2, k)  atol=2e-2
        @test gradient_greens(ExactGuevelDelhommeau(), e1, e2, k) ≈ gradient_greens(GFWu(), e1, e2, k)  atol=2e-2
    end
end
