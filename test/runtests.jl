using TestItems, TestItemRunner

@testsnippet RRModel begin
    struct RRModel <: RosenbluthSampleable
    end
end

@testsnippet RRModelSampleable begin
    function Rosenbluth.grow!(::RRModel)

    end
    function Rosenbluth.size(::RRModel)
        return 0
    end
    function Rosenbluth.atmosphere(::RRModel)
        return 1
    end

end

@testitem "Rosenbluth Sampleable cascades atmosphere" setup = [RRModel] begin
    let atmosphere_calls = 0
        function Rosenbluth.atmosphere(::RRModel)
            atmosphere_calls += 1
            return 1
        end

        @test Rosenbluth.positive_atmosphere(RRModel()) == 1
        @test atmosphere_calls == 1 # Atmosphere was called exactly once
    end
end

@testitem "Rosenbluth Sampleable defaults negative atmosphere" setup = [RRModel] begin
    @test Rosenbluth.negative_atmosphere(RRModel()) == 1
end

@testitem "sample throws when grow! not implemented" setup = [RRModel] begin
    function Rosenbluth.atmosphere(::RRModel)
        atmosphere_calls += 1
        return 1
    end
    function Rosenbluth.size(::RRModel)
        return 0
    end

    @test_throws ArgumentError sample(RRModel, 10, 10)
end
@testitem "sample throws when size not implemented" setup = [RRModel] begin
    function Rosenbluth.atmosphere(::RRModel)
        return 1
    end
    function Rosenbluth.grow!(::RRModel)

    end

    @test_throws ArgumentError sample(RRModel, 10, 10)
end
@testitem "sample throws when atmosphere not implemented" setup = [RRModel] begin
    function Rosenbluth.grow!(::RRModel)

    end
    function Rosenbluth.size(::RRModel)
        return 0
    end

    @test_throws ArgumentError sample(RRModel, 10, 10)
end

@testitem "sample does not throw when model is sampleable" setup = [RRModel, RRModelSampleable] begin
    # Stub out garm as the model is not actually well defined
    Rosenbluth.garm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true) = true

    @test_nowarn sample(RRModel, 10, 10)
end

@testitem "sample defaults to no pruning/enrichment" setup = [RRModel, RRModelSampleable] begin
    let garm_calls = 0, pegarm_calls = 0, growshrinkgarm_calls = 0
        function Rosenbluth.garm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            garm_calls += 1
            return true
        end
        function Rosenbluth.pegarm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            pegarm_calls += 1
            return true
        end
        function Rosenbluth.growshrinkgarm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            growshrinkgarm_calls += 1
            return true
        end

        sample(RRModel, 10, 10)

        @test garm_calls == 1
        @test pegarm_calls == 0
        @test growshrinkgarm_calls == 0
    end
end
@testitem "standard pruning/enrichment defaults to grow only if not shrink defined" setup = [RRModel, RRModelSampleable] begin
    let garm_calls = 0, pegarm_calls = 0, growshrinkgarm_calls = 0
        function Rosenbluth.garm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            garm_calls += 1
            return true
        end
        function Rosenbluth.pegarm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            pegarm_calls += 1
            return true
        end
        function Rosenbluth.growshrinkgarm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            growshrinkgarm_calls += 1
            return true
        end

        sample(RRModel, 10, 10, prune_enrich_method=:standard)

        @test garm_calls == 0
        @test pegarm_calls == 1
        @test growshrinkgarm_calls == 0
    end
end


@testitem "when shrink is defined, the sampler calls growshrinkgarm with standard p/e" setup = [RRModel, RRModelSampleable] begin
    function Rosenbluth.shrink!(::RRModel)
    end

    # stub out the actual samplers to check which one is called
    let garm_calls = 0, pegarm_calls = 0, growshrinkgarm_calls = 0
        function Rosenbluth.garm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            garm_calls += 1
            return true
        end
        function Rosenbluth.pegarm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            pegarm_calls += 1
            return true
        end
        function Rosenbluth.growshrinkgarm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            growshrinkgarm_calls += 1
            return true
        end

        sample(RRModel, 10, 10, prune_enrich_method=:standard)

        @test garm_calls == 0
        @test pegarm_calls == 0
        @test growshrinkgarm_calls == 1
    end
end

@testitem "Atmospheric Flattening throws when atm limits not defined" setup = [RRModel, RRModelSampleable] begin
    let garm_calls = 0, pegarm_calls = 0, growshrinkgarm_calls = 0
        function Rosenbluth.garm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            garm_calls += 1
            return true
        end
        function Rosenbluth.pegarm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            pegarm_calls += 1
            return true
        end
        function Rosenbluth.growshrinkgarm(::Type{RRModel}, max_size::Int, num_samples::Int; logging=true)
            growshrinkgarm_calls += 1
            return true
        end

        @test_throws ArgumentError sample(RRModel, 10, 10, prune_enrich_method=:atm_flat)
    end
end

@run_package_tests