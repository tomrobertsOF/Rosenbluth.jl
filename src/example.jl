using Revise
using Rosenbluth
using Rosenbluth.Models

import Rosenbluth: grow!, shrink!, atmosphere, positive_atmosphere, negative_atmosphere, size, atmosphericflattening

using BenchmarkTools

# @btime growshrinkgarm(SiteTree, 10, 1000)
# @btime Rosenbluth.pegarm(SiteTree, 100, 1000)

#s, w = Rosenbluth.growshrinkgarm(SiteTree, 10, 1000)

#w, s = Rosenbluth.sample(Models.SAW2D, 100, 10000, prune_enrich=true)

using Logging
debuglogger = ConsoleLogger(stderr, Logging.Debug)
Base.with_logger(debuglogger) do
    w, s = Rosenbluth.sample(Models.SAW2D, 20, 100000)
end

