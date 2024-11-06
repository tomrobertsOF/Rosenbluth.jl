using Revise
using Rosenbluth
using Rosenbluth.Models

using BenchmarkTools

@btime growshrinkgarm(Models.SiteTree, 100, 1000, logging=false);
@btime pegarm(Models.SiteTree, 100, 1000, logging=false);
