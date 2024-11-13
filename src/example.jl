using Revise
using Rosenbluth
using Rosenbluth.Models

results = Rosenbluth.garmsample(Models.SiteTree, 100, 50000; prune_enrich_method=Rosenbluth.PruneEnrichMethod.ATMOSPHERIC_FLATTENING)
vec.(sum.(results, dims=(2,3)))[1]

# target_size = 10;
# num_tours = 1000;

# modeltype = Models.SiteTree

# weights = zeros(Float64, target_size, Rosenbluth.max_aplus(modeltype, target_size), Rosenbluth.max_aminus(modeltype, target_size));
# samples = zeros(Int, target_size, Rosenbluth.max_aplus(modeltype, target_size), Rosenbluth.max_aminus(modeltype, target_size));

# results = (weights, samples)

# target_weights = zeros(Float64, target_size)
# target_samples = zeros(Int, target_size)

# targets = (target_weights, target_samples)

# update_results! = (model, results, weight) -> begin
# bin_index = [Rosenbluth.size(model), Rosenbluth.positive_atmosphere(model), Rosenbluth.negative_atmosphere(model)]
# results[1][bin_index...] += weight
# results[2][bin_index...] += 1
# end;

# update_targets! = (model, targets, weight) -> begin
#         n = Rosenbluth.size(model)
#         targets[1][n] += weight
#         targets[2][n] += 1
# end;

# get_target_weight = (model, targets, started_tours) -> begin
#     n = Rosenbluth.size(model)
#     return targets[1][n] / started_tours
# end;

# for t in 1:num_tours
# Rosenbluth.mostgenerictour!(modeltype, target_size, t, results, targets,
#     update_results!,
#     update_targets!,
#     get_target_weight,
# )
# end

w, s = (vec.(sum.(results, dims=(2,3))))
# display(targets ./ num_tours)