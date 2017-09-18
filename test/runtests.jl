#---------------------------------------------------------------
#Script functionality - comment out for testing/importing

@assert length(ARGS) >= 2

input_filename = ARGS[1];
output_filename = ARGS[2];


penalty_fns = Dict(
  "alexs" => alexs_penalty,
  # "sato"  => literal_paper_penalty,
  # "sebs"  => sebastians_penalty,
  # "LMA"   => local_max_additive_penalty,
  "LMM"   => local_max_multiplicative_penalty
)

if length(ARGS) < 3
  penalty_fn = alexs_penalty
else
  penalty_fn = penalty_fns[ARGS[3]]
end


println("Reading points");
pts = read_points_file( input_filename );

@time edges, nodes, roots, radii, dests = Teasar( pts, penalty_fn );

# massaging the edges so they'll fit in a MAT file
edges = Int[edges[i][j] for i=1:length(edges),j=1:2];

println("Saving results...");
#My convention - nodes and edges
# are both indices into the rows of pts
# results = Dict(
#   "p" => pts,
#   "e" => edges,
#   "n" => nodes,
#   "root" => roots,
#   "rad" => radii,
#   "dest" => dests);

#Alex's convention - nodes are the coordinates themselves,
# edges are indices into the rows of nodes

# further massaging to fix Alex's convention
e2 = zeros(size(edges));
for i in eachindex(nodes)
  e2[ edges .== nodes[i] ] = i
end

results = Dict(
  "p" => pts,
  "e" => e2,
  "n" => pts[nodes,:],
  "root" => roots,
  "rad" => radii,
  "dest" => dests);

MAT.matwrite(output_filename, results)


