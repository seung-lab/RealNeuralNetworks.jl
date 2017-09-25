using JLD
using TEASAR
using TEASAR.SWCs

data = load("/tmp/point_clouds.jld")
point_cloud = data["point_cloud"]
dbf = data["dbf"]

swc = skeletonize(point_cloud; dbf=dbf)
SWCs.save(swc, "/tmp/76880.swc")
