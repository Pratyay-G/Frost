function choi_mpo_to_mps(mpo::MPO, ham)
    N = length(mpo)
    sites = siteinds(mpo)

    mps_tensors = ITensor[]
    mpo_tensors = ITensor[]

    for n in 1:N
        W = mpo[n]

        s = sites[n][2]
        sP = dag(prime(s))
        s_fused = combiner(s, sP; tags = tags(s), dir = dir(s))

        W *= s_fused
        push!(mps_tensors, W)

        newsP = sim(sP)
        local_id = dag(op(newsP, "Id"))
        s_fused_bis = replaceinds(s_fused, sP => newsP)
        T = (ham[n] * local_id) * s_fused_bis
        T *= prime(dag(s_fused_bis))
        push!(mpo_tensors, dag(T))
    end

    return MPS(mps_tensors), MPO(mpo_tensors)
end
