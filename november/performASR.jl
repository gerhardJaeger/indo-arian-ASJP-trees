using RCall;
using DataFrames;
using CSV;
using StatsBase;
using StatsFuns;
using Distributions;
using Random
using LinearAlgebra



R"library(ape)";

const posterior = CSV.read("output/gamma11.p",delim="\t");
const treetrace = CSV.read("output/gamma11.t",delim='\t');



#__convention__
#- nodes are numbered in postorder:
#   - 1..n: tips
#   - n+1..2n-1: internal nodes, mother before daughters

data = CSV.read("coded.csv")


function rootASR(itI::Int)
    q = zeros((4,4))
    counter = 0
    for i in 1:4
        for j in 1:4
            if i!=j
                counter += 1
                q[i,j] = posterior[itI,Symbol(string("rates[",counter,"]"))]
            end
        end
    end
    for i in 1:4
        q[i,i] = -sum(q[i,:]);
    end

    treeS = treetrace[itI,:tree];
    treeR = R"read.tree(text=$treeS)"
    tree = convert(Dict,treeR)
    nInternal = tree["Nnode"]
    n = length(tree["tip.label"])
    nNodes = nInternal+n
    loglik = zeros((nNodes,4));
    branchLengths = zeros((nNodes,nNodes))
    for i in 1:length(tree["edge.length"])
        b = tree["edge"][i,:]
        l = tree["edge.length"][i]
        branchLengths[b[1],b[2]]=l/1000
    end
    for nd in 1:n
        loglik[nd,:] .-= Inf
        nm = tree["tip.label"][nd]
        s = data.code[data.taxon.==nm][1]
        loglik[nd,s] = 0
    end
    for nd in nNodes:-1:(n+1)
        daughters = tree["edge"][tree["edge"][:,1].==nd,2]
        ll = branchLengths[nd,daughters]
        for j in 1:length(daughters)
            d = daughters[j]
            loglik[nd,:] += mapslices(logsumexp,
                log.(exp(ll[j].*q')).+loglik[d,:], dims=1)'
        end
    end

    rootLL = loglik[n+1,:]
    exp.(rootLL .- logsumexp(rootLL))
end

rootPosterior = zeros((nrow(posterior), 4))

for i in 1:nrow(posterior)
    rootPosterior[i,:] = rootASR(i)
end

rootPosterior = DataFrame(rootPosterior)

@rlibrary bayesplot
@rlibrary ggplot2



mcmc_areas(rootPosterior) + theme(text = element_text(size=40))
