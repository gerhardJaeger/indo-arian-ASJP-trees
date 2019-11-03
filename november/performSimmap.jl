const N = 100;


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



function sampleAnnotatedEdges(q::Array{Float64,2},itI::Int64)
    eq = exp(1000*q)[1,:]
    treeS = treetrace[itI,:tree]
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
        s = data[:code][data[:taxon].==nm][1]
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
    states = zeros(Int,nNodes)
    for l in 1:n
        nm = tree["tip.label"][l]
        s = data[:code][data[:taxon].==nm][1]
        states[l]=s
    end
    edges = tree["edge"]
    for nd in (n+1):nNodes
        logAsr = loglik[n+1,:].+log.(eq)
        if nd>n+1
            mother = edges[edges[:,2].==nd,1][1]
            logAsr += exp(branchLengths[mother,nd]*q)[states[mother],:]
        end
        asr = exp.(logAsr)
        asr /= sum(asr)
        s = sample(aweights(asr))
        states[nd]=s
    end
    nEdges = length(tree["edge.length"])
    annEdges = zeros((nEdges,3))
    for i in 1:nEdges
        annEdges[i,:] = [states[edges[i,1]],states[edges[i,2]],
                         tree["edge.length"][i]/1000]
    end
    return annEdges
end;

function sampleAnnotatedEdges(itI::Int64)
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
    ae = sampleAnnotatedEdges(q,itI)
    return ae,q
end;


function collectMutations(itI::Int)
    ae,q = sampleAnnotatedEdges(itI);
    function getMutations(q::Array{Float64,2},initial::Int,final::Int,
                          len::Float64)
        n = length(q[1,:])
        Lambda = -diag(q)
        P = copy(q)
        for i in 1:n
            P[i,i]=0
            P[i,:] /= sum(P[i,:])
        end
        simFinal = 0
        History = []
        while simFinal!=final
            History = []
            current = initial
            remaining = len
            while true
                #global remaining
                #global current
                #println("remaining =\t" * string(remaining))
                d = Exponential(1/Lambda[current])
                wt = rand(d,1)[1]
                #println("wt =\t\t" * string(wt))
                #println(wt >= remaining)
                if wt >= remaining
                    push!(History,(current,remaining))
                    #println("break")
                    break
                end
                newState = sample(aweights(P[current,:]))
                push!(History,(current,wt))
                remaining -= wt
                current = newState
            end
            simFinal = current
            #println(simFinal)
        end
        stateSequence = [x[1] for x in History]
        mutations = [stateSequence[1:end-1] stateSequence[2:end]]
    end;
    i=1;
    mutations = getMutations(q,
                             Int(ae[i,1]),
                             Int(ae[i,2]),
                             ae[i,3]);
    for i in 2:size(ae)[1]
        mutations = vcat(mutations,
                         getMutations(q,
                                      Int(ae[i,1]),
                                      Int(ae[i,2]),
                                      ae[i,3]))
    end
    mutations = DataFrame(mutations)
    names!(mutations,[:source,:target])
    mutations[:iteration] = itI
    return mutations
end;

nGenerations = length(posterior.Iteration);


generations = randperm(nGenerations)[1:N];

mutations = collectMutations(generations[1]);

for g in 2:N
    global mutations
    print(g,"\n")
    mutations = vcat(mutations,
                     collectMutations(generations[g]));
end

begin
    n = 4;
    mMatrix = zeros((n,n));
    for i in 1:length(mutations[:,1])
        a = mutations[:source][i]
        b = mutations[:target][i]
        mMatrix[a,b] += 1;
    end
end

mMatrix/N


CSV.write("simmapCounts.csv",mutations);
