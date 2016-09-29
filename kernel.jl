include("tree.jl")

type KIFMM
    tree::Tree
    
    # surfaces variables
    nSurface
    stdUpEquivSurface
    stdUpCheckSurface
    stdDownEquivSurface
    stdDownCheckSurface
    stdNormal
    
    # kernel variables
    kernel
    
    # charge variables
    chargeTree
    
    KIFMM() = new(Tree(), 0, [], [], [], [], Nil, [])
end

function initialize!(k::KIFMM, nSurface, source, normal, target, charge, nSource, nTarget, rank, maxLevel, kernel::Function)
    k.tree = Tree()
    k.nSurface = nSurface
    # build tree
    t = k.tree
    populate!(t, source, normal, target, nSource, nTarget, rank, maxLevel)
    
    k.chargeTree = charge
    k.kernel = kernel
    
    k.stdUpEquivSurface = []
    k.stdUpCheckSurface = []
    k.stdDownEquivSurface = []
    k.stdDownCheckSurface = []
    k.stdNormal = []
    
    dr = 0.1
    
    getStdSurfaces!(√2 + dr, nSurface, k.stdUpEquivSurface)
    getStdSurfaces!(4.0 - √2 - 2 * dr, nSurface, k.stdUpCheckSurface)
    getStdSurfaces!(√2 + dr, nSurface, k.stdDownCheckSurface)
    getStdSurfaces!(4.0 - √2 - 2 * dr, nSurface, k.stdDownEquivSurface)
    getStdSurfaces!(1.0, nSurface, k.stdNormal)
end

function FMM!(k::KIFMM)
    upPass!(k, 1)
    potential = zeros(k.tree.nTarget)
    downPass!(k, 1, potential)
    return potential
end

function getStdSurfaces!(r::Float64, nSurface::Int, surface)
    θ = 2 * π/nSurface
    for l = 0:nSurface-1
        push!(surface, [r * cos(l*θ), r * sin(l*θ)])
    end
end

function getLocalSurface!(origin, center, radius, remote)
    nSurface = length(origin)
    for l = 1:nSurface
        push!(remote, [radius[1] * origin[l][1], radius[2] * origin[l][2]] + center)
    end
end

function upPass!(k::KIFMM, rootId)
    # not relevant to target, only related to source terms
    # S2M, M2M
    # bottom up, postorder traversal
    t = k.tree
    d = t.dict
    node = d[rootId]
    # at current node, initialize all surfaces
    
    # if there is no target points in this node then 
    # this is no Down-* surfaces
    getLocalSurface!(k.stdDownEquivSurface, node.center, node.radius, node.shiftedDownEquivalentSurface)
    getLocalSurface!(k.stdDownCheckSurface, node.center, node.radius, node.shiftedDownCheckSurface)
    node.downwardEquivalent = zeros(k.nSurface)

    # if there is no source points in this node then
    # there is no Up-* surfaces
    getLocalSurface!(k.stdUpEquivSurface, node.center, node.radius, node.shiftedUpEquivalentSurface)
    getLocalSurface!(k.stdUpCheckSurface, node.center, node.radius, node.shiftedUpCheckSurface)
    node.upwardEquivalent = zeros(k.nSurface)

    
    

    
    
    
    if (node.isLeaf)
        # at leaf, S2M
        getCharge!(k, rootId)
        # upward transfer to equivalent surface (target)
        contribution = getAccumulation(k, node.source, node.normal, node.shiftedUpCheckSurface, node.charge)
        node.upwardEquivalent = getInversion(k, node.shiftedUpEquivalentSurface, k.stdNormal, node.shiftedUpCheckSurface, contribution)
    else
        # at non-leaf, M2M
        contribution = zeros(k.nSurface)
        for l = 1:4
            upPass!(k, node.child[l])
            if (!d[node.child[l]].isEmpty) 
                contribution += getAccumulation(k, d[node.child[l]].shiftedUpEquivalentSurface, k.stdNormal, node.shiftedUpCheckSurface, d[node.child[l]].upwardEquivalent)
            end
        end
        node.upwardEquivalent = getInversion(k, node.shiftedUpEquivalentSurface, k.stdNormal, node.shiftedUpCheckSurface, contribution)
    end 
end

function downPass!(k::KIFMM, rootId, potential)
    # preorder traversal, first root, then children
    # M2L, L2L, L2T
    t = k.tree
    d = t.dict
    node = d[rootId]
    
 
    if (node.parent != -1)
        # V List
        contribution = zeros(k.nSurface)
        for l = 1: node.nVList
            if (!d[node.vList[l]].isEmpty)
                contribution += getAccumulation(k, d[node.vList[l]].shiftedUpEquivalentSurface, k.stdNormal, node.shiftedDownCheckSurface, d[node.vList[l]].upwardEquivalent)
            end
        end
        
        node.downwardEquivalent += getInversion(k, node.shiftedDownEquivalentSurface, k.stdNormal, node.shiftedDownCheckSurface, contribution)

        # X list
        contribution = zeros(k.nSurface)
        for l = 1: node.nXList
            if (!d[node.xList[l]].isEmpty)
                contribution += getAccumulation(k, d[node.xList[l]].shiftedUpEquivalentSurface, k.stdNormal, node.shiftedDownCheckSurface, d[node.xList[l]].upwardEquivalent)
            end
        end

        node.downwardEquivalent += getInversion(k, node.shiftedDownEquivalentSurface, k.stdNormal, node.shiftedDownCheckSurface, contribution)

        # L2L 
        # from parent
        parentNode = d[node.parent]
        contribution = getAccumulation(k, parentNode.shiftedDownEquivalentSurface, k.stdNormal, node.shiftedDownCheckSurface, parentNode.downwardEquivalent)
        node.downwardEquivalent += getInversion(k, node.shiftedDownEquivalentSurface, k.stdNormal, node.shiftedDownCheckSurface, contribution)        
    end
    
    if (node.isLeaf && node.nTarget != 0)
        node.potential = zeros(node.nTarget)
        # get U , W, Parent
        for l = 1:node.nUList
            # direct interactions
            getCharge!(k, node.uList[l])
            if (!d[node.uList[l]].isEmpty)
                node.potential += getAccumulation(k, d[node.uList[l]].source, d[node.uList[l]].normal, node.target, d[node.uList[l]].charge)
            end
        end
        
        for l = 1:node.nWList
            if (!d[node.wList[l]].isEmpty)
                node.potential += getAccumulation(k, d[node.wList[l]].shiftedUpEquivalentSurface, k.stdNormal, node.target, d[node.wList[l]].upwardEquivalent)
            end
        end
        
        # from down-* surface
        node.potential += getAccumulation(k, node.shiftedDownEquivalentSurface, k.stdNormal, node.target, node.downwardEquivalent)
        
        for l = 1:node.nTarget
            potential[node.targetIndex[l]] += node.potential[l]
        end
    end

    
    if !node.isLeaf
        for l = 1:4
            downPass!(k, node.child[l], potential)
        end
    end
    
end


function getCharge!(k::KIFMM, id)
    # lazy evaluation
    node = k.tree.dict[id]
    if (node.chargeComputed == true)
        
    else
        node.chargeComputed = true
        node.charge = zeros(node.nSource)
        for l = 1:node.nSource
            node.charge[l] = k.chargeTree[node.sourceIndex[l]]
        end
    end
end

function kernelEval(source, normal, target, f)
    sourceSize = length(source)
    targetSize = length(target)
    K = zeros(targetSize, sourceSize)
    for s = 1:sourceSize
        for t = 1: targetSize
            K[t, s] = f(source[s], normal[s], target[t])
        end
    end
    return K
end

function getAccumulation(k::KIFMM, source, normal, target, charge)
    K = kernelEval(source, normal, target, k.kernel)
    # K(target, source) * charge
    # density of length(target)
    # charge of lenth(source)
    return K * charge
end

function getInversion(k::KIFMM, source, normal, target, density)
    K = kernelEval(source, normal, target, k.kernel)
    return pinv(K)*density
end




