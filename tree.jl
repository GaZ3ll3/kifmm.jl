include("node.jl")

using PyPlot
using DataStructures

ϵ= 1e-12

type Tree
    dict::Dict{Int, Node}
    idMax
    root
    nSource
    nTarget
    
    rank
    maxLevel
    
    sourceTree
    normalTree
    targetTree
    
    center
    radius
    Tree() = new(Dict(), 0, -1, 0, 0, 0, 0, [], [], [], [Inf, Inf], [Inf, Inf])
end

function populate!(t::Tree, source, normal, target, nSource::Int, nTarget::Int, rank::Int, maxLevel::Int) 
    t.rank = rank
    t.nSource = nSource
    t.nTarget = nTarget
    t.maxLevel = 0
    t.sourceTree = source
    t.normalTree = normal
    t.targetTree = target
    
    getCenterRadius!(t, t.sourceTree)
    # set root's level and nodeIndex
    t.root = 1; t.dict[1] = Node(0, 0)
    t.idMax = t.root
    t.dict[t.root].nSource = t.nSource
    t.dict[t.root].nTarget = t.nTarget
    t.dict[t.root].center = t.center
    t.dict[t.root].radius = t.radius
    t.dict[t.root].sourceIndex = Array(1:nSource)
    t.dict[t.root].targetIndex = Array(1:nTarget)
    
    assignChildren!(t, t.root, maxLevel)
    buildTree!(t)
    checkTarget(t)
    
end

function getCenterRadius!(t::Tree, source)
    maxX = maximum(source[:, 1])
    minX = minimum(source[:, 1])
    maxY = maximum(source[:, 2])
    minY = minimum(source[:, 2])
    t.center = ([minX, minY] + [maxX, maxY])/2.0
    t.radius = ([maxX, maxY] - [minX, minY])/2.0
end

function assignChildren!(t::Tree, id::Int, maxLevel::Int)
    @assert(length(t.dict) >= id)
    @assert(t.root != -1)
    node = t.dict[id]
    # get target points
    for k = 1:node.nTarget
        push!(node.target, t.targetTree[node.targetIndex[k], :])
    end
    
    # no source points
    if (node.nSource == 0)
        node.isLeaf = true
        node.isEmpty = true
    end
    
    if (node.nSource != 0)
        for k = 1:node.nSource
            push!(node.source, t.sourceTree[node.sourceIndex[k], :])
            push!(node.normal, t.normalTree[node.sourceIndex[k], :])
        end
        
        # divide into subtrees
        if (node.nSource <= t.rank || node.nLevel == maxLevel)
            # at leaf
            node.isLeaf = true
            if t.maxLevel < node.nLevel
                t.maxLevel = node.nLevel
            end
        else
            # not a leaf, cut into children
            for k = 1:4
                t.idMax += 1
                push!(node.child, t.idMax)
                t.dict[t.idMax] = Node(node.nLevel + 1, k)
                t.dict[t.idMax].parent = id
                t.dict[t.idMax].center[1] = node.center[1] + (((k-1) & 1) - 0.5) * node.radius[1]
                t.dict[t.idMax].center[2] = node.center[2] + (((k-1) >> 1) - 0.5) * node.radius[2]
                t.dict[t.idMax].radius[1] = node.radius[1] * 0.5
                t.dict[t.idMax].radius[2] = node.radius[2] * 0.5
                t.dict[t.idMax].nSource = 0 # pre-allocate
                t.dict[t.idMax].nTarget = 0 # pre-allocate
            end
            
            # distribute all the particles, complete
            # sourceIndex, targetIndex, nSource, nTarget
            for k = 1:node.nSource
                # id is parent.
                parentSourceIndex = node.sourceIndex[k]
                firstBit = t.sourceTree[parentSourceIndex, 2] < node.center[2] ? 0:1
                secondBit = t.sourceTree[parentSourceIndex, 1] < node.center[1] ? 0:1
                childIndex = 2 * firstBit + secondBit
                cid = node.child[childIndex + 1]
                push!(t.dict[cid].sourceIndex, parentSourceIndex)
                t.dict[cid].nSource += 1
            end
            
            for k = 1:node.nTarget
                parentTargetIndex = node.targetIndex[k]
                firstBit = t.targetTree[parentTargetIndex, 2] < node.center[2] ? 0:1
                secondBit = t.targetTree[parentTargetIndex, 1] < node.center[1] ? 0:1
                childIndex = 2 * firstBit + secondBit
                cid = node.child[childIndex + 1]
                push!(t.dict[cid].targetIndex, parentTargetIndex)
                t.dict[cid].nTarget += 1
            end
            
            for k = 1: 4
                assignChildren!(t, node.child[k], maxLevel)
            end
        end
    end
end


function buildTree!(t::Tree)
    rootId = t.root
    d = t.dict
    rootCenter = d[rootId].center
    rootRadius = d[rootId].radius
    xMin, xMax = rootCenter[1] - rootRadius[1], rootCenter[1] + rootRadius[1]
    yMin, yMax = rootCenter[2] - rootRadius[2], rootCenter[2] + rootRadius[2]
    
    nodeQueue = Queue(Int)
    enqueue!(nodeQueue, rootId)
    while(!isempty(nodeQueue))
        frontId = dequeue!(nodeQueue)
        buildNode!(t, frontId, xMin, xMax, yMin, yMax)
        if length(d[frontId].child) != 0
            for k = 1:4
                enqueue!(nodeQueue, d[frontId].child[k])
            end
        end
    end
end


function buildNode!(t::Tree, id::Int, xMin, xMax, yMin, yMax)
    uSet = Set{Int}()
    vSet = Set{Int}()
    wSet = Set{Int}()
    xSet = Set{Int}()
    d = t.dict
    node = d[id]
    
    if (node.parent != -1) 
        parentNode = t.dict[node.parent]
        dx = node.radius[1]
        dy = node.radius[2]
        xStart = parentNode.center[1] - dx
        yStart = parentNode.center[2] - dy
        
        for col = -2:3
            for row = -2:3
                # in 6 x 6 grid
                currentX = xStart + 2 * col * dx
                currentY = yStart + 2 * row * dy
                if  (currentX >= xMin &&  currentX <= xMax && 
                     currentY >= yMin && currentY <= yMax && 
                    !(abs(node.center[1]-currentX) <ϵ  && abs(node.center[2] - currentY) < ϵ))
                    # within bounds and not equal to itself.
                    currentId = findNode(t, 1, [currentX, currentY])
                    adjacent = isAdjacent(t, id, currentId)
                    
                    currentNode = d[currentId]
                    if (currentNode.nLevel < node.nLevel)
                        # from coarse mesh
                        if (adjacent)
                            if currentNode.isLeaf
                                push!(uSet, currentId)
                            end
                        else
                            push!(xSet, currentId)
                        end
                    end
                    
                    if (currentNode.nLevel == node.nLevel)
                        if (!adjacent)
                            push!(vSet, currentId)
                        else
                            if (node.isLeaf)
                                rest = Queue(Int)
                                enqueue!(rest, currentId)
                                while(!isempty(rest))
                                    frontId = dequeue!(rest)
                                    frontNode = d[frontId]
                                    if (!isAdjacent(t, frontId, id))
                                        push!(wSet, frontId)
                                    else
                                        if (frontNode.isLeaf)
                                            push!(uSet, frontId)
                                        else
                                            
                                            for k = 1:4
                                                enqueue!(rest, frontNode.child[k])
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end 
        end  
    end   
    
    if (node.isLeaf)
        push!(uSet, id)
    end
    
    node.uList = Int[]
    node.vList = Int[]
    node.wList = Int[]
    node.xList = Int[]
    
    for u in uSet
        push!(node.uList, u)
    end
    
    for v in vSet
        push!(node.vList, v)
    end
    
    for w in wSet
        push!(node.wList, w)
    end
    
    for x in xSet
        push!(node.xList, x)
    end
    
    node.nUList = length(node.uList)
    node.nVList = length(node.vList)
    node.nWList = length(node.wList)
    node.nXList = length(node.xList)
end

function findNode(t::Tree, id::Int, center)
    # tree traversal
    node = t.dict[id]
    if norm(node.center - center, 1) < ϵ
        return id
    else
        if node.isLeaf
            return id
        else
            col = node.center[1] > center[1] ? 0 : 1
            row = node.center[2] > center[2] ? 0 : 1
            id = 2 * row + col
            return findNode(t, node.child[id + 1], center)
        end
    end         
end

function isAdjacent(t::Tree, aId, bId)
    nodeA = t.dict[aId]
    nodeB = t.dict[bId]
    diff = abs(nodeA.center - nodeB.center)
    r = nodeA.radius + nodeB.radius
   
    xAdjacent = (abs(diff[1]- r[1]) < ϵ) && (r[2] >= diff[2] - ϵ) # not far on the other direction
    yAdjacent = (abs(diff[2]- r[2]) < ϵ) && (r[1] >= diff[1] - ϵ) # not far on the other direction
    
    return xAdjacent || yAdjacent
end

    
    
function display(t::Tree, checknode=[])
    
    #scatter(t.sourceTree[:, 1], t.sourceTree[:, 2])
    lines = Any[]
    d = t.dict
    for i = 1 : length(d)
        cx = d[i].center[1]
        cy = d[i].center[2]
        rx = d[i].radius[1]
        ry = d[i].radius[2]
        xs = [cx - rx, cx + rx, cx + rx, cx - rx]
        ys = [cy - ry, cy - ry, cy + ry, cy + ry]
        push!(lines, collect(zip(xs, ys)))
    end
    line_seg = matplotlib[:collections][:LineCollection](lines)
    fig = figure("quadtree")
    ax = axes()
    ax[:add_collection](line_seg)
    #scatter(t.sourceTree[:, 1], t.sourceTree[:, 2], s = 3.8, alpha=0.5)
    #scatter(t.targetTree[:, 1], t.targetTree[:, 2], s = 5, alpha = 0.9)
    axis("image")
    if length(checknode) != 0
        for id in checknode
            checkList(t, id)
        end
    end
    
end

function checkList(t::Tree, nodeId)
    node = t.dict[nodeId]
    ucenterX = []
    ucenterY = []
    vcenterX = []
    vcenterY = []
    wcenterX = []
    wcenterY = []
    xcenterX = []
    xcenterY = []
    for id in node.uList
        nd = t.dict[id]
        push!(ucenterX, nd.center[1])
        push!(ucenterY, nd.center[2])
    end
    
    for id in node.vList
        nd = t.dict[id]
        push!(vcenterX, nd.center[1])
        push!(vcenterY, nd.center[2])
    end
    for id in node.wList
        nd = t.dict[id]
        push!(wcenterX, nd.center[1])
        push!(wcenterY, nd.center[2])
    end
    for id in node.xList
        nd = t.dict[id]
        push!(xcenterX, nd.center[1])
        push!(xcenterY, nd.center[2])
    end
    
    scatter(ucenterX, ucenterY, color = [1,0,0], s = 16.0, alpha = 0.7)
    scatter(vcenterX, vcenterY, color = [0,1.0,0.], s = 16.0, alpha = 0.7)
    scatter(wcenterX, wcenterY, color = [0,0,1], s = 16.0, alpha = 0.7)
    scatter(xcenterX, xcenterY, color = [0.5,0.5,0.5], s = 16.0, alpha = 0.7)
    scatter(node.center[1], node.center[2], color=[0,0,0], s = 16.0, alpha = 0.7)
end

function checkTarget(t::Tree)
    d, l = t.dict, length(t.dict)
    for i = 1:l
        node = d[i]
        ru = node.center + node.radius
        ld = node.center - node.radius
        nt = node.nTarget
        ti = node.targetIndex
        for j = 1: nt
            coordinate = t.targetTree[ti[j], :]
            if !((coordinate[1] <= ru[1] + ϵ && ϵ + coordinate[1] >= ld[1])
                && (coordinate[2] <= ru[2] + ϵ && ϵ + coordinate[2] >= ld[2]))
                @printf("(%6.4f, %6.4f) is outside of source box.\n", coordinate[1], coordinate[2])
            end
        end
    end
end
