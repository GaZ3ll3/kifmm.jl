type Node
    parent # parent Id of this node, -1 means root.
    child # child node Ids.
    uList # neighbour list (including itself), only leaf has
    vList # far field list
    wList # near field, only leaf has
    xList # far field list
    
    nLevel # level of this node
    nUList # number of neighbours
    nVList # number of V far field
    nWList # number of W near field
    nXList # number of X far field
    
    nodeIndex # numbering of this node in its parent's child
    
    center # center coordinate of this node
    radius # radius of this node
    
    nSource # source points number
    nTarget # target points  number
    
    sourceIndex # index for source for locating source
    targetIndex # index for target for locating target
    
    source # source list
    normal # normal list
    target # target list
    
    shiftedUpEquivalentSurface
    shiftedDownEquivalentSurface
    shiftedUpCheckSurface
    shiftedDownCheckSurface
    
    upwardEquivalent # associated with box
    downwardEquivalent # associated with box
    
    charge # associated with source
    potential # associated with target
    

    
    isLeaf
    isEmpty
    chargeComputed
    
    Node(nl, ni) = new(-1, Int[], Int[], Int[], Int[], Int[], nl, 0,0,0,0, ni, [Inf, Inf], [Inf, Inf], 0, 0, Int[], Int[], [],[],[], [],[],[], [], Float64[], Float64[], Float64[], Float64[], false, false, false)
end


    

    
    
    


