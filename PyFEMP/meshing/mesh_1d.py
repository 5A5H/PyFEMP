def msh_line(X0, X1, N, type='U1'):
    """
    msh_line(X0, X1, N, type='U1') -> XI, ELEM
    Function to return a regular mesh from a line.

    Input:
    X0           -> left coordinate
    X1           -> right coordinate
    N            -> regular mesh divisioning
    type = 'U1' -> specify mesh topology

    Output:
    XI           -> np.array with the nodal coordinates [...,[x,y],...]
    ELEM         -> np.array with the nodal indexes for each element [...,[n1,n2,n3,n4],...]
    """
    import numpy as np

    # generate nodes
    XI = np.linspace(X0, X1, N + 1)

    # generate elements
    if (type=='U1'): ELEM = np.array([np.array([nxi  ,nxi+1], dtype=np.uint) for nxi in range(N)], dtype=np.uint)

    else:
        raise NameError('Type :'+type+" not supported!")
        return 0, 0

    return XI, ELEM