# Class version of 1D Solver
import numpy as np
import matplotlib.pyplot as plt


class FEM_Simulation:
    '''
    Object that represents a 1D FEM Simulation.
    '''

    def __init__(self, Element):

        self.Element = Element

        # get initialize data from the element
        self.NoElementDim, \
        self.NoElementNodes, \
        self.ElementDofNames, \
        self.NoElementHistory, \
        self.ElementMaterialNames, \
        self.ElementPostNames = Element.Elmt_Init()

        self.NoElementMaterial = len(self.ElementMaterialNames)
        self.NoNodeDofs = len(self.ElementDofNames)

        # general program variables
        self.verbose = False
        self.verbose_system = True
        self.state = 0

        # general discretization variables
        self.time = 0.0
        self.dt = 1.0
        self.step = 1
        self.NoElements = 0
        self.NoNodes = 0
        self.NoDofs = 0
        self.XI = 0
        self.ELEM = 0

        # initialize fields for boundary conditions
        self.NBC = []
        self.EBC = []

        # prepare list of material parameter
        self.ElementMaterial = []

        # make some noise
        print("FEM Solver Instance Created")



    def Add_Mesh(self, NodesList, ElementConnectivity, verbose=False):
        ''' 
        Add_Mesh(self, NodesList, ElementConnectivity, verbose=False) -> void
        Sets a mesh based on a list of nodes and matrix of element connectivity.

        Input :
        NodeList            -> List of nodal coordinates [... , [x,y], ...]
        ElementConnectivity -> Matrix of nodal indexes per element [... , [n1, n2, n3], ...]
        '''
        # check input
        no_mesh_no, mesh_dim = NodesList.shape
        no_mesh_el, mesh_no_el = ElementConnectivity.shape
        if (verbose): print('Mesh NoNodes          : ',no_mesh_no)
        if (verbose): print('Mesh Dimension        : ',mesh_dim)
        if (verbose): print('Mesh NoElements       : ',no_mesh_el)
        if (verbose): print('Mesh Nodes per Element: ',mesh_no_el)

        if (self.NoElementDim != mesh_dim): raise NameError('Mesh dimension is not the same as elements.')
        if (self.NoElementNodes != mesh_no_el): raise NameError('Mesh is not compatible to element topology.')

        # process infos
        self.NoElements = no_mesh_el
        self.NoNodes = no_mesh_no
        self.NoDofs = no_mesh_no * self.NoNodeDofs
        if (verbose): print('Mesh Total Dofs       : ',self.NoDofs)

        self.XI = NodesList
        self.ELEM = ElementConnectivity

        if (self.verbose): print(' Finite Elemenmt Mesh Read!')
        self.state = 1




    def Add_Material(self, MaterialList, Option=None):
        '''Adds Material parameters as a list [....] for the next element without already specified material. With Option=All, all elements are set with the given list of parameters.'''
        if len(MaterialList) != self.NoElementMaterial:
            print(
                'Error: Number of material parameter does not fit element requirements!')
            print('       Requred parameters are :')
            print(*self.ElementMaterialNames, sep=", ")
            raise NameError('Error processing material parameters!')

        if Option == "All":
            self.ElementMaterial = []
            for i in range(self.NoElements):
                self.ElementMaterial.append(MaterialList)
            print(' Material set for All Elements')

        if Option == None:
            self.ElementMaterial.append(MaterialList)
            print(' Material set for Element %i' % len(self.ElementMaterial))

    
    
    def Add_EBC(self, NodeSelector, DofSelector, Value):
        '''Sets an essential boundary condition by NodeSelector, DofSelector, Value'''
        NodeList = self.SelectNodes(NodeSelector)
        AffectedDof = self.SelectDof(DofSelector)
        if AffectedDof >= self.NoNodeDofs:
            print("Error: Nodal degrees of freedom do not exceed %i" %
                  self.NoNodeDofs)
        for node in NodeList:
            self.EBC.append([node, AffectedDof, Value])

    
    
    def Add_NBC(self, NodeSelector, DofSelector, Value):
        '''Sets an essential boundary condition by NodeSelector, DofSelector, Value'''
        NodeList = self.SelectNodes(NodeSelector)
        AffectedDof = self.SelectDof(DofSelector)
        if AffectedDof >= self.NoNodeDofs:
            print("Error: Nodal degrees of freedom do not exceed %i" %
                  self.NoNodeDofs)
        for node in NodeList:
            self.NBC.append([node, AffectedDof, Value])

    def SelectDof(self, Input):
        '''Returns a single integer for the dof'''
        if isinstance(Input, int):
            return Input
        if isinstance(Input, str):
            for i, dofname in enumerate(self.ElementDofNames):
                if dofname == Input:
                    return i
        return 100


    def SelectNodes(self, Input):
        '''
        Returns a list containing the node number that fit the input.
        Input can be: 
            A single node index SelectNodes(0)
            A list of indexes node index SelectNodes([0,1,2])
            A conditional based on the dimension 1D: SelectNodes("x==0") 2D: SelectNodes("x==0 && y==0")
        '''
        Outlist = []

        # if input is a singe integer, check if there is a node for this integer and return it in a list
        if isinstance(Input, int):
            if Input <= self.NoElements:
                Outlist.append(Input)
            else:
                print('Error: %i is not a valid node number!' % Input)
                return

        # if input is a list, we check each entry for being an integer and proceed as before. if the entry is
        # an interger we append it to the output list
        if isinstance(Input, list):
            for i in range(len(Input)):
                if isinstance(Input[i], int):
                    if Input[i] <= self.NoElements:
                        Outlist.append(Input[i])
                    else:
                        print('Error: %i is not a valid node number!' %
                              Input[i])
                        return
                else:
                    print(
                        'Error: ', Input[i], " is not a valid node number! Integer required.")

        # if input is a string, it is supposed to be a conditional
        if isinstance(Input, str):

            # 1D - condition is x only
            if (self.NoElementDim==1):
                conditional = eval("lambda x: "+Input)
                Outlist = np.arange(self.NoNodes)[[conditional(x) for x in self.XI]]
            
            # 2D - condition is x and y
            elif (self.NoElementDim==2):
                conditional = eval("lambda x, y: "+Input)
                Outlist = np.arange(self.NoNodes)[[conditional(x,y) for x, y in self.XI]]

        return Outlist



    def Analysis(self):
        '''Enters into the Analysis phase. At least there must be finite elements and Materials'''
        if self.state < 1:
            self.state_report()
            return
        elif self.NoElements < 1:
            print('Error: No Elements! Use AddMesh')
            return
        elif len(self.ElementMaterial) != self.NoElements:
            print('Error: Not sufficent Material provided! Use AddMaterial')
            return
        self.h_t = np.zeros(self.NoElements * self.NoElementHistory)
        self.h_n = np.zeros(self.NoElements * self.NoElementHistory)

        print('Entering Analysis phase')
        if (self.verbose_system):
            print('---------------------------------')
            print('FE Setup Summary :')
            print('NoElementNodes   :', self.NoElementNodes)
            print('NoNodeDofs       :', self.NoNodeDofs)
            print('ElementDofNames  :', self.ElementDofNames)
            print('ElementPostNames :', self.ElementPostNames)
            print('NoElementHistory :', self.NoElementHistory)
            print('NoElements       :', self.NoElements)
            print('NoNodes          :', self.NoNodes)
            print('NoDofs           :', self.NoDofs)
            print('NoEssential BC   :', len(self.EBC))
            print('NoNatural BC     :', len(self.NBC))
            print('---------------------------------')

        self.state = 100

    def state_report(self):
        '''Gives hints to the user what to do next, based on a standard procedure.'''
        if self.state == 0:
            print('Input required: Call the AddMesh() function.')
        elif self.state == 1:
            print('state is 1')

    def NextStep(self, time=1, lambda_load=1):
        if self.state < 100:
            print('Error: Simulation has not entered analysis phase via Analysis().')
        if (self.verbose_system):
            print('\nCurrent Time : %f5.2' % time)
        self.dt = time - self.time
        if self.dt < 1e-16:
            print(
                'Error: Time given in NextStep is smaller than internal time: %f5.2' % time)
            return
        self.time = time
        self.step += 1
        self.h_n = self.h_t
        self.h_t = np.zeros(self.NoElements * self.NoElementHistory)
        self.lambda_load = lambda_load
        # apply EBC to DI
        for i in range(len(self.EBC)):
            ebc_dof_index = (self.EBC[i][0] * self.NoNodeDofs) + self.EBC[i][1]
            self.DI[ebc_dof_index] = self.lambda_load * self.EBC[i][2]
        return

    def CallElement(self, i):
        if i > self.NoElements:
            print('Error: Input exceeds number of elements. max is : %i8' %
                  self.NoElements)
        elmt_node_start = i * (self.NoElementNodes - 1)
        elmt_dof_start = i * (self.NoElementNodes - 1) * self.NoNodeDofs
        elmt_hist_start = i * (self.NoElementNodes - 1) * self.NoElementHistory

        Elmt_XI = self.XI[elmt_node_start:elmt_node_start+self.NoElementNodes]
        Elmt_UI = self.DI[elmt_dof_start:elmt_dof_start +
                          self.NoElementNodes*self.NoNodeDofs]
        Elmt_Hn = self.h_n[elmt_hist_start:elmt_hist_start +
                           self.NoElementHistory]
        Elmt_Ht = np.zeros(self.NoElementHistory)
        Elmt_Mat = self.ElementMaterial[i]
        if (self.verbose):
            print('Calling  : ', i)
            print('Elmt XI  : ', Elmt_XI)
            print('Elmt UI  : ', Elmt_UI)
            print('Elmt Mat : ', Elmt_Mat)
            print('Elmt Hn  : ', Elmt_Hn)
            print('Elmt Ht  : ', Elmt_Ht)
        # call element routine to get element vector / element matrix
        r_e, K_e = self.Element.Elmt_KS(
            Elmt_XI,
            Elmt_UI,
            Elmt_Hn,
            Elmt_Ht,
            Elmt_Mat,
            self.dt
        )
        self.h_t[elmt_hist_start:elmt_hist_start +
                 self.NoElementHistory] = Elmt_Ht
        if (self.verbose):
            print('Element Vector :', r_e)
            print('Element Matrix :')
            print(K_e)
        return r_e, K_e

    def Assemble(self):
        r = np.zeros(self.NoDofs)
        K = np.zeros((self.NoDofs, self.NoDofs))
        for e in range(self.NoElements):
            elmt_node_start = e * (self.NoElementNodes - 1)
            elmt_dof_start = e * (self.NoElementNodes - 1) * self.NoNodeDofs
            r_e, K_e = self.CallElement(e)

           # Assemble global vector and global matrix
            for i in range(self.NoElementNodes * self.NoNodeDofs):
                r[elmt_dof_start + i] += r_e[i]
                for j in range(self.NoElementNodes * self.NoNodeDofs):
                    K[elmt_dof_start + i, elmt_dof_start + j] += K_e[i][j]

        if (self.verbose):
            print('Global Vector :', r)
            print('Global Matrix :')
            print(K)
        return r, K

    def FormLinearSystem(self):
        r, K = self.Assemble()
        # Apply Natural Boundary Conditions
        r_NBC = np.zeros(self.NoDofs)
        for i in range(len(self.NBC)):
            nbc_dof_index = (self.NBC[i][0] * self.NoNodeDofs) + self.NBC[i][1]
            r_NBC[nbc_dof_index] -= self.lambda_load * self.NBC[i][2]
            if (self.verbose):
                print('Global Vector of NBCs :', r_NBC)

        # Apply Essential Boundary Conditions / Build reduction operator
        #    Build vector of Dof-Indexes to delete
        ebc_dof_indexes = []
        for i in range(len(self.EBC)):
            ebc_dof_index = (self.EBC[i][0] * self.NoNodeDofs) + self.EBC[i][1]
            ebc_dof_indexes.append(ebc_dof_index)
            self.DI[ebc_dof_index] = self.lambda_load * self.EBC[i][2]

        #    Build reduction operator
        I_full = np.diag(np.ones(self.NoDofs))
        I_red = np.delete(I_full, ebc_dof_indexes, 0)

        # Build reduced System of equations
        #    Sum up force vectors
        RHS = r + r_NBC

        # Reduce the System
        RHS = I_red.dot(RHS)
        LHS = I_red.dot(K.dot(I_red.T))
        if (self.verbose):
            print('Lin EqS RHS :', RHS)
            print('Lin EqS LHS :')
            print(LHS)

        return RHS, LHS, I_red

    def NewtonIteration(self):
        RHS, LHS, I_red = self.FormLinearSystem()
        Residual = np.sqrt(RHS.dot(RHS))
        if (self.verbose_system):
            print('      |R|    : %f10.8' % Residual)
        # Solve the linear System
        LHS_inv = np.linalg.inv(LHS)
        dDI = - LHS_inv.dot(RHS)
        if (self.verbose):
            print('Lin EqS Sol :', dDI)
        # Compute Norm of Solution Vector
        Norm_dDI = np.sqrt(dDI.dot(dDI))
        if (self.verbose_system):
            print('      |dDI|  : %f10.8' % Norm_dDI)

        # Update Global Vector of DOFs
        self.DI += dDI.dot(I_red)
        if (self.verbose):
            print('Current vector of unknowns :', self.DI)

        return np.sqrt(Residual*Norm_dDI)

    def Mod_Material(self, MaterialList, Element_I):
        ''' Modifies an already defined Material of element I'''
        if len(MaterialList) != self.NoElementMaterial:
            print(
                'Error: Number of material parameter does not fit element requirements!')
            print('       Requred parameters are :')
            print(*self.ElementMaterialNames, sep=", ")
            return
        if Element_I > len(self.ElementMaterial) - 1:
            print('Error: Material of Element %i is not defined yet.' % Element_I)
            return

        self.ElementMaterial[Element_I] = MaterialList
        return

    def NodalDof(self, NodeSelector, DofSelector):
        '''Return current DoF by NodeSelector, DofSelector'''
        NodeList = self.SelectNodes(NodeSelector)
        AffectedDof = self.SelectDof(DofSelector)
        if AffectedDof >= self.NoNodeDofs:
            print("Error: Nodal degrees of freedom do not exceed %i" %
                  self.NoNodeDofs)
        for node in NodeList:
            dof_index = (node * self.NoNodeDofs) + AffectedDof
            return self.DI[dof_index]

    def CallElementPost(self, i, PostName):
        if i > self.NoElements:
            print('Error: Input exceeds number of elements. max is : %i8' %
                  self.NoElements)
        elmt_node_start = i * (self.NoElementNodes - 1)
        elmt_dof_start = i * (self.NoElementNodes - 1) * self.NoNodeDofs
        elmt_hist_start = i * (self.NoElementNodes - 1) * self.NoElementHistory

        Elmt_XI = self.XI[elmt_node_start:elmt_node_start+self.NoElementNodes]
        Elmt_UI = self.DI[elmt_dof_start:elmt_dof_start +
                          self.NoElementNodes*self.NoNodeDofs]
        Elmt_Hn = self.h_n[elmt_hist_start:elmt_hist_start +
                           self.NoElementHistory]
        Elmt_Ht = np.zeros(self.NoElementHistory)
        Elmt_Mat = self.ElementMaterial[i]
        # call element routine to get element vector / element matrix
        X1, X2, P1, P2 = self.Element.Elmt_Post(
            Elmt_XI,
            Elmt_UI,
            Elmt_Hn,
            Elmt_Ht,
            Elmt_Mat,
            self.dt,
            PostName
        )
        return X1, X2, P1, P2

    def PostProcessing(self, PostName):
        XI = np.zeros(self.NoElementNodes*self.NoElements)
        PI = np.zeros(self.NoElementNodes*self.NoElements)
        for e in range(self.NoElements):
            X1, X2, P1, P2 = self.CallElementPost(e, PostName)
            XI[e*2+0] = X1
            XI[e*2+1] = X2
            PI[e*2+0] = P1
            PI[e*2+1] = P2
        return XI, PI
