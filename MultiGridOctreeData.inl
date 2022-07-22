
#include "Octree.h"

#define FORCE_UNIT_NORMALS 1
#define USE_DOT_RATIOS 1

#define MEMORY_ALLOCATOR_BLOCK_SIZE 1<<12
#define ITERATION_POWER 1.0/3

const float EPSILON=float(1e-6);
const float ROUND_EPS=float(1e-5);

// VertexData
long long EdgeIndex(const OctNode* node,const int& eIndex,const int& maxDepth,int idx[DIMENSION]){
    int o,i1,i2;
    int d,off[3];
    node->depthAndOffset(d,off);
    for(int i=0;i<DIMENSION;i++){idx[i]=BinaryNode<float>::CornerIndex(maxDepth+1,d+1,off[i]<<1,1);}
    Cube::FactorEdgeIndex(eIndex,o,i1,i2);
    switch(o){
        case 0:
            idx[1]=BinaryNode<float>::CornerIndex(maxDepth+1,d,off[1],i1);
            idx[2]=BinaryNode<float>::CornerIndex(maxDepth+1,d,off[2],i2);
            break;
        case 1:
            idx[0]=BinaryNode<float>::CornerIndex(maxDepth+1,d,off[0],i1);
            idx[2]=BinaryNode<float>::CornerIndex(maxDepth+1,d,off[2],i2);
            break;
        case 2:
            idx[0]=BinaryNode<float>::CornerIndex(maxDepth+1,d,off[0],i1);
            idx[1]=BinaryNode<float>::CornerIndex(maxDepth+1,d,off[1],i2);
            break;
    };
    return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}
long long VertexData::FaceIndex(const OctNode* node,const int& fIndex,const int& maxDepth,int idx[DIMENSION]){
    int dir,offset;
    Cube::FactorFaceIndex(fIndex,dir,offset);
    int d,o[3];
    node->depthAndOffset(d,o);
    for(int i=0;i<DIMENSION;++i)
       idx[i]=BinaryNode<float>::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);
    idx[dir]=BinaryNode<float>::CornerIndex(maxDepth+1,d,o[dir],offset);
    return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}

long long VertexData::CenterIndex(const OctNode *node, const int &maxDepth, int idx[DIMENSION]) {
    int d,o[3];
    node->depthAndOffset(d,o);
    for(int i=0;i<DIMENSION;++i){
        idx[i]=BinaryNode<float>::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);
    }
    return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}

long long VertexData::CornerIndex(const OctNode* node, const int& cIndex, const int& maxDepth, int idx[DIMENSION]){
    int x[DIMENSION];
    Cube::FactorCornerIndex(cIndex,x[0],x[1],x[2]);
    int d,o[3];
    node->depthAndOffset(d,o);
    for(int i=0;i<DIMENSION;++i){
        idx[i]=BinaryNode<float>::CornerIndex(maxDepth+1,d,o[i],x[i]);
    }
    return (long long)(idx[0]) | (long long)(idx[1])<<15 | (long long)(idx[2])<<30;
}


// SortedTreeNodes
SortedTreeNodes::SortedTreeNodes(void) {
    nodeCount=NULL;
    treeNodes=NULL;
    maxDepth=0;
}

SortedTreeNodes::~SortedTreeNodes(void) {
    if(nodeCount) delete [] nodeCount;
    if(treeNodes) delete [] treeNodes;
}

void SortedTreeNodes::set(OctNode& root, const bool& setIndex) {
    if(nodeCount) delete [] nodeCount;
    if(treeNodes) delete [] treeNodes;

    maxDepth=root.maxDepth()+1;
    // record the starting address of each depth, save in [0,maxDepth()+1]
    nodeCount=new int[maxDepth+1];
    // create an array contains pointer to all nodes
    treeNodes=new OctNode*[root.nodes()];

    OctNode* temp=root.nextNode();
    int i,cnt=0;
    while(temp){
        treeNodes[cnt++]=temp;
        temp=root.nextNode(temp);
    }
    // small depth node at front
    qsort(treeNodes,cnt,sizeof(const OctNode*),OctNode::CompareForwardPointerDepths);
    for(i=0;i<=maxDepth;++i){nodeCount[i]=0;}
    for(i=0;i<cnt;++i){
        if(setIndex) treeNodes[i]->nodeData.nodeIndex=i;
        nodeCount[treeNodes[i]->depth()+1]++;
    }
    for(i=1;i<=maxDepth;++i)
        nodeCount[i]+=nodeCount[i-1];
//    for(i=1;i<=maxDepth;++i)
//        printf("%d\n",nodeCount[i]);
}

// SortedTreeLeaves

SortedTreeLeaves::SortedTreeLeaves(void){
    treeLeaves=NULL;
}

SortedTreeLeaves::~SortedTreeLeaves(void){
    if(treeLeaves) delete [] treeLeaves;
}

void SortedTreeLeaves::set(OctNode& root){
    if(treeLeaves) delete [] treeLeaves;
    leafCount = root.leaves();
    treeLeaves=new OctNode*[leafCount];

    OctNode* temp=root.nextLeaf();
    int cnt=0;
    while(temp){
        treeLeaves[cnt++] = temp;
        temp = root.nextLeaf(temp);
    }
    qsort(treeLeaves,cnt,sizeof(const OctNode*),OctNode::CompareBackwardPointerDepths);
}

void SortedTreeLeaves::set(OctNode& root,const int& maxDepth){
    if(treeLeaves) delete [] treeLeaves;
    leafCount = root.leaves();
    treeLeaves=new OctNode*[leafCount];

    OctNode* temp=root.nextLeaf();
    int cnt=0;
    while(temp){
        if(temp->depth()<=maxDepth) treeLeaves[cnt++] = temp;
        temp = root.nextLeaf(temp);
    }
    qsort(treeLeaves,cnt,sizeof(const OctNode*),OctNode::CompareBackwardPointerDepths);
}

// Octree
template<int Degree>
float Octree<Degree>::GetLaplacian(const int idx[DIMENSION]) const
{
#if USE_DOT_RATIOS
    return float(fData.dotTable[idx[0]]*fData.dotTable[idx[1]]*fData.dotTable[idx[2]]*(fData.d2DotTable[idx[0]]+fData.d2DotTable[idx[1]]+fData.d2DotTable[idx[2]]));
#else // !USE_DOT_RATIOS
    double dot[3];
    dot[0]=fData.dotTable[idx[0]];
    dot[1]=fData.dotTable[idx[1]];
    dot[2]=fData.dotTable[idx[2]];
    return Real(
            fData.d2DotTable[idx[0]]*dot[1]*dot[2]+
            fData.d2DotTable[idx[1]]*dot[0]*dot[2]+
            fData.d2DotTable[idx[2]]*dot[0]*dot[1]
    );
#endif // USE_DOT_RATIOS
}

template <int Degree>
float Octree<Degree>::GetDivergence(const int idx[DIMENSION],const Point3D<float>& normal) const
{
#if USE_DOT_RATIOS
    double dot=fData.dotTable[idx[0]]*fData.dotTable[idx[1]]*fData.dotTable[idx[2]];
    return float(dot*(fData.dDotTable[idx[0]]*normal.coords[0]+fData.dDotTable[idx[1]]*normal.coords[1]+fData.dDotTable[idx[2]]*normal.coords[2]));
#else // !USE_DOT_RATIOS
    double dot[DIMENSION];
	dot[0]=fData.dotTable[idx[0]];
	dot[1]=fData.dotTable[idx[1]];
	dot[2]=fData.dotTable[idx[2]];
	return Real(
		fData.dDotTable[idx[0]]*normal.coords[0]*dot[1]*dot[2]+
		fData.dDotTable[idx[1]]*normal.coords[1]*dot[0]*dot[2]+
		fData.dDotTable[idx[2]]*normal.coords[2]*dot[0]*dot[1]
	);
#endif // !USE_DOT_RATIOS
}


// DivergenceFunction
template<int Degree>
void Octree<Degree>::DivergenceFunction::Function(OctNode *node1, OctNode *node2) {
    for(int i=0;i<3;++i)
        scratch[i]=index[i]+int(node1->off[i]);
    /**     SetLaplacianWeights():
     *      $scratch has been set to be the index of [node2][node1].
     *      $node1 is changing.                                      */
    node1->nodeData.value+=ot->GetDivergence(scratch,normal);
}

// LaplacianProjectionFunction
template<int Degree>
void Octree<Degree>::LaplacianProjectionFunction::Function(OctNode *node1, OctNode *node2) {
    for(int i=0;i<3;++i)
        scratch[i]=index[i]+int(node1->off[i]);
    /**     $scratch denote [this][node1].
     *      $node1->nodeData.value is still divergence.
     *      $this->value is the Solution of node1 to Lx=v.
     *
     *      Paper:
     *      Furthermore, there is an inherent multiresolution structure ,
     *      so we use an approach similar to the multigrid approach in [GKS02],
     *      solving the restriction Ld of L to the space spanned by the depth d functions (using a conjugate gradient solver)
     *      and projecting the fixed-depth solution back  to update the residual.*/
    node1->nodeData.value-=float(ot->GetLaplacian(scratch)*value);
}

template<int Degree>
int Octree<Degree>::LaplacianMatrixFunction::Function(OctNode *node1, OctNode *node2) {
    float temp;
    int d1=int(node1->d);
    int x1,y1,z1;
    x1=int(node1->off[0]);
    y1=int(node1->off[1]);
    z1=int(node1->off[2]);

    /**     this->d2 is $depth in GetFixedDepthLaplacian().
     *      $node1 start from &tree, $dDepth can't be smaller than 0.
     *      $dDepth > 0 means $node1 is not deep enough and
     *      $node1 will dive deeper and call this function again    */
    int dDepth=d2-d1;
    int d;
    d=(x2>>dDepth)-x1;
    if(d < 0) return 0;
    if(!dDepth){
        /**     $node1 is at $depth.
         *      Now need to make sure that every pair of two nodes
         *      is only counted once, because matrix is symmetric.
         *      for() in GetFixedDepthLaplacian() enumerate all nodes,
         *      and this function only count the node smaller than enumerated node. */
        if(!d){
            d=y2-y1;
            if(d < 0) return 0;
            if(!d){
                d=z2-z1;
                if(d < 0) return 0;
            }
        }
        for(int i=0;i<3;++i)
            scratch[i]=index[i]+(int)(node1->off[i]);

        temp=ot->GetLaplacian(scratch);
        if(node1==node2)
            temp/=2;

        if(fabs(temp)>EPSILON){
            rowElements[elementCount].Value=temp;
            /**     N = node1 offset in SortedNodes array with $depth   */
            rowElements[elementCount].N=node1->nodeData.nodeIndex-offset;
            ++elementCount;
        }
        return 0;
    }
    return 1;
}

template<int Degree>
int Octree<Degree>::RestrictedLaplacianMatrixFunction::Function(const OctNode* node1, const OctNode* node2){
    int d1,d2,off1[3],off2[3];
    node1->depthAndOffset(d1,off1);
    node2->depthAndOffset(d2,off2);
    int dDepth=d2-d1;
    int d;
    d=(off2[0]>>dDepth)-off1[0];
    if(d<0) return 0;

    if(!dDepth){
        if(!d){
            d=off2[1]-off1[1];
            if(d<0) return 0;
            if(!d){
                d=off2[2]-off1[2];
                if(d<0) return 0;
            }

            if(!OctNode::Overlap2(depth,offset,0.5,
                                  d1,off1,radius))
                return 0;

            for(int i=0;i<3;++i){
                scratch[i]=index[i]+BinaryNode<float>::Index(d1,off1[i]);
            }
            float temp=ot->GetLaplacian(scratch);
            if(node1==node2) temp/=2;

            if(fabs(temp)>EPSILON){
                rowElements[elementCount].Value=temp;
                rowElements[elementCount].N=node1->nodeData.nodeIndex;
                elementCount++;
            }
            return 0;
        }
    }
    return 1;
}

template<int Degree>
void Octree<Degree>::PointIndexValueFunction::Function(const OctNode *node) {
//    cnt++;
    int idx[DIMENSION];
    for(int i=0;i<3;++i)
        idx[i]=index[i]+int(node->off[i])*res2;

//    for(int i=0;i<3;++i){
//        printf("index[%d]:%d\n",i,index[i]);
//    }
//    for(int i=0;i<3;++i){
//        printf("off[%d]:%d\n",i,node->off[i]);
//    }
//    printf("%f,%f\n",node->nodeData.value,float(valueTables[idx[0]] * valueTables[idx[1]] * valueTables[idx[2]]));

    /**     Fo = Fo,x * Fo,y * Fo,z     */
    value += node->nodeData.value * float(valueTables[idx[0]] * valueTables[idx[1]] * valueTables[idx[2]]);
}

template<int Degree>
void Octree<Degree>::AdjacencyCountFunction::Function(const OctNode* node1, const OctNode* node2){
    adjacencyCount++;
}

// AdjacencySetFunction
template<int Degree>
void Octree<Degree>::AdjacencySetFunction::Function(const OctNode* node1, const OctNode* node2){
    adjacencies[adjacencyCount++]=node1->nodeData.nodeIndex;
}

// RefineFunction
template<int Degree>
void Octree<Degree>::RefineFunction::Function(OctNode* node1,const OctNode* node2){
    if(!node1->children && node1->depth()<this->depth)
        node1->initChildren();
//    cnt++;
}

template<int Degree>
int Octree<Degree>::GetFixedDepthLaplacian(SparseSymmetricMatrix<float>& matrix, const int& depth, const SortedTreeNodes& sNodes){
    LaplacianMatrixFunction mf;
    mf.ot=this;

    /**     the start of node with $depth       */
    mf.offset=sNodes.nodeCount[depth];
    /**     2 < myRadius << 3     */
    float myRadius=int(2*radius-ROUND_EPS)+ROUND_EPS;

    /**     matrix.rows = number of nodes at $depth     */
    matrix.Resize(sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth]);

    /**     call for max possible memory        */
    mf.rowElements=(MatrixEntry<float> *)malloc(sizeof(MatrixEntry<float>) * matrix.rows );
    for(int i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];++i){
        mf.elementCount=0;
        /**     mf.d2 is supposed to be $depth,
         *      every one is at the same depth  */
        mf.d2=int(sNodes.treeNodes[i]->d);

        mf.x2=int(sNodes.treeNodes[i]->off[0]);
        mf.y2=int(sNodes.treeNodes[i]->off[1]);
        mf.z2=int(sNodes.treeNodes[i]->off[2]);

        /**     input $temp node index information to $mf   */
        mf.index[0]=mf.x2*fData.res;
        mf.index[1]=mf.y2*fData.res;
        mf.index[2]=mf.z2*fData.res;

        /**     radius1:    1.500001
         *      radius2:    0.5         */
        OctNode::ProcessTerminatingNodeAdjacentNodes(sNodes.treeNodes[i],myRadius-float(0.5),
                                                     &tree,float(0.5),
                                                     &mf);
        /**     Set row [0, number of nodes with $depth] to be each round $elementCount */
        matrix.SetRowSize(i-sNodes.nodeCount[depth],mf.elementCount);
        memcpy(matrix.m_ppElements[i-sNodes.nodeCount[depth]],
               mf.rowElements, sizeof(MatrixEntry<float>) * mf.elementCount);
    }
    free(mf.rowElements);
}

template<int Degree>
int Octree<Degree>::GetRestrictedFixedDepthLaplacian(SparseSymmetricMatrix<float>& matrix,
                                                      const int* entries,const int& entryCount,
                                                      const OctNode* rNode,const float& radius,
                                                      const SortedTreeNodes& sNodes)
{
    int i;
    RestrictedLaplacianMatrixFunction mf;

    float myRadius=int(2*radius-ROUND_EPS)+ROUND_EPS;   // 4.001

    mf.ot=this;
    mf.radius=radius;   // 2.001

    /**     node with $depth-$startingDepth */
    rNode->depthAndOffset(mf.depth,mf.offset);
    matrix.Resize(entryCount);
    mf.rowElements=(MatrixEntry<float>*)malloc(sizeof(MatrixEntry<float>) * matrix.rows);

    /**     RestrictedLaplacianMatrixFunction::Function()
     *      need to use index as Laplacian matrix index  */
    for(i=0;i<entryCount;++i)
        sNodes.treeNodes[entries[i]]->nodeData.nodeIndex=i;
    for(i=0;i<entryCount;++i){
        mf.elementCount=0;
        for(int j=0;j<3;++j){
            mf.index[j]=int(sNodes.treeNodes[entries[i]]->off[j]) * fData.res;
        }
        OctNode::ProcessTerminatingNodeAdjacentNodes(sNodes.treeNodes[entries[i]],myRadius-float(0.5),
                                                     &tree,float(0.5),
                                                     &mf);
        matrix.SetRowSize(i,mf.elementCount);
        memcpy(matrix.m_ppElements[i],
               mf.rowElements,
               sizeof(MatrixEntry<float>) * mf.elementCount);

    }
    /**     reset back to original index    */
    for(i=0;i<entryCount;++i)
        sNodes.treeNodes[entries[i]]->nodeData.nodeIndex=entries[i];

    free(mf.rowElements);
    return 1;
}

template<int Degree>
int Octree<Degree>::SolveFixedDepthMatrix(const int& depth, const SortedTreeNodes& sNodes){
    int i,iter=0;
    Vector<double> V,Solution;
    /**     generating lower triangular matrix  */
    SparseSymmetricMatrix<float> matrix;

    /**     current process depth is fixed,
     *      extract all nodes' divergence with this depth    */
    V.Resize(sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth]);
    for(i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];++i) {
        V[i - sNodes.nodeCount[depth]] = sNodes.treeNodes[i]->nodeData.value;
//        printf("%f\n",sNodes.treeNodes[i]->nodeData.value);
    }



    /**     empty the allocator     */
    SparseSymmetricMatrix<float>::Allocator.rollBack();
    GetFixedDepthLaplacian(matrix,depth,sNodes);

//    for(i=0;i<matrix.rows;++i){
//        for(int j=0;j < matrix.rowSizes[i];++j){
//            printf("[%d][%d]=%f\n",i,matrix.m_ppElements[i][j].N,matrix.m_ppElements[i][j].Value);
//        }
//    }

    iter+=SparseSymmetricMatrix<float>::Solve(matrix,V,
                                              int(pow(matrix.rows,ITERATION_POWER)),
                                              Solution,double(EPSILON),1);

    for(i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];++i){
        sNodes.treeNodes[i]->nodeData.value = float(Solution[i-sNodes.nodeCount[depth]]);
//        printf("Solution:%f\n",sNodes.treeNodes[i]->nodeData.value);
    }

    float myRadius,myRadius1,myRadius2;
    /**     myRadius = 1.0001   */
    myRadius=float(radius+ROUND_EPS-0.5);

    myRadius1=myRadius2=radius;

    /**     myRadius1 = 1.4449
     *      myRadius2 = 1.0001  */
    myRadius1=float( int(2*myRadius1+1-ROUND_EPS) )/2-ROUND_EPS;
    myRadius2=float( int(2*myRadius2  -ROUND_EPS) )/2+ROUND_EPS;

    /**     myRadius  = 1.0001 * width
     *      myRadius1 = 1.4449  * width
     *      myRadius2 = 1.0001  * width     */
    myRadius /=(1<<depth);
    myRadius1/=(1<<depth);
    myRadius2/=(1<<depth);

    float dx,dy,dz;
    int x1,x2,y1,y2,z1,z2;
    if(depth<sNodes.maxDepth-1){
        LaplacianProjectionFunction pf;
        OctNode* node1,* node2;
        pf.ot=this;
        int idx1,idx2;
        int offset=sNodes.nodeCount[depth];

        /**     update the residual of node1 by
         *      projecting node2 solution      */
        for(i=0;i<matrix.rows;++i){
            idx1=i;
            node1=sNodes.treeNodes[idx1+offset];
            if(!node1->children) continue;
            x1=int(node1->off[0]);
            y1=int(node1->off[1]);
            z1=int(node1->off[2]);
            for(int j=0;j<matrix.rowSizes[i];++j){
                idx2=matrix.m_ppElements[i][j].N;
                node2=sNodes.treeNodes[idx2+offset];
                x2=int(node2->off[0]);
                y2=int(node2->off[1]);
                z2=int(node2->off[2]);

                pf.value=Solution[idx2];
                pf.index[0]= x2 * fData.res;
                pf.index[1]= y2 * fData.res;
                pf.index[2]= z2 * fData.res;
                dx=float(x2-x1)/(1<<depth);
                dy=float(y2-y1)/(1<<depth);
                dz=float(z2-z1)/(1<<depth);

                if(fabs(dx) < myRadius && fabs(dy) < myRadius && fabs(dz) < myRadius) {
                    /**     When node1 and node2 is neighbor    */
                    node1->processNodeNodes(node2, &pf, 0);
                } else {
                    /**     radius1:    1.4449  * width
                     *      radius2:    1.0001  * width
                     *      width2:     width
                     *      (dx, dy, dz) is the vector(node1->node2)    */
                    OctNode::ProcessNodeAdjacentNodes(dx,dy,dz,
                                                      node2,myRadius1,
                                                      node1,myRadius2, float(1.0/(1<<depth)),
                                                      &pf,0);
                }
            }
        }

        /**     update the residual of node2 by
         *      projecting node1 solution      */
        for(i=0;i<matrix.rows;++i){
            idx1=i;
            node1=sNodes.treeNodes[idx1+offset];
            x1=int(node1->off[0]);
            y1=int(node1->off[1]);
            z1=int(node1->off[2]);
            pf.value=Solution[idx1];
            pf.index[0]= x1 * fData.res;
            pf.index[1]= y1 * fData.res;
            pf.index[2]= z1 * fData.res;
            for(int j=0;j<matrix.rowSizes[i];++j){
                idx2=matrix.m_ppElements[i][j].N;
                node2=sNodes.treeNodes[idx2+offset];
                if(idx1!=idx2 && node2->children){
                    x2=int(node2->off[0]);
                    y2=int(node2->off[1]);
                    z2=int(node2->off[2]);
                    dx=float(x1-x2)/(1<<depth);
                    dy=float(y1-y2)/(1<<depth);
                    dz=float(z1-z2)/(1<<depth);
                    if(fabs(dx) < myRadius && fabs(dy) < myRadius && fabs(dz) < myRadius) {
                        node2->processNodeNodes(node1,&pf,0);
                    }else{
                        OctNode::ProcessNodeAdjacentNodes(dx,dy,dz,
                                                          node1,myRadius1,
                                                          node2,myRadius2,float(1.0/(1<<depth)),
                                                          &pf,0);
                    }
                }
            }
        }
    }

    return iter;
}

template<int Degree>
int Octree<Degree>::SolveFixedDepthMatrix(const int& depth, const int& startingDepth, const SortedTreeNodes& sNodes){
    if(startingDepth>=depth)
        return SolveFixedDepthMatrix(depth,sNodes);

    int i,j,d,iter=0;
    SparseSymmetricMatrix<float> matrix;
    AdjacencySetFunction asf;
    AdjacencyCountFunction acf;
    Vector<float> Values;
    Vector<double> SubValues,SubSolution;

    Values.Resize(sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth]);

    for(i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];++i){
        Values[i-sNodes.nodeCount[depth]]=sNodes.treeNodes[i]->nodeData.value;
        sNodes.treeNodes[i]->nodeData.value=0;
    }

    float myRadius,myRadius2;

    myRadius=2*radius-float(0.5);   //2.5
    myRadius=int(myRadius-ROUND_EPS)+ROUND_EPS; //2.0001
    myRadius2=float(radius+ROUND_EPS-0.5);  //1.0001

    int idx1,idx2;
    /**     start from depth-startingDepth  */
    d=depth-startingDepth;
    for(i=sNodes.nodeCount[d];i<sNodes.nodeCount[d+1];++i) {
        idx1=i;
        OctNode* temp;
        /**     Author: Get all of the entries associated to the subspace   */
        acf.adjacencyCount=0;
        temp=sNodes.treeNodes[idx1]->nextNode();
        /**     Count the sub-node at $depth in the sub-tree of node1 */
        while(temp) {
            if (temp->depth() == depth) {
                acf.Function(temp, temp);
                temp = sNodes.treeNodes[idx1]->nextBranch(temp);
            } else {
                temp = sNodes.treeNodes[idx1]->nextNode(temp);
            }
        }
        for(j=sNodes.nodeCount[d];j<sNodes.nodeCount[d+1];++j){
            idx2=j;
            if(idx1==idx2) continue;
            /**     Count how many sub-nodes at $depth in the sub-tree of node2 is close to node1  */
            OctNode::ProcessFixedDepthNodeAdjacentNodes(sNodes.treeNodes[idx1],float(0.5),
                                                        sNodes.treeNodes[idx2],myRadius,
                                                        depth,&acf);
        }
        if(!acf.adjacencyCount) continue;
        asf.adjacencies=new int[acf.adjacencyCount];
        asf.adjacencyCount=0;

        temp=sNodes.treeNodes[idx1]->nextNode();
        /**     Record the global index of the nodes in sub-tree of node1 at $depth  */
        while(temp) {
            if (temp->depth() == depth) {
                asf.Function(temp,temp);
                temp = sNodes.treeNodes[idx1]->nextBranch(temp);
            } else {
                temp = sNodes.treeNodes[idx1]->nextNode(temp);
            }
        }

        for(j=sNodes.nodeCount[d];j<sNodes.nodeCount[d+1];++j){
            idx2=j;
            if(idx1==idx2) continue;
            OctNode::ProcessFixedDepthNodeAdjacentNodes(sNodes.treeNodes[idx1],float(0.5),
                                                        sNodes.treeNodes[idx2],myRadius,
                                                        depth,&asf);
        }

        /**     Author: Get the associated vector   */
        SubValues.Resize(asf.adjacencyCount);
        for(j=0;j<asf.adjacencyCount;++j){
            SubValues[j]=Values[asf.adjacencies[j]-sNodes.nodeCount[depth]];
        }
        SubSolution.Resize(asf.adjacencyCount);
        for(j=0;j<asf.adjacencyCount;++j) {
            SubSolution[j]=sNodes.treeNodes[asf.adjacencies[j]]->nodeData.value;
        }

        /**     Author: Get the associated matrix   */
        SparseSymmetricMatrix<float>::Allocator.rollBack();
        /**     myRadius:   2.001       */
        GetRestrictedFixedDepthLaplacian(matrix,asf.adjacencies,asf.adjacencyCount,
                                         sNodes.treeNodes[i],myRadius,sNodes);

        /**     Author: Solve the matrix    */
        iter+=SparseSymmetricMatrix<float>::Solve(matrix,SubValues,
                                                  int(pow(matrix.rows,ITERATION_POWER)),
                                                  SubSolution,double(EPSILON),0);

        LaplacianProjectionFunction lpf;
        lpf.ot=this;

        /**     Author: Update the solution for all nodes in the sub-tree    */
        for(j=0;j<asf.adjacencyCount;++j){
            temp=sNodes.treeNodes[asf.adjacencies[j]];
            /**     Initial :
             *      depth > depth - startingDepth    */
            while(temp->depth() > sNodes.treeNodes[i]->depth())
                temp=temp->parent;
            idx2=temp->nodeData.nodeIndex;
            if(idx2 >= sNodes.treeNodes[i]->nodeData.nodeIndex) // = idx1 ?
                /**     Solution saved in nodeData.value    */
                sNodes.treeNodes[asf.adjacencies[j]]->nodeData.value=float(SubSolution[j]);
        }

        /**     Author: Update the values in the next depth    */
        int x1,x2,y1,y2,z1,z2;
        float dx,dy,dz;
        float sRadius=radius/(1<<depth);
        /**     leaf don't update next depth    */
        if(depth < sNodes.maxDepth-1){
            OctNode* node1,* node2;
            /**     use asf.adjacencies to transform the local index to global index    */
            for(j=0;j<matrix.rows;++j){
                idx1=asf.adjacencies[j];
                node1=sNodes.treeNodes[idx1];
                if(!node1->children) continue;

                x1=int(node1->off[0]);
                y1=int(node1->off[1]);
                z1=int(node1->off[2]);

                for(int k=0;k<matrix.rowSizes[j];++k){
                    idx2=asf.adjacencies[matrix.m_ppElements[j][k].N];
                    node2=sNodes.treeNodes[idx2];
                    temp=node2;
                    while(temp->depth()>d) // d=depth-startingDepth;
                        temp=temp->parent;
                    /**     they don't have same ancestor    */
                    if(temp!=sNodes.treeNodes[i]) continue;

                    x2=int(node2->off[0]);
                    y2=int(node2->off[1]);
                    z2=int(node2->off[2]);

                    lpf.value=float(SubSolution[matrix.m_ppElements[j][k].N]);
                    lpf.index[0]= x2 * fData.res;
                    lpf.index[1]= y2 * fData.res;
                    lpf.index[2]= z2 * fData.res;

                    dx=float(x2-x1)/(1<<depth);
                    dy=float(y2-y1)/(1<<depth);
                    dz=float(z2-z1)/(1<<depth);

                    if(fabs(dx)<myRadius2 && fabs(dy)<myRadius2 && fabs(dz)<myRadius2 ) {
                        node1->processNodeNodes(node2, &lpf, 0);
                    }else{
                        OctNode::ProcessNodeAdjacentNodes(dx,dy,dz,
                                                          node2,sRadius,
                                                          node1,sRadius,float(1.0/(1<<depth)),
                                                          &lpf,0);
                    }
                }
            }

            for(j=0;j<matrix.rows;++j){
                idx1=asf.adjacencies[j];
                node1=sNodes.treeNodes[idx1];
                temp=node1;
                while(temp->depth()>d)
                    temp=temp->parent;
                if(temp!=sNodes.treeNodes[i]) continue;

                x1=int(node1->off[0]);
                y1=int(node1->off[1]);
                z1=int(node1->off[2]);

                lpf.value = float(SubSolution[j]);
                lpf.index[0]= x1 * fData.res;
                lpf.index[1]= y1 * fData.res;
                lpf.index[2]= z1 * fData.res;

                for(int k=0;k<matrix.rowSizes[j];++k){
                    idx2=asf.adjacencies[matrix.m_ppElements[j][k].N];
                    node2=sNodes.treeNodes[idx2];
                    if(!node2->children)
                        continue;

                    if(idx1!=idx2) {
                        x2 = int(node2->off[0]);
                        y2 = int(node2->off[1]);
                        z2 = int(node2->off[2]);

                        dx = float(x1 - x2) / (1 << depth);
                        dy = float(y1 - y2) / (1 << depth);
                        dz = float(z1 - z2) / (1 << depth);

                        if (fabs(dx) < myRadius2 && fabs(dy) < myRadius2 && fabs(dz) < myRadius2) {
                            node2->processNodeNodes(node1, &lpf, 0);
                        } else {
                            OctNode::ProcessNodeAdjacentNodes(dx, dy, dz,
                                                              node1, sRadius,
                                                              node2, sRadius, float(1.0 / (1 << depth)),
                                                              &lpf, 0);
                        }
                    }
                }
            }
        }

        delete [] asf.adjacencies;
    }
    return iter;

}

template<int Degree>
void Octree<Degree>::SetIsoSurfaceCorners(const float& isoValue,const int subdivisionDepth,const int& fullDepthIso){
    int i,j;
    unordered_map<long long, float> values;
    float cornerValues[Cube::CORNERS];

    PointIndexValueFunction cf;
    Point3D<float> position;
    long long key;

    SortedTreeNodes *sNodes=new SortedTreeNodes();
    sNodes->set(tree,0);

    OctNode* temp=tree.nextNode();
    while(temp){
        temp->nodeData.isoNode=NULL;
        temp=tree.nextNode(temp);
    }

    NodeData::UseIndex=0;

    cf.valueTables=fData.valueTables;
    cf.res2=fData.res2;

    for(i=0;i<sNodes->nodeCount[subdivisionDepth];++i){
        temp=sNodes->treeNodes[i];
        if(!temp->children){
            for(j=0;j<Cube::CORNERS;++j){
                cf.value=0;
                /**     Get the corner index in valueTable  */
                VertexData::CornerIndex(temp,j,fData.depth,cf.index);
                for(int k=0;k<DIMENSION;++k)
                    position.coords[k]=BinaryNode<float>::CornerIndexPosition(cf.index[k],fData.depth+1);

                OctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);

                cornerValues[j]=cf.value;
            }

            /**     check if this node is leaf and if there's surface in this node  */
            if(temp->depth() < fData.depth || MarchingCubes::HasRoots(cornerValues,isoValue)){
                temp->nodeData.isoNode = new IsoNodeData();
                memcpy(temp->nodeData.isoNode->cornerValues,
                       cornerValues,
                       sizeof(float) * Cube::CORNERS);
            }
        }
    }

    for(i=sNodes->nodeCount[subdivisionDepth];i<sNodes->nodeCount[subdivisionDepth+1];++i){
        temp=sNodes->treeNodes[i]->nextLeaf();
        while(temp){
            /**     assign the corner value of a leaf node with depth >= subdivisionDepth   */
            for(j=0;j<Cube::CORNERS;++j){
                key=VertexData::CornerIndex(temp,j,fData.depth,cf.index);
                if(values.find(key)!=values.end())
                    cornerValues[j]=values[key];
                else{
                    cf.value=0;
                    for(int k=0;k<DIMENSION;++k)
                        position.coords[k]=BinaryNode<float>::CornerIndexPosition(cf.index[k],fData.depth+1);
                    OctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);
                    values[key]=cf.value;
                    cornerValues[j]=cf.value;
                }
            }
            if(temp->depth()<fData.depth || MarchingCubes::HasRoots(cornerValues,isoValue)){
                temp->nodeData.isoNode=new IsoNodeData();
                memcpy(temp->nodeData.isoNode->cornerValues,
                       cornerValues,
                       sizeof(float) * Cube::CORNERS);
            }
            temp=sNodes->treeNodes[i]->nextLeaf(temp);
        }
        values.clear();
    }

    delete sNodes;

    if(subdivisionDepth){
        PreValidate(isoValue,fData.depth,subdivisionDepth);
    }

    temp=tree.nextLeaf();
    while(temp){
        if(subdivisionDepth) Validate(temp,isoValue,fData.depth,fullDepthIso,subdivisionDepth);
        else                 Validate(temp,isoValue,fData.depth,fullDepthIso);
        temp=tree.nextLeaf(temp);
    }

}

template<int Degree>
int Octree<Degree>::IsBoundaryFace(const OctNode* node,const int& faceIndex,const int& subdivideDepth){
    int dir,offset,d,o[3],idx;

    if(subdivideDepth<0) return 0;
    if(node->depth()<=subdivideDepth) return 1;
    Cube::FactorFaceIndex(faceIndex,dir,offset);
    node->depthAndOffset(d,o);

    idx=(int(o[dir])<<1) + (offset<<1);
    return !( idx % (2 << (int(node->depth()) - subdivideDepth)) );
}

template<int Degree>
void Octree<Degree>::PreValidate(OctNode* node,const float& isoValue,const int& maxDepth,const int& subdivideDepth){
    int sub=0;
    for(int i=0;i<Cube::NEIGHBORS && !sub;++i){
        OctNode* neighbor=node->faceNeighbor(i);
        if(neighbor && neighbor->children && IsBoundaryFace(node,i,subdivideDepth))
            sub=1;
    }

    if(sub){
        Subdivide(node,isoValue,maxDepth);
        for(int i=0;i<Cube::NEIGHBORS;++i){
            if(IsBoundaryFace(node,i,subdivideDepth)){
                OctNode* neighbor=node->faceNeighbor(i);
                while(neighbor && !neighbor->children){
                    PreValidate(neighbor,isoValue,maxDepth,subdivideDepth);
                    neighbor=node->faceNeighbor(i);
                }
            }
        }
    }
}

template<int Degree>
void Octree<Degree>::PreValidate(const float& isoValue,const int& maxDepth,const int& subdivideDepth){
    OctNode* temp;

    temp=tree.nextLeaf();
    while(temp){
        PreValidate(temp,isoValue,maxDepth,subdivideDepth);
        temp=tree.nextLeaf(temp);
    }
}

template<int Degree>
void Octree<Degree>::Validate(OctNode* node,const float& isoValue,const int& maxDepth,const int& fullDepthIso,const int& subdivideDepth){
    int i,sub=0;
    OctNode* treeNode=node;
    OctNode* neighbor;
    if(node->depth()>=maxDepth){return;}

    // Check if full-depth extraction is enabled and we have an iso-node that is not at maximum depth
    if(!sub && fullDepthIso && MarchingCubes::HasRoots(node->nodeData.isoNode->cornerValues,isoValue))
        sub=1;

    // Check if the node has faces that are ambiguous and are adjacent to finer neighbors
    for(i=0;i<Cube::NEIGHBORS && !sub;i++){
        neighbor=treeNode->faceNeighbor(i);
        if(neighbor && neighbor->children)
            if(MarchingCubes::IsAmbiguous(node->nodeData.isoNode->cornerValues,isoValue,i) ||
            IsBoundaryFace(node,i,subdivideDepth))
            sub=1;
    }

    // Check if the node has edges with more than one root
    for(i=0;i<Cube::EDGES && !sub;i++)
        if(EdgeRootCount(node,i,isoValue,maxDepth)>1)
            sub=1;

    for(i=0;i<Cube::NEIGHBORS && !sub;i++){
        neighbor=node->faceNeighbor(i);
        if(	neighbor && neighbor->children &&
               !MarchingCubes::HasRoots(node->nodeData.isoNode->cornerValues,isoValue,i) &&
               InteriorFaceRootCount(neighbor,Cube::FaceReflectFaceIndex(i,i),isoValue,maxDepth))
            sub=1;
    }
    if(sub){
        Subdivide(node,isoValue,maxDepth);
        for(i=0;i<Cube::NEIGHBORS;i++){
            neighbor=treeNode->faceNeighbor(i);
            if(neighbor && !neighbor->children)
                Validate(neighbor,isoValue,maxDepth,fullDepthIso,subdivideDepth);
        }
        for(i=0;i<Cube::EDGES;i++){
            neighbor=treeNode->edgeNeighbor(i);
            if(neighbor && !neighbor->children)
                Validate(neighbor,isoValue,maxDepth,fullDepthIso,subdivideDepth);
        }
        for(i=0;i<Cube::CORNERS;i++)
            if(!node->children[i].children)
                Validate(&node->children[i],isoValue,maxDepth,fullDepthIso,subdivideDepth);
    }
}

template<int Degree>
void Octree<Degree>::Validate(OctNode* node,const float& isoValue,const int& maxDepth,const int& fullDepthIso){
    int i,sub=0;
    OctNode* leaf=node;
    OctNode* neighbor;

    /**     leaf at maximum depth    */
    if(node->depth()>=maxDepth) return;

    /**     Author: Check if full-depth extraction is enabled
     *      and we have an iso-node that is not at maximum depth     */
    if(!sub && fullDepthIso && MarchingCubes::HasRoots(node->nodeData.isoNode->cornerValues,isoValue))
        sub=1;

    /**     Author: Check if the node has faces that are ambiguous
     *      and are adjacent to finer neighbors                      */
    for(i=0;i < Cube::NEIGHBORS && !sub;++i){
        neighbor=leaf->faceNeighbor(i);
        if(neighbor && neighbor->children &&
           MarchingCubes::IsAmbiguous(node->nodeData.isoNode->cornerValues,isoValue,i) )
            sub=1;
    }

    /**     Author: Check if the node has edges with more than one root */
    for(i=0;i<Cube::EDGES && !sub;++i){
        if(EdgeRootCount(node,i,isoValue,maxDepth)>1)
            sub=1;
    }

    for(i=0;i<Cube::NEIGHBORS && !sub;++i){
        neighbor=node->faceNeighbor(i);
        if(neighbor && neighbor->children &&
           !MarchingCubes::HasRoots(node->nodeData.isoNode->cornerValues,isoValue,i) &&
           InteriorFaceRootCount(neighbor,Cube::FaceReflectFaceIndex(i,i) , isoValue,maxDepth))
            sub=1;
    }

    if(sub){
        Subdivide(node,isoValue,maxDepth);
        for(i=0;i<Cube::NEIGHBORS;++i){
            neighbor=leaf->faceNeighbor(i);
            if(neighbor && !neighbor->children)
                Validate(neighbor,isoValue,maxDepth,fullDepthIso);
        }
        for(i=0;i<Cube::EDGES;++i){
            neighbor=leaf->edgeNeighbor(i);
            if(neighbor && !neighbor->children)
                Validate(neighbor,isoValue,maxDepth,fullDepthIso);
        }
        for(i=0;i<Cube::CORNERS;++i){
            if(!node->children[i].children)
                Validate(&neighbor->children[i],isoValue,maxDepth,fullDepthIso);
        }
    }

}

template<int Degree>
void Octree<Degree>::Subdivide(OctNode* node,const float& isoValue,const int& maxDepth){
    int i,j,c[4];
    float value,cornerValues[Cube::CORNERS];
    float cornerValues2[Cube::CORNERS][Cube::CORNERS];
    PointIndexValueFunction cf;
    Point3D<float> position;
    cf.valueTables=fData.valueTables;
    cf.res2=fData.res2;

    for(i=0;i<Cube::CORNERS;++i)
        cornerValues[i]=node->nodeData.isoNode->cornerValues[i];

    node->initChildren();
    for(i=0;i<Cube::CORNERS;++i)
        node->children[i].nodeData.isoNode=NULL;

    for(i=0;i<Cube::CORNERS;++i)
        cornerValues2[i][i]=cornerValues[i];

    /**     Set center corner   */
    cf.value=0;
    VertexData::CenterIndex(node,maxDepth,cf.index);
    for(i=0;i<DIMENSION;++i)
        /**     transform from index to length in res2 = 1<<(maxDepth+1)     */
        position.coords[i]=BinaryNode<float>::CornerIndexPosition(cf.index[i],maxDepth+1);
    OctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);
    /**     the value at center position   */
    value=cf.value;
    for(i=0;i<Cube::CORNERS;++i)
        cornerValues2[i][Cube::AntipodalCornerIndex(i)]=value;

    for(i=0;i<Cube::NEIGHBORS;++i){
        int dir,offset,e;
        Cube::FactorFaceIndex(i,dir,offset);
        cf.value=0;
        VertexData::FaceIndex(node,i,maxDepth,cf.index);
        for(i=0;i<DIMENSION;++i)
            position.coords[i]=BinaryNode<float>::CornerIndexPosition(cf.index[i],maxDepth+1);
        OctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);
        value=cf.value;
        Cube::FaceCorners(i,c[0],c[1],c[2],c[3]);
        e=Cube::EdgeIndex(dir,0,0);
        for(j=0;j<4;++j){
            cornerValues2[c[j]][Cube::EdgeReflectCornerIndex(c[j],e)]=value;
        }
    }

    for(i=0;i<Cube::EDGES;++i){
        int o,i1,i2,f;
        Cube::FactorEdgeIndex(i,o,i1,i2);
        cf.value=0;
        VertexData::EdgeIndex(node,i,maxDepth,cf.index);
        for(i=0;i<DIMENSION;++i)
            position.coords[i]=BinaryNode<float>::CornerIndexPosition(cf.index[i],maxDepth+1);
        OctNode::ProcessPointAdjacentNodes(position,&tree,radius,&cf);
        value=cf.value;
        Cube::EdgeIndex(i,c[0],c[1]);
        f=Cube::FaceIndex(o,0);
        for(j=0;j<2;++j)
            cornerValues2[c[j]][Cube::FaceReflectCornerIndex(c[j],f)]=value;
    }

    for(i=0;i<Cube::CORNERS;++i){
        if(node->children[i].depth()<maxDepth || MarchingCubes::HasRoots(cornerValues2[i],isoValue)){
            node->children[i].nodeData.isoNode=new IsoNodeData();
            memcpy(node->children[i].nodeData.isoNode->cornerValues,
                   cornerValues2[i],
                   sizeof(float) * Cube::CORNERS);
        }
    }

    if(node->nodeData.isoNode){
        node->nodeData.isoNode=NULL;
    }

}

template<int Degree>
int Octree<Degree>::InteriorFaceRootCount(const OctNode* node,const int& faceIndex,const float& isoValue,const int maxDepth){
    int c1,c2,e1,e2,dir,off,cnt=0;
    int corners[Cube::CORNERS/2];
    if(node->children){
        Cube::FaceCorners(faceIndex,corners[0],corners[1],corners[2],corners[3]);
        Cube::FactorFaceIndex(faceIndex,dir,off);
        c1=corners[0];
        c2=corners[3];
        /**     get edge on face    */
        switch(dir){
            case 0:
                e1=Cube::EdgeIndex(1,off,1);
                e2=Cube::EdgeIndex(2,off,1);
                break;
            case 1:
                e1=Cube::EdgeIndex(0,off,1);
                e2=Cube::EdgeIndex(2,1,off);
                break;
            case 2:
                e1=Cube::EdgeIndex(0,1,off);
                e2=Cube::EdgeIndex(1,1,off);
                break;
        };
        /**     check inner face root   */
        cnt+= EdgeRootCount(&node->children[c1],e1,isoValue,maxDepth)+
              EdgeRootCount(&node->children[c1],e2,isoValue,maxDepth);
        switch(dir){
            case 0:
                e1=Cube::EdgeIndex(1,off,0);
                e2=Cube::EdgeIndex(2,off,0);
                break;
            case 1:
                e1=Cube::EdgeIndex(0,off,0);
                e2=Cube::EdgeIndex(2,0,off);
                break;
            case 2:
                e1=Cube::EdgeIndex(0,0,off);
                e2=Cube::EdgeIndex(1,0,off);
                break;
        };
        cnt+= EdgeRootCount(&node->children[c2],e1,isoValue,maxDepth)+
              EdgeRootCount(&node->children[c2],e2,isoValue,maxDepth);
        for(int i=0;i<Cube::CORNERS/2;++i){
            if(node->children[corners[i]].children){
                cnt+=InteriorFaceRootCount(&node->children[corners[i]],faceIndex,isoValue,maxDepth);
            }
        }
    }
    return cnt;
}


template<int Degree>
int Octree<Degree>::EdgeRootCount(const OctNode* node,const int& edgeIndex,const float& isoValue,const int& maxDepth){
    int f1,f2,c1,c2;
    const OctNode* temp;
    Cube::FacesAdjacentToEdge(edgeIndex,f1,f2);

    int eIndex=edgeIndex;
    const OctNode* finest=node;
    if(node->depth()<maxDepth){
        temp=node->faceNeighbor(f1);
        if(temp && temp->children){
            finest=temp;
            /**     the original edge's new index on f1 faceNeighbor cube  */
            eIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f1);
        }else{
            temp=node->faceNeighbor(f2);
            if(temp && temp->children){
                finest=temp;
                eIndex=Cube::FaceReflectEdgeIndex(edgeIndex,f2);
            }else{
                temp=node->edgeNeighbor(edgeIndex);
                if(temp && temp->children){
                    finest=temp;
                    /**     the original edge's new index on edgeNeighbor cube  */
                    eIndex=Cube::EdgeReflectEdgeIndex(edgeIndex);
                }
            }
        }
    }

    Cube::EdgeCorners(eIndex,c1,c2);
    if(finest->children) {
        return EdgeRootCount(&finest->children[c1], eIndex, isoValue, maxDepth) +
               EdgeRootCount(&finest->children[c2], eIndex, isoValue, maxDepth);
    } else{
        if(!finest->nodeData.isoNode){
            return 0;
        }
        float* cValues=finest->nodeData.isoNode->cornerValues;
        if( (cValues[c1]<=isoValue && cValues[c2]<=isoValue) || (cValues[c1]>=isoValue && cValues[c2]>=isoValue) )
            return 0;
    }
    return 1;
}


template<int Degree>
int Octree<Degree>::NonLinearUpdateWeightContribution(OctNode* node, const Point3D<float>& position) {
    int i,j,k;
    /**     Get the neighbor of $node,
     *      initChildren is called and all elements is valid.   */
    OctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);

    /**     dx[0] denote dim x,
     *      dx[1] denote dim y,
     *      dx[2] denote dim z.                 */
    double x,dx[DIMENSION][3];
    double width;
    /**     Get the center and width of $node   */
    Point3D<float> center;
    float w;
    node->centerAndWidth(center,w);
    width=w;

    for(i=0;i<DIMENSION;++i){
        /**     x denote the distance to index 0 in this dimension,
         *      x bigger, $position is more close to neighbor with index 0 in this dimension */
        x=(center.coords[i]-position.coords[i]-width)/width;
        dx[i][0]=1.125+1.500*x+0.500*x*x;

        /**     x denote the distance to index 1 in this dimension,
         *      abs(x) bigger, less contribution to neighbor with index 1 in this dimension */
        x=(center.coords[i]-position.coords[i])/width;
        dx[i][1]=0.750        -      x*x;

        dx[i][2]=1.0-dx[i][0]-dx[i][1];
    }

    for(i=0;i<3;++i){
        for(j=0;j<3;++j){
            for(k=0;k<3;++k){
                if(neighbors.neighbors[i][j][k])
                    neighbors.neighbors[i][j][k]->nodeData.centerWeightContribution += float(dx[0][i] * dx[1][j] * dx[2][k]);
            }
        }
    }
    return 0;
}


template<int Degree>
float Octree<Degree>::NonLinearGetSampleWeight(OctNode* node, const Point3D<float>& position) {
    float weight=0;
    double x,dx[DIMENSION][3];
    int i,j,k;
    OctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);
    double width;
    Point3D<float> center;
    float w;
    node->centerAndWidth(center,w);
    width=w;

    for(i=0;i<DIMENSION;++i){
        /**     x denote the distance to index 0 in this dimension,
         *      x bigger, $position is more close to neighbor with index 0 in this dimension */
        x=(center.coords[i]-position.coords[i]-width)/width;
        dx[i][0]=1.125+1.500*x+0.500*x*x;

        /**     x denote the distance to index 1 in this dimension,
         *      abs(x) bigger, less contribution to neighbor with index 1 in this dimension */
        x=(center.coords[i]-position.coords[i])/width;
        dx[i][1]=0.750        -      x*x;

        dx[i][2]=1.0-dx[i][0]-dx[i][1];
    }

    for(i=0;i<DIMENSION;++i){
        for(j=0;j<DIMENSION;++j){
            for(k=0;k<DIMENSION;++k){
                if(neighbors.neighbors[i][j][k])
                    weight+=float(dx[0][i]*dx[1][j]*dx[2][k]*
                                  neighbors.neighbors[i][j][k]->nodeData.centerWeightContribution);
            }
        }
    }
    return float(1.0/weight);
}

template<int Degree>
void Octree<Degree>::NonLinearGetSampleDepthAndWeight(OctNode* node, const Point3D<float>& position,
                                                      const float& samplesPerNode,
                                                      float& depth,
                                                      float& weight)
{
    OctNode* temp=node;
    weight=float(1.0)/NonLinearGetSampleWeight(temp,position);
    if(weight>=samplesPerNode+1)
        depth=float(temp->depth()+log( weight / (samplesPerNode+1) ) / log(double(1<<(DIMENSION-1) ) ) );
    else{
        float oldAlpha,newAlpha;
        oldAlpha=newAlpha=weight;
        while(newAlpha<(samplesPerNode+1) && temp->parent){
            temp=temp->parent;
            oldAlpha=newAlpha;
            newAlpha=float(1.0)/NonLinearGetSampleWeight(temp,position);
        }
        depth=float(temp->depth()+log( newAlpha / (samplesPerNode+1) ) / log( newAlpha / oldAlpha ) );
    }
    weight=float(pow(double(1<<(DIMENSION-1) ), -double(depth) ) );
}

template<int Degree>
int Octree<Degree>::NonLinearSplatOrientedPoint(OctNode* node, const Point3D<float>& position, const Point3D<float>& normal) {
    double x,dxdydz,dx[DIMENSION][3];
    int i,j,k;
    OctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);
    double width;
    Point3D<float> center;
    float w;
    node->centerAndWidth(center,w);
    width=w;

    for(i=0;i<DIMENSION;++i){
        /**     x denote the distance to index 0 in this dimension,
         *      x bigger, $position is more close to neighbor with index 0 in this dimension */
        x=(center.coords[i]-position.coords[i]-width)/width;
        dx[i][0]=1.125+1.500*x+0.500*x*x;

        /**     x denote the distance to index 1 in this dimension,
         *      abs(x) bigger, less contribution to neighbor with index 1 in this dimension */
        x=(center.coords[i]-position.coords[i])/width;
        dx[i][1]=0.750        -      x*x;

        dx[i][2]=1.0-dx[i][0]-dx[i][1];
    }

    for(i=0;i<DIMENSION;++i){
        for(j=0;j<DIMENSION;++j){
            for(k=0;k<DIMENSION;++k){
                if(neighbors.neighbors[i][j][k]){
                    dxdydz=dx[0][i]*dx[1][j]*dx[2][k];
                    int idx=neighbors.neighbors[i][j][k]->nodeData.nodeIndex;
                    if(idx<0){
                        Point3D<float> n;
                        n.coords[0]=n.coords[1]=n.coords[2]=0;
                        idx=neighbors.neighbors[i][j][k]->nodeData.nodeIndex=int(normals->size());
                        normals->push_back(n);
                    }
                    (*normals)[idx].coords[0]+=float(normal.coords[0]*dxdydz);
                    (*normals)[idx].coords[1]+=float(normal.coords[1]*dxdydz);
                    (*normals)[idx].coords[2]+=float(normal.coords[2]*dxdydz);
//                    printf("%lf\n",float(normal.coords[0]*dxdydz));
                }
            }
        }
    }
    return 0;
}

template<int Degree>
void Octree<Degree>::NonLinearSplatOrientedPoint(const Point3D<float>& position,
                                                 const Point3D<float>& normal,
                                                 const int& splatDepth,
                                                 const float& samplesPerNode,
                                                 const int& minDepth,
                                                 const int& maxDepth)
{
    double dx;
    Point3D<float> n;
    OctNode* temp;
    int i;
    double width;
    Point3D<float> myCenter;
    float myWidth;
    myCenter.coords[0]=myCenter.coords[1]=myCenter.coords[2]=float(0.5);
    myWidth=float(1.0);

    temp=&tree;
    /**     reach the nearest node to $position point at splatDepth */
    while(temp->depth()<splatDepth){
        if(!temp->children){
            printf("Can't reach splatDepth.\nOctree<Degree>::NonLinearSplatOrientedPoint error.\n");
            return;
        }
        int cIndex=OctNode::CornerIndex(myCenter,position);
        temp=&temp->children[cIndex];
        myWidth/=2;
        for(i=0;i<DIMENSION;++i){
            if(cIndex & (1<<i)) myCenter.coords[i]+=myWidth/2;
            else                myCenter.coords[i]-=myWidth/2;
        }
    }

    float alpha,newDepth;
    NonLinearGetSampleDepthAndWeight(temp,position,samplesPerNode,newDepth,alpha);

    if(newDepth<minDepth)
        newDepth=float(minDepth);
    if(newDepth>maxDepth)
        newDepth=float(maxDepth);
    int topDepth=int(ceil(newDepth));

    dx=1.0-(topDepth-newDepth);
    if(topDepth<=minDepth){
        topDepth=minDepth;
        dx=1;
    }
    else if(topDepth>maxDepth){
        topDepth=maxDepth;
        dx=1;
    }

    /**     make $temp be the nearest node at $topDepth to $position point  */
    while(temp->depth()>topDepth) temp=temp->parent;
    while(temp->depth()<topDepth) {
        if(!temp->children) temp->initChildren();
        int cIndex=OctNode::CornerIndex(myCenter,position);
        temp=&temp->children[cIndex];
        myWidth/=2;
        for(i=0;i<DIMENSION;++i){
            if(cIndex & (1<<i)) myCenter.coords[i]+=myWidth/2;
            else                myCenter.coords[i]-=myWidth/2;
        }
    }

    /**     scale the original normal   */
    width=1.0/(1<<temp->depth());
    for(i=0;i<DIMENSION;++i) {
        n.coords[i] = normal.coords[i] * alpha
                      / float(pow(width, 3))
                      * float(dx);
    }

    /**     update the contribution to normals[temp] at $topDepth   */
    NonLinearSplatOrientedPoint(temp,position,n);
    /**     In case that newDepth is between [topDepth-1, topDepth],
     *      update the contribution to normals[temp->parent] at $topDepth-1 */
    if(fabs(1.0-dx)>EPSILON){
        dx=float(1.0-dx);
        temp=temp->parent;
        width=1.0/(1<<temp->depth());

        for(i=0;i<DIMENSION;++i) {
            n.coords[i] = normal.coords[i] * alpha
                          / float(pow(width,3))
                          * float(dx);
        }
        NonLinearSplatOrientedPoint(temp,position,n);
    }
}

template<int Degree>
int Octree<Degree>::HasNormals(OctNode *node, const float &epsilon) {
    int hasNormals=0;
    const int& idx=node->nodeData.nodeIndex;
    const Point3D<float>& normal=(*normals)[idx];
    if(idx >= 0 && Length(normal) > epsilon)
        hasNormals=1;
    if(node->children)
        for(int i=0;i<Cube::CORNERS && !hasNormals; ++i){
            hasNormals|=HasNormals(&node->children[i],epsilon);
        }
    return hasNormals;
}

// Octree public function
template<int Degree>
Octree<Degree>::Octree(void) {
    radius=0;
    postNormalSmooth=0;
}

template<int Degree>
void Octree<Degree>::setFunctionData(const PPolynomial<Degree>& ReconstructionFunction,
                     const int& maxDepth,
                     const int& normalize,
                     const float& normalSmooth)
{
    radius=float(fabs(ReconstructionFunction.polys[0].start));
    if(normalSmooth>0) postNormalSmooth=normalSmooth;
    fData.set(maxDepth,ReconstructionFunction,normalize,USE_DOT_RATIOS);
}

template<int Degree>
int Octree<Degree>::setTree(char* fileName,
                            const int& maxDepth,
                            const int& binary,
                            const int& kernelDepth,
                            const float& samplesPerNode,
                            const float& scaleFactor,
                            Point3D<float>& center,
                            float& scale,
                            const int& resetSamples)
{
    FILE* fp;
    Point3D<float> min,max,position,normal,myCenter;
    float myWidth;
    int i,cnt=0;
    float c[2*DIMENSION];
    OctNode* temp;
    int splatDepth=0;

    NodeData::UseIndex=1;
    neighborKey.set(maxDepth);
    splatDepth=std::max(splatDepth,kernelDepth);
    if(binary) fp=fopen(fileName,"rb");
    else fp=fopen(fileName,"r");
    if(!fp) return 0;


    /**     Get the min and max in each dimension   */
    while(1){
        if(binary){
            if(fread(c,sizeof(float),2*DIMENSION,fp)!=2*DIMENSION)
                break;
        } else{
            if(fscanf(fp," %f %f %f %f %f %f ",&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])!=2*DIMENSION)
                break;
        }
        for(i=0;i<DIMENSION;++i){
            if(!cnt || c[i]<min.coords[i]) min.coords[i]=c[i];
            if(!cnt || c[i]>max.coords[i]) max.coords[i]=c[i];
        }
        ++cnt;
    }

    /**     center and scale used for normalization is ready     */
    for(i=0;i<DIMENSION;++i){
        if(!i || scale<max.coords[i]-min.coords[i]) scale=float(max.coords[i]-min.coords[i]);
        center.coords[i]=float(max.coords[i]+min.coords[i])/2;
    }
    scale*=scaleFactor;
    for(i=0;i<DIMENSION;++i)
        center.coords[i]-=scale/2;


    /**     Initialize weight contribution of the node in this tree     */
    if(splatDepth>0){
        /**     set position indicator  */
        fseek(fp,0,SEEK_SET);
        cnt=0;
        while(1){
            if(binary){
                if(fread(c,sizeof(float),2*DIMENSION,fp)!=2*DIMENSION)
                    break;
            }else{
                if(fscanf(fp," %f %f %f %f %f %f ",&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])!=2*DIMENSION)
                    break;
            }
            /**     project point to [0,1]^3    */
            for(i=0;i<DIMENSION;++i)
                position.coords[i]=(c[i]-center.coords[i])/scale;

            /**     center at (0.5, 0.5, 0.5)   */
            myCenter.coords[0]=myCenter.coords[1]=myCenter.coords[2]=float(0.5);
            myWidth=float(1.0);

            /**     kick out the point which coords is smaller than 0 or bigger than 1  */
            for(i=0;i<DIMENSION;++i)
                if(position.coords[i]<myCenter.coords[i]-myWidth/2 || position.coords[i]>myCenter.coords[i]+myWidth/2)
                    break;
            if(i!=DIMENSION) continue;

            /**     Now, the point is fine to use.
             *      Update the $position contribution to every neighbor node with depth < splatDepth.
             *      myCenter and myWidth use to track the trace to its subnode. */
            temp=&tree;
            int d=0;
            while(d<splatDepth){
                NonLinearUpdateWeightContribution(temp,position);
                if(!temp->children) temp->initChildren();
                /**     cIndex = cornerIndex close to position  */
                int cIndex=OctNode::CornerIndex(myCenter,position);
                /**     step to the nearest subnode,
                 *      adjust the myCenter and myWidth     */
                temp=&temp->children[cIndex];
                myWidth/=2;
                for(i=0;i<DIMENSION;++i) {
                    if (cIndex & (1<<i)) myCenter.coords[i] += myWidth / 2;
                    else                 myCenter.coords[i] -= myWidth / 2;
                }
                ++d;
            }
            NonLinearUpdateWeightContribution(temp,position);
            ++cnt;
        }
    }

    normals=new std::vector<Point3D<float> >();
    fseek(fp,0,SEEK_SET);
    cnt=0;

    /**     Update normals in each related node */
    while(1){
        if(binary){
            if(fread(c,sizeof(float),2*DIMENSION,fp)!=2*DIMENSION)
                break;
        }else{
            if(fscanf(fp," %f %f %f %f %f %f ",&c[0],&c[1],&c[2],&c[3],&c[4],&c[5])!=2*DIMENSION)
                break;
        }
        for(i=0;i<DIMENSION;++i){
            position.coords[i]=(c[i]-center.coords[i])/scale;
            normal.coords[i]=c[DIMENSION+i];
        }
        myCenter.coords[0]=myCenter.coords[1]=myCenter.coords[2]=float(0.5);
        myWidth=float(1.0);

        for(i=0;i<DIMENSION;++i)
            if(position.coords[i]<myCenter.coords[i]-myWidth/2 || position.coords[i]>myCenter.coords[i]+myWidth/2)
                break;
        if(i!=DIMENSION) continue;

#if FORCE_UNIT_NORMALS
        float len=float(Length(normal));
        if(len>EPSILON)
            len=1.0f/len;
        len*=(2<<maxDepth);
        for(i=0;i<DIMENSION;++i){
            normal.coords[i]*=len;
        }
#endif  // FORCE_UNIT_NORMALS

        /**     Now, position and normal has been normalized    */

        if(resetSamples && samplesPerNode>0 && splatDepth){
            /**     Update the normals of the node at depth defined by $samplesPerNode
             *      and nodes' weight contribution initialized  */
            NonLinearSplatOrientedPoint(position,normal,splatDepth,samplesPerNode,1,maxDepth);
        }else{

            /**     calculate the weight $alpha */
            float alpha=1;
            temp=&tree;
            if(splatDepth){
                int d=0;
                /**     reach the splatDepth node   */
                while(d<splatDepth){
                    int cIndex=OctNode::CornerIndex(myCenter,position);
                    temp=&temp->children[cIndex];
                    myWidth/=2;
                    for(i=0;i<DIMENSION;++i){
                        if(cIndex & (1<<DIMENSION)) myCenter.coords[i]+=myWidth/2;
                        else                        myCenter.coords[i]-=myWidth/2;
                    }
                    d++;
                }
                alpha=NonLinearGetSampleWeight(temp,position);
            }

            /**     scale normal    */
            for(i=0;i<DIMENSION;++i)
                normal.coords[i]*=alpha;

            /**     Update the normals of the node only at maxDepth  */
            int d=0;
            while(d<maxDepth){
                if(!temp->children) temp->initChildren();
                int cIndex=OctNode::CornerIndex(myCenter,position);
                temp=&temp->children[cIndex];
                myWidth/=2;
                for(i=0;i<DIMENSION;++i){
                    if(cIndex & (1<<DIMENSION)) myCenter.coords[i]+=myWidth/2;
                    else                        myCenter.coords[i]-=myWidth/2;
                }
                d++;
            }
            NonLinearSplatOrientedPoint(temp,position,normal);
        }
    }
    fclose(fp);
    return cnt;
}

template<int Degree>
void Octree<Degree>::ClipTree(void) {
    OctNode* temp;
    temp=tree.nextNode();
    while(temp){
        if(temp->children){
            int hasNormals=0;
            for(int i=0;i<Cube::CORNERS && !hasNormals;++i){
                hasNormals=HasNormals(&temp->children[i],EPSILON);
            }
            if(!hasNormals)
                temp->children=NULL;
        }
        temp=tree.nextNode(temp);
    }
}

template<int Degree>
void Octree<Degree>::finalize1(const int& refineNeighbors) {
    OctNode* temp;

    if(refineNeighbors >= 0){
        RefineFunction rf;
        temp=tree.nextNode();
        while(temp){
            const int& idx=temp->nodeData.nodeIndex;
            const Point3D<float>& normal=(*normals)[idx];
            /**     if current node is valid    */
            if(idx >= 0 && Length(normal) > EPSILON){
                rf.depth=temp->depth()-refineNeighbors;
                /**     node1:      temp
                 *      radius1:    2*radius
                 *      node2:      &tree
                 *      radius2:    0.5
                 *      depth:      rf.depth
                 *      F:          rf           */
                OctNode::ProcessMaxDepthNodeAdjacentNodes(temp,2*radius,&tree,float(0.5),temp->depth()-refineNeighbors,&rf);
            }
            temp=tree.nextNode(temp);
        }
//        printf("!%d\n",rf.cnt);
    }
}

template<int Degree>
void Octree<Degree>::SetLaplacianWeights(void) {
    OctNode* temp;

    fData.setDotTables(fData.DOT_FLAG | fData.D_DOT_FLAG);
    DivergenceFunction df;

    df.ot=this;
    temp=tree.nextNode();

    /**     Since that normals is added to nodes with splatDepth,
     *      use them to update the divergence of all nodes with different depths in nodeData.value  */
    while(temp){
        const int& idx=temp->nodeData.nodeIndex;
        const Point3D<float>& normal=(*normals)[idx];
        /**     Check if $temp has valid normal     */
        if(idx < 0 || Length(normal) <= EPSILON) {
            temp=tree.nextNode(temp);
            continue;
        }
        df.normal=normal;
        /**     Set the first dimension in fData.dotTable and fData.valueTable.
         *      temp->off denote the index of corresponding baseFunctions to $temp node.    */
        for(int i=0;i<DIMENSION;++i) {
            df.index[i] = int(temp->off[i]) * fData.res;
        }

        /**     radius:     1.5     */

        /**     $df member index has been set to be the index of $temp node's baseFunction.
         *      The sub-node of &tree which is near to &temp.
         *      Their nodeData.value was added with the divergence of &temp in them.
         *      normals of the nodes is used to calculate this value.               */
        OctNode::ProcessNodeAdjacentNodes(temp,radius,&tree,radius,&df);
        temp=tree.nextNode(temp);
    }
    fData.clearDotTables(fData.D_DOT_FLAG);
    temp=tree.nextNode();

    /**     Save the normal length as weight of node
     *      in temp->nodeData.centerWeightContribution     */
    while(temp){
        const int& idx=temp->nodeData.nodeIndex;
        if(idx<0)
            temp->nodeData.centerWeightContribution=0;
        else {
            const Point3D<float>& normal=(*normals)[idx];
            temp->nodeData.centerWeightContribution = float(Length(normal));
        }
        temp=tree.nextNode(temp);
    }

    /**     normals has been transformed to weight    */
    delete normals;
    normals=NULL;
}

template<int Degree>
void Octree<Degree>::finalize2(const int& refineNeighbors){
    OctNode* temp;
    if(refineNeighbors>=0){
        RefineFunction rf;
        temp=tree.nextNode();
        while(temp){
            const float& divergence=temp->nodeData.value;
            if(fabs(divergence)>EPSILON){
                rf.depth=temp->depth()-refineNeighbors;
                OctNode::ProcessMaxDepthNodeAdjacentNodes(temp,2*radius,&tree,float(0.5),rf.depth,&rf);
            }
            temp=tree.nextNode(temp);
        }
    }
}

template<int Degree>
int Octree<Degree>::LaplacianMatrixIteration(const int& subdivideDepth){
    int i,iter=0;
    SortedTreeNodes sNodes;
    double t;
    fData.setDotTables(fData.D2_DOT_FLAG);
    /**     nodeData.nodeIndex of all nodes in &tree will be refresh.
     *      They will be the index in SortedTreeNodes       */
    sNodes.set(tree,true);

    SparseMatrix<float>::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);

    /**     eliminate the &tree node divergence     */
    sNodes.treeNodes[0]->nodeData.value=0;
    for(i=1;i<sNodes.maxDepth;++i){
        if(subdivideDepth > 0)
            iter+=SolveFixedDepthMatrix(i,subdivideDepth,sNodes);
        else
            iter+=SolveFixedDepthMatrix(i,sNodes);
    }
    SparseMatrix<float>::Allocator.reset();
    fData.clearDotTables(fData.DOT_FLAG | fData.D_DOT_FLAG | fData.D2_DOT_FLAG);
    return iter;
}

template<int Degree>
float Octree<Degree>::GetIsoValue(void){
    const OctNode* temp;
    float isoValue,weightSum,w;
    float myRadius;
    PointIndexValueFunction cf;
    Point3D<float> center;
    float width;

    fData.setValueTables(fData.VALUE_FLAG,0);
    cf.valueTables=fData.valueTables;
    cf.res2=fData.res2;
    myRadius=radius;
    isoValue=weightSum=0;
    temp=tree.nextNode();
    while(temp){
        /**     length of normal, only maxDepth node with normal has WeightContribution     */
        w=temp->nodeData.centerWeightContribution;
        if(w>EPSILON){
            cf.value=0;
            VertexData::CenterIndex(temp,fData.depth,cf.index);
            temp->centerAndWidth(center,width);
            /**     For those nodes that is close to center enough,
             *      sum up (x * Fo) = (x * Fo,x * Fo,y * Fo,z)          */
            OctNode::ProcessPointAdjacentNodes(center,&tree,myRadius,&cf);
//            printf("%f %f\n",cf.value,w);
            isoValue+=cf.value*w;
            weightSum+=w;
        }
        temp=tree.nextNode(temp);
    }
//    printf("cnt:%d\n",cf.cnt);
    return isoValue/weightSum;
}

//template<int Degree>
//void Octree<Degree>::GetMCIsoTriangles(const float& isoValue,CoredMeshData* mesh,const int& fullDepthIso){
//    double t;
//    OctNode* temp;
//    SortedTreeLeaves sLeaves;
//
//    unordered_map<long long,int> roots;
//    unordered_map<long long,Point3D<float> > *normalHash=new unordered_map<long long,Point3D<float> >();
//
//    SetIsoSurfaceCorners(isoValue,0,fullDepthIso);
//
//}

