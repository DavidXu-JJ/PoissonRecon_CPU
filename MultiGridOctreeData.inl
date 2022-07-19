
#include "Octree.h"

#define FORCE_UNIT_NORMALS 1
#define USE_DOT_RATIOS 1

#define MEMORY_ALLOCATOR_BLOCK_SIZE 1<<12
#define ITERATION_POWER 1.0/3

const float EPSILON=float(1e-6);
const float ROUND_EPS=float(1e-5);

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
    for(i=0;i<maxDepth;++i){nodeCount[i]=0;}
    for(i=0;i<cnt;++i){
        if(setIndex) treeNodes[i]->nodeData.nodeIndex=i;
        nodeCount[treeNodes[i]->depth()+1]++;
    }
    for(i=1;i<=maxDepth;++i)
        nodeCount[i]+=nodeCount[i-1];

}

// Octree
template<int Degree>
float Octree<Degree>::GetLaplacian(const int idx[DIMENSION]) const
{
#if USE_DOT_RATIOS
    return Real(fData.dotTable[idx[0]]*fData.dotTable[idx[1]]*fData.dotTable[idx[2]]*(fData.d2DotTable[idx[0]]+fData.d2DotTable[idx[1]]+fData.d2DotTable[idx[2]]));
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
    return Real(dot*(fData.dDotTable[idx[0]]*normal.coords[0]+fData.dDotTable[idx[1]]*normal.coords[1]+fData.dDotTable[idx[2]]*normal.coords[2]));
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
     *      $scratch has been set to be the index of <node2, node1>.
     *      $node1 is changing.                                      */
    node1->nodeData.value+=ot->GetDivergence(scratch,normal);
}

// LaplacianProjectionFunction
template<int Degree>
void Octree<Degree>::LaplacianProjectionFunction::Function(OctNode *node1, OctNode *node2) {
    for(int i=0;i<3;++i)
        scratch[i]=index[i]+int(node1->off[i]);
    /**     $scratch contains node1 and node2 index information     */
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

// RefineFunction
template<int Degree>
void Octree<Degree>::RefineFunction::Function(OctNode* node1,const OctNode* node2){
    if(!node1->children && node1->depth()<this->depth)
        node1->initChildren();
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
        OctNode const * const temp=sNodes.treeNodes[i];
        mf.elementCount=0;
        /**     mf.d2 is supposed to be $depth,
         *      every one is at the same depth  */
        mf.d2=int(temp->d);

        mf.x2=int(temp->off[0]);
        mf.y2=int(temp->off[1]);
        mf.z2=int(temp->off[2]);

        /**     input $temp node index information to $mf   */
        mf.index[0]=mf.x2*fData.res;
        mf.index[1]=mf.y2*fData.res;
        mf.index[2]=mf.z2*fData.res;

        /**     radius1:    1.500001
         *      radius2:    0.5         */
        OctNode::ProcessTerminatingNodeAdjacentNodes(temp,myRadius-float(0.5),
                                                     &tree,float(0.5),
                                                     &mf);
        /**     Set row [0, number of nodes with $depth] to be each round $elementCount */
        matrix.SetRowSize(i-sNodes.nodeCount[depth],mf.elementCount);
        memcpy(matrix.m_ppElements[i-sNodes.nodeCount[depth]],
               mf.rowElements, sizeof(MatrixEntry<float>) * mf.elementCount);
    }
    free(mf.rowElements);
}

//template<int Degree>
//int Octree<Degree>::SolveFixedDepthMatrix(const int& depth, const SortedTreeNodes& sNodes){
//    int i,iter=0;
//    Vector<double> V,Solution;
//    /**     generating lower triangular matrix  */
//    SparseSymmetricMatrix<float> matrix;
//    float myRadius,myRadius1,myRadius2;
//    float dx,dy,dz;
//    int x1,x2,y1,y2,z1,z2;
//
//    /**     current process depth is fixed,
//     *      extract all nodes with this depth   */
//    V.Resize(sNodes.nodeCount[depth+1]-sNodes.nodeCount[depth]);
//    for(i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];++i)
//        V[i-sNodes.nodeCount[depth]]=sNodes.treeNodes[i]->nodeData.value;
//    /**     empty the allocator     */
//    SparseSymmetricMatrix<float>::Allocator.rollBack();
//    GetFixedDepthLaplacian(matrix,depth,sNodes);
//    iter+=SparseSymmetricMatrix<float>::Solve(matrix,V,
//                                              int(pow(matrix.rows,ITERATION_POWER)),
//                                              Solution,double(EPSILON),1);
//
//    for(i=sNodes.nodeCount[depth];i<sNodes.nodeCount[depth+1];++i){
//        sNodes.treeNodes[i]->nodeData.value = float(Solution[i-sNodes.nodeCount[depth]]);
//    }
//
//}

//template<int Degree>
//int Octree<Degree>::SolveFixedDepthMatrix(const int& depth, const int& startingDepth, const SortedTreeNodes& sNodes){
//
//}



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
                    neighbors.neighbors[i][j][k]->nodeData.centerWeightContribution+=float(dx[0][i]*dx[1][j]*dx[2][k]);
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
            if(cIndex & (1<<DIMENSION)) myCenter.coords[i]+=myWidth/2;
            else                        myCenter.coords[i]-=myWidth/2;
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
            if(cIndex & (1<<DIMENSION)) myCenter.coords[i]+=myWidth/2;
            else                        myCenter.coords[i]-=myWidth/2;
        }
    }

    /**     scale the original normal   */
    width=1.0/(1<<temp->depth());
    for(i=0;i<DIMENSION;++i) {
        n.coords[i] = normal.coords[i] * alpha
                      / float(pow(width, 3)
                      * float(dx));
    }

    /**     update the contribution to normals[temp] at $topDepth   */
    NonLinearSplatOrientedPoint(temp,position,n);
    /**     In case that newDepth is between [topDepth-1, topDepth],
     *      update the contribution to normals[temp->parent] at $topDepth-1 */
    if(fabs(1.0-dx)>EPSILON){
        dx=float(1.0/dx);
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
            NonLinearUpdateWeightContribution(position,normal,splatDepth,samplesPerNode,1,maxDepth);
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
            if(idx>=0 && Length(normal) > EPSILON){
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
         *      The subnode of &tree which is near to &temp.
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
