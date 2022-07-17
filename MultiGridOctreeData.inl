
#define FORCE_UNIT_NORMALS 1
#define USE_DOT_RATIOS 1

#include "Octree.h"

// VertexData
long long VertexData::CenterIndex(const OctNode *node, const int &maxDepth, int idx[DIMENSION]) {
    int d,o[3];
    node->depthAndOffset(d,o);
    for(int i=0;i<DIMENSION;++i){
        idx[i]=BinaryNode<float>::CornerIndex(maxDepth+1,d+1,o[i]<<1,1);
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

void SortedTreeNodes::set(OctNode& root, const int& setIndex) {
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

template<int Degree>
void Octree<Degree>::DivergenceFunction::Function(OctNode *node1, OctNode *node2) {
    for(int i=0;i<3;++i)
        scratch[i]=index[i]+int(node1->off[i]);
    node1->nodeData.value+=ot->GetDivergence(scratch,normal);
}

template<int Degree>
void Octree<Degree>::LaplacianProjectionFunction::Function(OctNode *node1, OctNode *node2) {
    for(int i=0;i<3;++i)
        scratch[i]=index[i]+int(node1->off[i]);
    node1->nodeData.value-=float(ot->GetLaplacian(scratch)*value);
}

/*
template<int Degree>
int Octree<Degree>::LaplacianMatrixFunction::Function(OctNode *node1, OctNode *node2) {
    float temp;
    int d1=int(node1->d);
    int x1,y1,z1;
    x1=int(node1->off[0]);
    y1=int(node1->off[1]);
    z1=int(node1->off[2]);

    int dDepth=d2-d1;
    int dcoords;
    dcoords=(x2>>dDepth)-x1;

}
*/

template<int Degree>
bool Octree<Degree>::NonLinearUpdateWeightContribution(OctNode* node, const Point3D<float>& position) {
    int i,j,k;
    /**     Get the neighbor of $node       */
    OctNode::Neighbors& neighbors=neighborKey.setNeighbors(node);

    double x,dx[DIMENSION][3];
    double width;
    /**     Get the center and width of $node   */
    Point3D<float> center;
    float w;
    node->centerAndWidth(center,w);

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
