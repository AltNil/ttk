#include <Debug.h>
#include <AbstractTriangulation.h>
#include <iostream>
#include <vector>

namespace ttk {
  
  class implicitQuadrangulation final : public ttk::AbstractTriangulation {

  public:
    implicitQuadrangulation(std::array<int,3> extent_) : extent(extent_){
      dx=extent[0];
      dxm1=extent[0]-1;
      dy=extent[1];
      dym1=extent[1]-1;
      dz=extent[2];
      dzm1=extent[2]-1;
    };

    int preconditionVertexNeighbors() final {
      return 1;
    };
    ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId& vertexId)
      const final {
        return 0; // TODO machen
    };
    int getVertexNeighbor(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId
    ) const final {
      return 1; // TODO machen
    };


    ttk::SimplexId getVertexEdgeNumber(const ttk::SimplexId& vertexId) const final {
		  return 0;
    };

    int getVertexEdge(
      const SimplexId &vertexId,
      const int &localEdgeId,
      SimplexId &edgeId
    ) const final {
      return 1;
    };

    int preconditionVertexEdges() final {
      return 1;
    };

    /**
     * Returns Quadrants instead of Triangles
     */
    ttk::SimplexId getVertexTriangleNumber(const ttk::SimplexId& vertexId) const final {
      return 0;
    };

    int getVertexTriangle(
      const SimplexId &vertexId,
      const int &localTriangleId,
      SimplexId &triangleId
    ) const final {
      return 1;
    };

    int preconditionVertexTriangles() final {
      return 1;
    };

    ttk::SimplexId getVertexStarNumber(const ttk::SimplexId& vertexId) const final {
      return 0;
    };

    int getVertexStar(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId
    ) const final {
      return 1;
    };

    int preconditionVertexStars() final {
      return 1;
    };





  private:
    inline int getCaseID(const ttk::SimplexId& vertexId) const{
      // Calculate the coordinates of the given vertexId
      auto xyDim = dx*dy;
      int z = vertexId/xyDim;
      int xy = (vertexId-z*xyDim);
      int y = xy/dx;
      int x = xy-y*dx;

      // case Clasification by coordinates
      int caseID = 0;
      
      caseID += (!bool(x)); // +1 if x == 0
      caseID += (!bool(x-dxm1))*2; // +2 if x == extent[0]-1
      caseID += (!bool(y))*3; // +3 if y == 0
      caseID += (!bool(y-dym1))*6; // +6 if y == extent[1]-1
      caseID += (!bool(z))*9; // +9 if z == 0
      caseID += (!bool(z-dzm1))*18; // +18 if z == extent[2]-1
      if(caseID>=27) return 26; // can not be the case unless Freudenthal3D is used in 2D scenario
      return caseID;
    };

    const std::array<int,3> extent;
    int dx;
    int dxm1;
    int dy;
    int dym1;
    int dz;
    int dzm1;

    
  };
} // namespace ttk
