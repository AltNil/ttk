#include <Debug.h>
#include <AbstractTriangulation.h>
#include <iostream>
#include <vector>

namespace ttk {
  
  class ImplicitQuadrangulation final : public ttk::AbstractTriangulation {

  public:
    ImplicitQuadrangulation(std::array<int,3> extent_) : extent(extent_){
        dx=extent[0];
        dxm1=extent[0]-1;
        dy=extent[1];
        dym1=extent[1]-1;
        dz=extent[2];
        dzm1=extent[2]-1;



        for (std::vector<std::array<int,3>> caseId : lutNeighborOffset){
            std::vector<int> shifts;
            for (std::array<int,3> coordDeltas: caseId){
                shifts.push_back(coordDeltas[0]+coordDeltas[1]*dx+(dx*dy*coordDeltas[2]));
            }
            lutNeighborIndexOffset.push_back(shifts);
        }

        int offset = getNumberOfVertices();
        for (std::vector<std::array<int,3>> caseId : lutNeighborOffset){
            std::vector<int> shifts;
            for (std::array<int,3> coordDeltas: caseId){
                int off = 0;
                if (coordDeltas[1] != 0){
                    off += offset;
                }
                if (coordDeltas[2] != 0){
                    off += offset;
                }
                int indexDelta = coordDeltas[0];
                if (coordDeltas[0] == -1){
                    indexDelta += coordDeltas[0];
                }
                if (coordDeltas[1] == -1){
                    indexDelta += coordDeltas[1]*dx;
                }
                if (coordDeltas[2] == -1){
                    indexDelta += dx*dy*coordDeltas[2];
                }
                
                shifts.push_back(indexDelta+off);
            }
            lutEdgeIndexOffset.push_back(shifts);
        }

        //Quadrats

        for (std::vector<std::array<int,3>> caseId : lutQuadradOffset){
            std::vector<int> shifts;
            for (std::array<int,3> coordDeltas: caseId){
                //xy direction -> + 0
                int off = 0;
                //xz direction
                if (coordDeltas[0] != 0 && coordDeltas[2] != 0){
                    off += offset;
                }
                //yz direction
                if (coordDeltas[1] != 0 && coordDeltas[2] != 0){
                    off += offset;
                }
                int indexDelta = coordDeltas[0];
                if (coordDeltas[0] == -1){
                    indexDelta += coordDeltas[0];
                }
                if (coordDeltas[1] == -1){
                    indexDelta += coordDeltas[1]*dx;
                }
                if (coordDeltas[2] == -1){
                    indexDelta += dx*dy*coordDeltas[2];
                }
                
                shifts.push_back(indexDelta+off);
            }
            lutQuadratIndexOffset.push_back(shifts);
        }

        //cube
        for (std::vector<std::array<int,3>> caseId : lutCubeOffest){
            std::vector<int> shifts;
            for (std::array<int,3> coordDeltas: caseId){
                int indexDelta = 0;
                if (coordDeltas[0] == -1){
                    indexDelta += coordDeltas[0];
                }
                if (coordDeltas[1] == -1){
                    indexDelta += coordDeltas[1]*dx;
                }
                if (coordDeltas[2] == -1){
                    indexDelta += dx*dy*coordDeltas[2];
                }
                shifts.push_back(indexDelta);
            }
            lutCubeIndexOffset.push_back(shifts);
        }
    };


    ttk::SimplexId getNumberOfVertices() const final {
      return this->extent[0]*this->extent[1]*this->extent[2];
    };


    int preconditionVertexNeighbors() final {
      return 1;
    };
    ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId& vertexId)
      const final {
        return lutVertexNeighborNumber[getCaseID(vertexId)];
    };
    int getVertexNeighbor(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId
    ) const final {
      neighborId = SimplexId(vertexId + lutNeighborIndexOffset[getCaseID(vertexId)][localNeighborId]);
      return 1;
    };


    ttk::SimplexId getVertexEdgeNumber(const ttk::SimplexId& vertexId) const final {
        return lutVertexNeighborNumber[getCaseID(vertexId)];
    };

    int getVertexEdge(
      const SimplexId &vertexId,
      const int &localEdgeId,
      SimplexId &edgeId
    ) const final {
        edgeId = SimplexId(vertexId + lutEdgeIndexOffset[getCaseID(vertexId)][localEdgeId]);
        return 1;
    };

    int preconditionVertexEdges() final {
        return 1;
    };

    /**
     * Returns Quadrants instead of Triangles
     */
    ttk::SimplexId getVertexTriangleNumber(const ttk::SimplexId& vertexId) const final {
        return lutQuadratNumber[getCaseID(vertexId)];
    };

    int getVertexTriangle(
      const SimplexId &vertexId,
      const int &localTriangleId,
      SimplexId &triangleId
    ) const final {
        triangleId = SimplexId(vertexId + lutQuadratIndexOffset[getCaseID(vertexId)][localTriangleId]);
        return 1;
    };

    int preconditionVertexTriangles() final {
        return 1;
    };

    ttk::SimplexId getVertexStarNumber(const ttk::SimplexId& vertexId) const final {
        return lutCubeNumber[getCaseID(vertexId)];
    };

    int getVertexStar(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId
    ) const final {
        starId = SimplexId(vertexId + lutCubeIndexOffset[getCaseID(vertexId)][localStarId]);
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

    std::array<int, 27> lutVertexNeighborNumber = {6,5,5,5,4,4,5,4,4,5,4,4,4,3,3,4,3,3,5,4,4,4,3,3,4,3,3};
    const std::vector<std::vector<std::array<int,3>>> lutNeighborOffset = {
      {{-1,0,0},{0,0,-1},{0,-1,0},{1,0,0},{0,0,1},{0,1,0}},
      {{0,0,-1},{0,-1,0},{1,0,0},{0,0,1},{0,1,0}},
      {{0,0,1},{0,1,0},{-1,0,0},{0,0,-1},{0,-1,0}},
      {{1,0,0},{0,0,1},{0,1,0},{-1,0,0},{0,0,-1}},
      {{0,0,-1},{1,0,0},{0,0,1},{0,1,0}},
      {{0,0,1},{0,1,0},{-1,0,0},{0,0,-1}},
      {{-1,0,0},{0,0,-1},{0,-1,0},{1,0,0},{0,0,1}},
      {{0,0,-1},{0,-1,0},{1,0,0},{0,0,1}},
      {{-1,0,0},{0,0,-1},{0,-1,0},{0,0,1}},
      {{0,-1,0},{1,0,0},{0,0,1},{0,1,0},{-1,0,0}},
      {{0,-1,0},{1,0,0},{0,0,1},{0,1,0}},
      {{0,-1,0},{0,0,1},{0,1,0},{-1,0,0}},
      {{1,0,0},{0,0,1},{0,1,0},{-1,0,0}},
      {{1,0,0},{0,0,1},{0,1,0}},
      {{0,0,1},{0,1,0},{-1,0,0}},
      {{-1,0,0},{0,-1,0},{1,0,0},{0,0,1}},
      {{0,-1,0},{1,0,0},{0,0,1}},
      {{-1,0,0},{0,-1,0},{0,0,1}},
      {{0,1,0},{-1,0,0},{0,0,-1},{0,-1,0},{1,0,0}},
      {{0,1,0},{0,0,-1},{0,-1,0},{1,0,0}},
      {{0,1,0},{-1,0,0},{0,0,-1},{0,-1,0}},
      {{1,0,0},{0,1,0},{-1,0,0},{0,0,-1}},
      {{1,0,0},{0,1,0},{0,0,-1}},
      {{0,1,0},{-1,0,0},{0,0,-1}},
      {{-1,0,0},{0,0,-1},{0,-1,0},{1,0,0}},
      {{0,0,-1},{0,-1,0},{1,0,0}},
      {{-1,0,0},{0,0,-1},{0,-1,0}}
    };
    std::vector<std::vector<int>> lutNeighborIndexOffset;   

    std::vector<std::vector<int>> lutEdgeIndexOffset;

    std::array<int, 27> lutQuadratNumber = {12,6,8,5,5,5,8,5,5,8,5,5,5,3,3,5,3,3,8,6,6,5,3,3,5,3,3};
    std::vector<std::vector<std::array<int,3>>> lutQuadradOffset = {
        {{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0},{-1,0,-1},{1,0,-1},{1,0,1},{-1,0,1},{0,-1,-1},{0,-1,1},{0,1,1},{0,1,-1}},
        {{1,-1,0},{1,1,0},{1,0,-1},{1,0,1},{0,1,1},{0,1,-1}},
        {{-1,-1,0},{-1,1,0},{-1,0,-1},{-1,0,1},{0,-1,-1},{0,-1,1},{0,1,1},{0,1,-1}},
        {{1,1,0},{-1,1,0},{1,0,1},{-1,0,1},{0,1,1}},
        {{1,1,0},{1,0,-1},{1,0,1},{0,1,1},{0,1,-1}},
        {{-1,1,0},{-1,0,-1},{-1,0,1},{0,1,1},{0,1,-1}},
        {{-1,-1,0},{1,-1,0},{-1,0,-1},{1,0,-1},{1,0,1},{-1,0,1},{0,-1,-1},{0,-1,1}},
        {{1,-1,0},{1,0,-1},{1,0,1},{0,-1,-1},{0,-1,1}},
        {{-1,-1,0},{-1,0,-1},{-1,0,1},{0,-1,-1},{0,-1,1}},
        {{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0},{1,0,1},{-1,0,1},{0,-1,1},{0,1,1}},
        {{1,-1,0},{1,1,0},{1,0,1},{0,-1,1},{0,1,1}},
        {{1,1,0},{-1,1,0},{-1,0,1},{0,-1,1},{0,1,1}},
        {{-1,1,0},{1,1,0},{-1,0,1},{1,0,1},{0,1,1}},
        {{1,1,0},{1,0,1},{0,1,1}},
        {{-1,1,0},{-1,0,1},{0,1,1}},
        {{-1,-1,0},{1,-1,0},{1,0,1},{-1,0,1},{0,-1,1}},
        {{1,-1,0},{1,0,1},{0,-1,1}},
        {{-1,-1,0},{-1,0,1},{0,-1,1}},
        {{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0},{-1,0,-1},{1,0,-1},{0,-1,-1},{0,1,-1}},
        {{1,-1,0},{1,1,0},{1,0,-1},{1,0,1},{0,-1,-1},{0,1,-1}},
        {{-1,-1,0},{-1,1,0},{-1,0,-1},{-1,0,1},{0,-1,-1},{0,1,-1}},
        {{1,1,0},{-1,1,0},{-1,0,-1},{1,0,-1},{0,1,-1}},
        {{1,1,0},{1,0,-1},{0,1,-1}},
        {{-1,1,0},{-1,0,-1},{0,1,-1}},
        {{-1,-1,0},{1,-1,0},{-1,0,-1},{1,0,-1},{0,-1,-1}},
        {{1,-1,0},{1,0,-1},{0,-1,-1}},
        {{-1,-1,0},{-1,0,-1},{0,-1,-1}}
    };
    std::vector<std::vector<int>> lutQuadratIndexOffset;
  
    std::array<int,27> lutCubeNumber = {8,4,4,4,2,2,4,2,2,4,2,2,2,1,1,2,1,1,4,2,2,2,1,1,2,1,1};
    std::vector<std::vector<std::array<int,3>>> lutCubeOffest = {
        {{-1,-1,-1},{1,-1,-1},{1,-1,1},{-1,-1,1},{-1,1,-1},{1,1,-1},{1,1,1},{-1,1,1}},
        {{1,-1,-1},{1,-1,1},{1,1,1},{1,1,-1}},
        {{-1,-1,-1},{-1,-1,1},{-1,1,1},{-1,1,-1}},
        {{-1,1,-1},{1,1,-1},{1,1,1},{-1,1,1}},
        {{1,1,-1},{1,1,1}},
        {{-1,1,-1},{-1,1,1}},
        {{-1,-1,-1},{1,-1,-1},{1,-1,1},{-1,-1,1}},
        {{1,-1,-1},{1,-1,1}},
        {{-1,-1,-1},{-1,-1,1}},
        {{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}},
        {{1,-1,1},{1,1,1}},
        {{-1,-1,1},{-1,1,1}},
        {{1,1,1},{-1,1,1}},
        {{1,1,1}},
        {{-1,1,1}},
        {{-1,-1,1},{1,-1,1}},
        {{1,-1,1}},
        {{-1,-1,1}},
        {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1}},
        {{1,-1,-1},{1,1,-1}},
        {{-1,-1,-1},{-1,1,-1}},
        {{-1,1,-1},{1,1,-1}},
        {{1,1,-1}},
        {{-1,1,-1}},
        {{-1,-1,-1},{1,-1,-1}},
        {{1,-1,-1}},
        {{-1,-1,-1}}
    };
    std::vector<std::vector<int>> lutCubeIndexOffset;

  }; // class ImplicitQuadrangulation
} // namespace ttk
