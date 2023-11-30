#include <Debug.h>
#include <AbstractTriangulation.h>
#include <iostream>
#include <vector>

namespace ttk {

class Freudenthal3D final : public ttk::AbstractTriangulation {

  public:

    Freudenthal3D(std::array<int,3> extent_) : extent(extent_){

      // precalculate the implicite index offsets for neighbors
      std::vector<std::vector<int>> deltas;
      for(std::vector<std::array<int,3>> caseId : lutNeghborOffset){
        std::vector<int> shifts;
        for(std::array<int,3> coordDeltas: caseId){
          shifts.push_back(coordDeltas[0]+coordDeltas[1]*extent[0]+(extent[0]*extent[1]*coordDeltas[2]));
        }
        lutIndexOffset.push_back(shifts);
      }

      // precalculation for the edge calculations
      int ex = extent[0];
      int ey = extent[1];
      int ez = extent[2];
      d1 = (ex-1)*ey*ez;
      d2 = d1 + ex*(ey-1)*ez;
      d3 = d2 + ex*ey*(ez-1);
      d4 = d3 + (ex-1)*(ey-1)*ez;
      d5 = d4 + ex*(ey-1)*(ez-1);
      d6 = d5 + (ex-1)*ey*(ez-1);
    }

    ttk::SimplexId getNumberOfVertices() const final {
      return this->extent[0]*this->extent[1]*this->extent[2];
    };


    /**
     * Returns the number of neighbors for a specified vertex
     * \param vertexId The Id of the vertex for which the number of neighbors should be searched
     * \return Number vor neigbors the given vertex has
    */
    ttk::SimplexId getVertexNeighborNumber(const ttk::SimplexId& vertexId)
     const final {
      // return the neigbor by the caseID of the given vertex
      return lutNumNeighbour3d[getCaseID(vertexId)];
    };

    /**
     * get a specified neighbor for a given vertex
     * \param vertexId The id of the vertex from which the neighbbor should be searched for
     * \param localNeighborId The local nubmer of the neighbor to be returned (Should be between 0 and vertexNeighborNumber-1)
     * \param [out] neighborId The output in which the explicit neighborId should be stored in
     * \return 1 if successfull
     */
    inline int getVertexNeighbor(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId
    ) const final {
      // vertexid + offset of the local neighbor dependant of the case the vertex is in
      neighborId = SimplexId(vertexId+lutIndexOffset[getCaseID(vertexId)][localNeighborId]);
      return 1;
    };

    int preconditionVertexNeighbors() final {
      return 1;
    };

    /**
     * Get the number of adjacent vertices from a given vertex
     * \param vertexId The id of the vertex the number of edges should be looked up for
     * \return Number of adjacent edges
     */
    ttk::SimplexId getVertexEdgeNumber(const SimplexId& vertexId) const {
      // As every Neighbor vertex induces one edge with the given vertex, it is sufficient to return the number of neighbors.
      return lutNumNeighbour3d[getCaseID(vertexId)];
    };

    /**
     * Get the id of an specific edge where one endpoint is the given vertex
     * \param vertexId The id of the vertex where the edge should start
     * \param localEdgeId The local number of the edge to be searched for. Must be between 0 and vertexEdgeNumber -1.
     * \param [out] edgeId The output in which the explicit edge id should be stored in
     * \return 1 if successfull, -1 else
     */
    int getVertexEdge ( 	const SimplexId &  	vertexId,
      const int &  	localEdgeId,
      SimplexId &  	edgeId 
    ) 		const {
      
      // prefetch dimensions as those are needed multiple times
      int ex = extent[0];
      int ey = extent[1];

      // Calculate the coordinates of the given vertexId
      auto xyDim = ex*ey;
      int z = vertexId/xyDim;
      int xy = (vertexId-z*xyDim);
      int y = xy/ex;
      int x = xy-y*ex;
      
      //search for the direction of the edge and tarnsform it to an id
      auto coordOffsets = lutNeghborOffset[getCaseID(vertexId)][localEdgeId];
      int lutIndex = coordOffsets[0]+1+3*(coordOffsets[1]+1)+9*(coordOffsets[2]+1); // generate unique id from coords

      // calculate the explicit id of the searched edge with starting coordinates x, y and z from the starting point(given vertex)
      // the pdf cases refer to the different layer, in which the edges are enumerated through the grid.
      switch (lutIndex)
      {
      case 14: // pdf case 1 positive
        edgeId = x + y*(ex-1)+(z*(ex-1)*ey);
        break;
      case 16: // pdf case 2 positive
        edgeId = d1 + x + (y*ex) + (z*ex*(ey-1));
        break;
      case 22: // pdf case 3 positive
        edgeId = d2 + x + (y*ex) + (z*ex*ey);
        break;
      case 15: // pdf case 4 positive
        edgeId = d3 + (x-1) + y*(ex-1) + z * (ex-1)*(ey-1);
        break;
      case 25: // pdf case 5 positive
        edgeId = d4 + x + y*ex + z * ex * (ey-1);
        break;
      case 21: //pdf case 6 positive
        edgeId = d5 + (x-1) + y * (ex-1) + z * (ex-1) * ey;
        break;
      case 24: // pdf case 7 positive
        edgeId = d6 + (x-1) + y * (ex-1) + z * (ex-1) * (ey-1);
        break;
      case 12: // pdf case 1 negative
        edgeId = (x-1)+y*(ex-1)+(z*(ex-1)*ey);
        break;
      case 10: // pdf case 2 negative
        edgeId = d1 + x + ((y-1)*ex) + (z*ex*(ey-1));
        break;
      case 4: // pdf case 3 negative
        edgeId = d2 + x + (y*ex) + ((z-1)*ex*ey);
        break;
      case 11: // pdf case 4 negative
        edgeId = d3 + (x) + (y-1)*(ex-1) + z * (ex-1)*(ey-1);
        break;
      case 1: // pdf case 5 negative
        edgeId = d4 + x + (y-1) * ex + (z-1) * ex * (ey-1);
        break;
      case 5: //pdf case 6 negative
        edgeId = d5 + (x) + y * (ex-1) + (z-1) * (ex-1) * ey;  
        break;
      case 2: // pdf case 7 negative
        edgeId = d6 + (x) + (y-1) * (ex-1) + (z-1) * (ex-1) * (ey-1);
        break;
      
      default:
        return -1; //something went wrong
        break;
      }

      return 1;
    };


    int preconditionVertexTriangles() final {
      return 1;
    };
    ttk::SimplexId getVertexTriangleNumber (const SimplexId &vertexId) const {
      return -1;
    };
    int getVertexTriangle(
      const SimplexId &vertexId,
      const int &localTriangleId,
      SimplexId &triangleId
    ) const {
      return 1;
    };


    int preconditionVertexStars() final {
      return 1;
    }
    ttk::SimplexId getVertexStarNumber (const SimplexId &vertexId) const {
      return -1;
    };
    int getVertexStar(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId
    ) const {
      return 1;
    };





  private:
    /**
     * Transpose the ipmlicite index of an vertex into a case identifyer which denotes
     * on which border the vertex is positioned.
     * \param vertexId The Id of the vertex for which the caseId should be searched
     * \return CaseId as an integer between 0 and 26
     */
    inline int getCaseID(const ttk::SimplexId& vertexId) const{
      int ex = extent[0];
      int ey = extent[1];
      int ez = extent[2];

      // Calculate the coordinates of the given vertexId
      auto xyDim = ex*ey;
      int z = vertexId/xyDim;
      int xy = (vertexId-z*xyDim);
      int y = xy/ex;
      int x = xy-y*ex;

      // case Clasification by coordinates
      int caseID = 0;
      
      caseID += (!bool(x)); // +1 if x == 0
      caseID += (!bool(x-(ex-1)))*2; // +2 if x == extent[0]-1
      caseID += (!bool(y))*3; // +3 if y == 0
      caseID += (!bool(y-(ey-1)))*6; // +6 if y == extent[1]-1
      caseID += (!bool(z))*9; // +9 if z == 0
      caseID += (!bool(z-(ez-1)))*18; // +18 if z == extent[2]-1
      if(caseID>=27) return 26; // can not be the case unless Freudenthal3D is used in 2D scenario
      return caseID;
    };

    /*
    int getVertexId(const std::array<int,3>& coordinates)  const{
      auto result =  coordinates[0]+coordinates[1]*extent[0]+(extent[0]*extent[1]*coordinates[2]);
      return result;
    };*/


    const std::array<int,3> extent;
    const int lutNumNeighbour3d[27] = {14,10,10,10,6,8,10,8,6,10,6,8,8,4,7,6,4,4,10,8,6,6,4,4,8,7,4};
    const std::vector<std::vector<std::array<int,3>>> lutNeghborOffset = {
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 1, 1}, {0, 1, 1}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {0, 0, -1}, {0, -1, -1}, {0, -1, 0}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {1, 0, 0}, {1, 0, -1}, {0, 0, -1}},
      {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 1}, {0, 0, -1}, {1, 0, -1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {0, 0, -1}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {-1, 0, 0}, {-1, 0, 1}, {0, 0, 1}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 0, 1}},
      {{0, -1, 0}, {-1, 0, 0}, {-1, 0, 1}, {0, 0, 1}, {0, -1, -1}, {0, 0, -1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}},
      {{0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {0, -1, 0}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}, {1, 0, 0}},
      {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 1}},
      {{-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {-1, 0, 1}, {0, 0, 1}, {-1, 1, 1}, {0, 1, 1}},
      {{0, -1, 0}, {-1, 0, 0}, {-1, 0, 1}, {0, 0, 1}, {1, -1, 0}, {1, 0, 0}},
      {{0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 0, 1}},
      {{0, -1, 0}, {-1, 0, 0}, {-1, 0, 1}, {0, 0, 1}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {0, 1, 0}},
      {{0, 0, -1}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {0, -1, -1}, {0, -1, 0}},
      {{0, 0, -1}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 0, -1}, {1, 0, 0}},
      {{0, 0, -1}, {1, 0, -1}, {1, 0, 0}, {0, 1, 0}},
      {{0, 0, -1}, {-1, 0, 0}, {-1, 1, 0}, {0, 1, 0}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}, {-1, 0, 0}},
      {{0, -1, -1}, {1, -1, -1}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {1, -1, 0}, {1, 0, 0}},
      {{0, -1, -1}, {0, 0, -1}, {0, -1, 0}, {-1, 0, 0}}
    };
    std::vector<std::vector<int>> lutIndexOffset;
    int d1;
    int d2;
    int d3;
    int d4;
    int d5;
    int d6;
  };
}