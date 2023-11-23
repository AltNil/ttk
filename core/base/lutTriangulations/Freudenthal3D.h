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


    ttk::SimplexId getVertexEdgeNumber(const SimplexId& vertexId) const {
      return lutNumNeighbour3d[getCaseID(vertexId)];
    };

    int getVertexEdge ( 	const SimplexId &  	vertexId,
      const int &  	localEdgeId,
      SimplexId &  	edgeId 
    ) 		const {

      //todo get calculation procedure of edge
      
      int ex = extent[0];
      int ey = extent[1];
      int ez = extent[2];

      // Calculate the coordinates of the given vertexId
      auto xyDim = ex*ey;
      int z = vertexId/xyDim;
      int xy = (vertexId-z*xyDim);
      int y = xy/ex;
      int x = xy-y*ex;
      
      auto coordOffsets = lutNeghborOffset[getCaseID(vertexId)][localEdgeId];
      int lutIndex = coordOffsets[0]+1+3*(coordOffsets[1]+1)+9*(coordOffsets[2]+1); // generate unique id from coords

      int d1 = (ex-1)*ey*ez;
      int d2 = d1 + ex*(ey-1)*ez;
      int d3 = d2 + ex*ey*(ez-1);
      int d4 = d3 + (ex-1)*(ey-1)*ez;
      int d5 = d4 + ex*(ey-1)*(ez-1);
      int d6 = d5 + (ex-1)*ey*(ez-1);
      switch(lutIndex)
      {
        case 12:
          x -=1;
          lutIndex = 14;
          break;
        case 10:
          y -= 1;
          lutIndex = 16;
          break;
        case 4:
          z -= 1;
          lutIndex = 22;
          break;
        case 11:
          x += 1;
          y -= 1;
          lutIndex = 15;
          break;
        case 1:
          y -= 1;
          z -= 1;
          lutIndex = 25;
          break;
        case 5:
          x += 1;
          z -= 1;
          lutIndex = 21;
          break;
        case 2:
          x += 1;
          y -= 1;
          z -= 1;
          lutIndex = 24;
          break;
      }
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
        //edgeId = -1;
        edgeId = d6 + (x-1) + y * (ex-1) + z * (ex-1) * (ey-1);
        break;/*
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
        edgeId = d6 + (x-1) + (y-1) * ex + (z-1) * (ex-1) * (ey-1);
        break;*/
      
      default:
        return -1; //something went wrong
        break;
      }


      //todo calculate edge from coord
      /*switch (calcId)
      {
      case 0: // pdf case 1 positive
        edgeId = x + y*(ex-1)+(z*(ex-1)*ey);
        break;
      case 1: // pdf case 2 positive
        edgeId = 18 + x + (y*ex) + (z*ex*(ey-1));
        break;
      case 2: // pdf case 3 positive
        edgeId = 36 + x + (y*ex) + (z*ex*ey);
        break;
      case 3: // pdf case 4 positive
        edgeId = 54 + (x-1) + y*(ex-1) + z * (ex-1)*(ey-1);
        break;
      case 4: // pdf case 5 positive
        edgeId = 66 + x + y*ex + z * ex * (ey-1);
        break;
      case 5: //pdf case 6 positive
        edgeId = 78 + (x-1) + y * (ex-1) + z * (ex-1) * ey;
        break;
      case 6: // pdf case 7 positive
        edgeId = 90 + (x-1) + y * ex + z * (ex-1) * (ey-1);
        break;
      case 7: // pdf case 1 negative
        edgeId = (x-1)+y*(ex-1)+(z*(ex-1)*ey);
        break;
      case 8: // pdf case 2 negative
        edgeId = 18 + x + ((y-1)*ex) + (z*ex*(ey-1));
        break;
      case 9: // pdf case 3 negative
        edgeId = 36 + x + (y*ex) + ((z-1)*ex*ey);
        break;
      case 10: // pdf case 4 negative
        edgeId = 54 + (x) + (y-1)*(ex-1) + z * (ex-1)*(ey-1);
        break;
      case 11: // pdf case 5 negative
        edgeId = 66 + x + (y-1) * ex + (z-1) * ex * (ey-1);
        break;
      case 12: //pdf case 6 negative
        edgeId = 78 + (x) + y * (ex-1) + (z-1) * (ex-1) * ey;
        break;
      case 13: // pdf case 7 negative
        edgeId = 90 + (x) + (y-1) * ex + (z-1) * (ex-1) * (ey-1);
        break;
      
      default:
        return -1; //something went wrong
        break;
      }*/
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
    const std::vector<std::vector<int>> lutEdgeOffset = {
      {54, 68, 27, 58, 12, 46, -4, 36, 73, -5, 47, 15, 83, 63},
      {54, 68, 27, 58, 12, 46, -4, 15, 63, 73},
      {-5, 47, 15, 36, 73, 83, 63, 27, 54, 12},
      {-5, 47, 15, 36, 73, 83, 63, -4, 58, 27},
      {-4, 15, 73, 63, 27, 58},
      {-5, 47, 15, 36, 73, 83, 63, 27},
      {54, 68, 27, 58, 12, 46, -4, -5, 36, 73},
      {54, 68, 27, 58, 12, 46, -4, 73},
      {12, -5, 36, 73, 54, 27},
      {-5, 47, 15, 36, 73, 83, 63, 12, 46, -4},
      {12, 46, -4, 73, 15, 63},
      {-5, 47, 15, 36, 73, 83, 63, 12},
      {-5, 47, 15, 36, 73, 83, 63, -4},
      {-4, 15, 73, 63},
      {-5, 47, 15, 36, 73, 83, 63},
      {12, -5, 36, 73, 46, -4},
      {12, 46, -4, 73},
      {12, -5, 36, 73},
      {54, 68, 27, 58, 12, 46, -4, -5, 47, 15},
      {54, 68, 27, 58, 12, 46, -4, 15},
      {27, -5, 47, 15, 54, 12},
      {27, -5, 47, 15, 58, -4},
      {27, 58, -4, 15},
      {27, -5, 47, 15},
      {54, 68, 27, 58, 12, 46, -4, -5},
      {54, 68, 27, 58, 12, 46, -4},
      {54, 27, 12, -5}
    };
  };
}