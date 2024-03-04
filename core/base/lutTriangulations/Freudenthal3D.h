#include <Debug.h>
#include <AbstractTriangulation.h>
#include <iostream>
#include <vector>

namespace ttk {

class Freudenthal3D final : public ttk::AbstractTriangulation {

  public:

    Freudenthal3D(std::array<int,3> extent_) : extent(extent_){
      dx=extent[0];
      dxm1=extent[0]-1;
      dy=extent[1];
      dym1=extent[1]-1;
      dz=extent[2];
      dzm1=extent[2]-1;

      /* 
       * precalculate the implicite index offsets for neighbors
       */
      //std::vector<std::vector<int>> deltas;
      for(std::vector<std::array<int,3>> caseId : lutNeighborOffset){
        std::vector<int> shifts;
        for(std::array<int,3> coordDeltas: caseId){
          shifts.push_back(coordDeltas[0]+coordDeltas[1]*extent[0]+(extent[0]*extent[1]*coordDeltas[2]));
        }
        lutNeighborIndexOffset.push_back(shifts);
      }

      /*
       * precalculation for the edge calculations
       */
      edgeDimension1 = dxm1*dy*dz;
      edgeDimension2 = edgeDimension1 + dx*dym1*dz;
      edgeDimension3 = edgeDimension2 + dx*dy*dzm1;
      edgeDimension4 = edgeDimension3 + dxm1*dym1*dz;
      edgeDimension5 = edgeDimension4 + dx*dym1*dzm1;
      edgeDimension6 = edgeDimension5 + dxm1*dy*dzm1;

      /*
       * precalculations for the triangle calculations
       */
      //deltas.clear();
      for (std::vector<std::array<int,3>> caseId : lutTriangleOffset3d){
        std::vector<int> shifts;
        for (std::array<int,3> coordDeltas: caseId){
          shifts.push_back(coordDeltas[0]+coordDeltas[1]*extent[0]+(extent[0]*extent[1]*coordDeltas[2]));
        }
        lutTriangleIndexOffset.push_back(shifts);
      }
 dxm1*2*dym1*dz;
      plane2 = plane1 + dxm1*2*(dy)*dzm1;
      plane3 = plane2 + (dx*dym1*dzm1*2);
      plane4 = plane3 + dxm1*dym1*dzm1*2;
      plane5 = plane4 + dxm1*dym1*dzm1*2;

      /*
       * precalculations for the star calculations
       */
      //deltas.clear();
      for (std::vector<std::array<int,3>> caseId : lutTetraOffset3d){
        std::vector<int> shifts;
        for (std::array<int,3> coordDeltas: caseId){
          shifts.push_back(coordDeltas[0]+coordDeltas[1]*extent[0]+(extent[0]*extent[1]*coordDeltas[2]));
        }
        lutTetraIndexOffsets.push_back(shifts);
      }
      plane1 = dxm1*2*dym1*dz;
      plane2 = plane1 + dxm1*2*(dy)*dzm1;
      plane3 = plane2 + (dx*dym1*dzm1*2);
      plane4 = plane3 + dxm1*dym1*dzm1*2;
      plane5 = plane4 + dxm1*dym1*dzm1*2;

      /*
       * precalculations for the star calculations
       */
      //deltas.clear();
      for (std::vector<std::array<int,3>> caseId : lutTetraOffset3d){
        std::vector<int> shifts;
        for (std::array<int,3> coordDeltas: caseId){
          shifts.push_back(coordDeltas[0]+coordDeltas[1]*extent[0]+(extent[0]*extent[1]*coordDeltas[2]));
        }
        lutTetraIndexOffsets.push_back(shifts);
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
      return lutNumNeighbor3d[getCaseID(vertexId)];
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
      neighborId = SimplexId(vertexId+lutNeighborIndexOffset[getCaseID(vertexId)][localNeighborId]);
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
      return lutNumNeighbor3d[getCaseID(vertexId)];
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
      //int dx = extent[0];
      //int dy = extent[1];

      // Calculate the coordinates of the given vertexId
      auto xyDim = dx*dy;
      int z = vertexId/xyDim;
      int xy = (vertexId-z*xyDim);
      int y = xy/dx;
      int x = xy-y*dx;
      
      //search for the direction of the edge and tarnsform it to an id
      int lutIndex = lutEdgeDirection3d[getCaseID(vertexId)][localEdgeId]; // generate unique id from coords

      // calculate the explicit id of the searched edge with starting coordinates x, y and z from the starting point(given vertex)
      // the pdf cases refer to the different layer, in which the edges are enumerated through the grid.
      switch (lutIndex)
      {
      case 14: // pdf case 1 positive
        edgeId = x + y*dxm1+(z*dxm1*dy);
        break;
      case 16: // pdf case 2 positive
        edgeId = edgeDimension1 + x + (y*dx) + (z*dx*dym1);
        break;
      case 22: // pdf case 3 positive
        edgeId = edgeDimension2 + x + (y*dx) + (z*dx*dy);
        break;
      case 15: // pdf case 4 positive
        edgeId = edgeDimension3 + (x-1) + y*dxm1 + z * dxm1*dym1;
        break;
      case 25: // pdf case 5 positive
        edgeId = edgeDimension4 + x + y*dx + z * dx * dym1;
        break;
      case 21: //pdf case 6 positive
        edgeId = edgeDimension5 + (x-1) + y * dxm1 + z * dxm1 * dy;
        break;
      case 24: // pdf case 7 positive
        edgeId = edgeDimension6 + (x-1) + y * dxm1 + z * dxm1 * dym1;
        break;
      case 12: // pdf case 1 negative
        edgeId = (x-1)+y*dxm1+(z*dxm1*dy);
        break;
      case 10: // pdf case 2 negative
        edgeId = edgeDimension1 + x + ((y-1)*dx) + (z*dx*dym1);
        break;
      case 4: // pdf case 3 negative
        edgeId = edgeDimension2 + x + (y*dx) + ((z-1)*dx*dy);
        break;
      case 11: // pdf case 4 negative
        edgeId = edgeDimension3 + (x) + (y-1)*dxm1 + z * dxm1*dym1;
        break;
      case 1: // pdf case 5 negative
        edgeId = edgeDimension4 + x + (y-1) * dx + (z-1) * dx * dym1;
        break;
      case 5: //pdf case 6 negative
        edgeId = edgeDimension5 + (x) + y * dxm1 + (z-1) * dxm1 * dy;  
        break;
      case 2: // pdf case 7 negative
        edgeId = edgeDimension6 + (x) + (y-1) * dxm1 + (z-1) * dxm1 * dym1;
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
      return lutVertexTriangles3d[getCaseID(vertexId)];
    };

    int getVertexTriangle (
      const SimplexId &vertexId,
      const int &localTriangleId,
      SimplexId &trinagleId
    ) const {
      trinagleId = SimplexId(-1);
      //return 1;

      //int dx = extent[0];
      //int dy = extent[1];

      int caseId = getCaseID(vertexId);

      int lutIndex = lutTriangleDirection3d[caseId][localTriangleId];
      //std::array<int,3> vertexOffset = lutTriangleOffset3d[caseId][localTriangleId];

      int vId = vertexId + lutTriangleIndexOffset[caseId][localTriangleId];


      // Calculate the coordinates of the given vertexId
      auto xyDim = dx*dy;
      int z = (vId/xyDim);
      int xy = (vId-z*xyDim);
      int y = xy/dx;
      int x = xy-y*dx;

      // the cases are expecting the smallest point to be the origin of calculation
      int res = 0;

      switch (lutIndex){
        case 446:
          res = 2*(x + y*dxm1 + z*dxm1*dym1);
          break;
        case 447:
          res = 2*(x + y*dxm1 + z*dxm1*dym1)-1;
          break;

        case 608:
          res = plane1 + 2*(x + y*dxm1 + z*dxm1*(dy));
          break;
        case 609:
          res = plane1 + 2*(x + y*dxm1 + z*dxm1*(dy))-1;
          break;

        case 691: 
          res = plane2 + 2*(x + y*(dx) + z*(dx)*dym1);
          break;
        case 697:
          res = plane2 + 2*(x + y*(dx) + z*(dx)*dym1)+1;
          break;

        case 659:
          res = plane3 + 2*((x-1) + y*dxm1 + z*dxm1*dym1);
          break;
        case 669:
          res = plane3 + 2*((x-1) + y*dxm1 + z*dxm1*dym1)+1;
          break;

        case 689:
          res = plane4 + 2*((x) + y*dxm1 + z*dxm1*dym1);
          break;
        case 699:
          res = plane4 + 2*((x) + y*dxm1 + z*dxm1*dym1)-1;
          break;

        case 717:
          res = plane5 + 2*((x-1) + y*dxm1 + z*dxm1*dym1);
          break;
        case 724:
          res = plane5 + 2*((x-1) + y*dxm1 + z*dxm1*dym1)+1;
          break;
      
      default:
        // something went wrong
        return -1;
        break;
      }

      trinagleId = SimplexId(res);
      return 1;
    };


    ttk::SimplexId getVertexStarNumber(const ttk::SimplexId& vertexId) const final {
      return lutVertexTetrahedrons3d[getCaseID(vertexId)];
    };

    int getVertexStar(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId
    ) const final {
      //int dx = extent[0];
      //int dy = extent[1];

      int caseId = getCaseID(vertexId);
      
      int lutIndex = lutTetrahedronDirection3d[caseId][localStarId];
      int vId = vertexId + lutTetraIndexOffsets[caseId][localStarId];

      auto xyDim = dx*dy;
      int z = (vId/xyDim);
      int xy = (vId - z*xyDim);
      int y = xy/dx;
      int x = xy-y*dx;

      starId = 6*(x+y*dxm1+z*dxm1*dym1)+lutIndex;

      return 1;
    };

    int preconditionVertexStars() final {
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

    /*
    int getVertexId(const std::array<int,3>& coordinates)  const{
      auto result =  coordinates[0]+coordinates[1]*extent[0]+(extent[0]*extent[1]*coordinates[2]);
      return result;
    };*/
    int dx;
    int dxm1;
    int dy;
    int dym1;
    int dz;
    int dzm1;

    const std::array<int,3> extent;
    

    const int lutNumNeighbor3d[27] = {14,10,10,10,6,8,10,8,6,10,6,8,8,4,7,6,4,4,10,8,6,6,4,4,8,7,4};
    //std::array<int, 27> lutNumNeighbor3d = {6,5,5,5,4,4,5,4,4,5,4,4,4,3,3,4,3,3,5,4,4,4,3,3,4,3,3}; //option from quadrangulation to make comparisons
    
    const std::vector<std::vector<std::array<int,3>>> lutNeighborOffset = {
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
    std::vector<std::vector<int>> lutNeighborIndexOffset;

    // edge tables
    int edgeDimension1;
    int edgeDimension2;
    int edgeDimension3;
    int edgeDimension4;
    int edgeDimension5;
    int edgeDimension6;

    const std::vector<std::vector<int>> lutEdgeDirection3d = {
      {1, 2, 4, 5, 10, 11, 14, 21, 22, 12, 15, 16, 24, 25},
      {1, 2, 4, 5, 10, 11, 14, 16, 25, 22},
      {12, 15, 16, 21, 22, 24, 25, 4, 1, 10},
      {12, 15, 16, 21, 22, 24, 25, 14, 5, 4},
      {14, 16, 22, 25, 4, 5},
      {12, 15, 16, 21, 22, 24, 25, 4},
      {1, 2, 4, 5, 10, 11, 14, 12, 21, 22},
      {1, 2, 4, 5, 10, 11, 14, 22},
      {10, 12, 21, 22, 1, 4},
      {12, 15, 16, 21, 22, 24, 25, 10, 11, 14},
      {10, 11, 14, 22, 16, 25},
      {12, 15, 16, 21, 22, 24, 25, 10},
      {12, 15, 16, 21, 22, 24, 25, 14},
      {14, 16, 22, 25},
      {12, 15, 16, 21, 22, 24, 25},
      {10, 12, 21, 22, 11, 14},
      {10, 11, 14, 22},
      {10, 12, 21, 22},
      {1, 2, 4, 5, 10, 11, 14, 12, 15, 16},
      {1, 2, 4, 5, 10, 11, 14, 16},
      {4, 12, 15, 16, 1, 10},
      {4, 12, 15, 16, 5, 14},
      {4, 5, 14, 16},
      {4, 12, 15, 16},
      {1, 2, 4, 5, 10, 11, 14, 12},
      {1, 2, 4, 5, 10, 11, 14},
      {1, 4, 10, 12}
    };

    int plane1;
    int plane2;
    int plane3;
    int plane4;
    int plane5;

    const int lutVertexTriangles3d[27] = {36,21,21,21,9,15,21,15,9,21,9,15,15,5,12,9,5,5,21,15,9,9,5,5,15,12,5};
    //std::array<int, 27> lutVertexTriangles3d = {12,6,8,5,5,5,8,5,5,8,5,5,5,3,3,5,3,3,8,6,6,5,3,3,5,3,3}; //option from quadrangulation to make comparisons

    const std::vector<std::vector<int>> lutTriangleDirection3d = {
      {447, 717, 659, 699, 691, 697, 724, 609, 669, 689, 608, 446, 659, 669, 691, 697, 689, 699, 717, 724, 446, 447, 608, 609, 447, 659, 691, 669, 697, 446, 699, 609, 608, 689, 724, 717},
      {659, 669, 691, 697, 689, 699, 717, 724, 446, 447, 608, 609, 669, 697, 446, 608, 691, 689, 697, 717, 691},
      {447, 717, 659, 699, 691, 697, 724, 609, 669, 689, 608, 446, 447, 659, 691, 699, 691, 697, 609, 724, 697},
      {447, 717, 659, 699, 691, 697, 724, 609, 669, 689, 608, 446, 609, 724, 697, 608, 669, 446, 609, 608, 689},
      {608, 669, 697, 446, 609, 608, 691, 689, 697},
      {609, 446, 447, 724, 697, 717, 659, 699, 691, 697, 724, 609, 669, 689, 608},
      {659, 669, 691, 697, 689, 699, 717, 724, 446, 447, 608, 609, 699, 447, 609, 659, 608, 691, 609, 717, 608},
      {659, 669, 691, 697, 689, 699, 717, 724, 446, 447, 608, 609, 717, 608, 691},
      {699, 447, 691, 697, 609, 659, 608, 691, 609},
      {447, 717, 659, 699, 691, 697, 724, 609, 669, 689, 608, 446, 447, 659, 691, 446, 447, 717, 608, 446, 689},
      {446, 447, 717, 608, 691, 446, 691, 689, 697},
      {447, 659, 608, 691, 609, 447, 717, 659, 699, 691, 697, 724, 669, 689, 446},
      {447, 717, 659, 699, 691, 697, 724, 609, 669, 689, 608, 446, 446, 608, 689},
      {446, 608, 691, 689, 697},
      {447, 717, 659, 699, 691, 697, 724, 609, 669, 689, 608, 446},
      {447, 659, 608, 691, 609, 446, 447, 717, 608},
      {446, 447, 717, 608, 691},
      {447, 659, 608, 691, 609},
      {659, 669, 691, 697, 689, 699, 717, 724, 446, 447, 608, 609, 699, 447, 609, 446, 447, 724, 697, 669, 446},
      {659, 669, 691, 697, 689, 699, 717, 724, 446, 447, 608, 609, 669, 697, 446},
      {699, 447, 691, 697, 609, 446, 447, 724, 697},
      {609, 446, 447, 724, 697, 608, 669, 446, 609},
      {608, 669, 697, 446, 609},
      {609, 446, 447, 724, 697},
      {699, 447, 691, 697, 609, 659, 669, 689, 699, 717, 724, 446, 447, 608, 609},
      {659, 669, 691, 697, 689, 699, 717, 724, 446, 447, 608, 609},
      {699, 447, 691, 697, 609}
    };

    const std::vector<std::vector<std::array<int,3>>> lutTriangleOffset3d = {
      {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {1, -1, -1}, {1, -1, -1}, {0, -1, -1}, {0, -1, -1}, {0, -1, -1}, {1, -1, -1}, {1, -1, -1}, {1, -1, -1}, {0, -1, 0}, {1, -1, 0}, {0, 0, -1}, {1, 0, -1}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {1, 0, -1}, {0, 0, -1}, {0, 0, 0}, {0, -1, -1}, {0, 0, -1}, {0, 0, 0}, {0, 0, 0}, {0, 0, -1}, {1, -1, 0}},
      {{1, -1, -1}, {1, -1, -1}, {0, -1, -1}, {0, -1, -1}, {0, -1, -1}, {1, -1, -1}, {1, -1, -1}, {1, -1, -1}, {0, -1, 0}, {1, -1, 0}, {0, 0, -1}, {1, 0, -1}, {1, 0, -1}, {0, 0, -1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {1, -1, 0}, {0, -1, 0}},
      {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, -1}, {0, -1, -1}, {0, -1, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}},
      {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 0, 0}, {1, 0, -1}, {0, 0, 0}, {0, 0, 0}},
      {{0, 0, -1}, {1, 0, -1}, {0, 0, -1}, {0, 0, 0}, {1, 0, -1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
      {{0, 0, -1}, {-1, 0, 0}, {0, 0, 0}, {0, 0, -1}, {0, 0, -1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {-1, 0, 0}, {-1, 0, 0}},
      {{1, -1, -1}, {1, -1, -1}, {0, -1, -1}, {0, -1, -1}, {0, -1, -1}, {1, -1, -1}, {1, -1, -1}, {1, -1, -1}, {0, -1, 0}, {1, -1, 0}, {0, 0, -1}, {1, 0, -1}, {0, -1, -1}, {0, -1, 0}, {0, 0, -1}, {0, -1, 0}, {-1, 0, 0}, {0, -1, 0}, {0, 0, 0}, {1, -1, 0}, {0, 0, 0}},
      {{1, -1, -1}, {1, -1, -1}, {0, -1, -1}, {0, -1, -1}, {0, -1, -1}, {1, -1, -1}, {1, -1, -1}, {1, -1, -1}, {0, -1, 0}, {1, -1, 0}, {0, 0, -1}, {1, 0, -1}, {1, -1, 0}, {0, 0, 0}, {0, -1, 0}},
      {{0, -1, -1}, {0, -1, 0}, {0, -1, -1}, {0, -1, -1}, {0, 0, -1}, {0, -1, 0}, {-1, 0, 0}, {0, -1, 0}, {0, 0, 0}},
      {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {1, -1, 0}, {1, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
      {{0, -1, 0}, {1, -1, 0}, {1, -1, 0}, {0, 0, 0}, {0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
      {{0, -1, 0}, {0, -1, 0}, {-1, 0, 0}, {0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {-1, 0, 0}, {-1, 0, 0}},
      {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
      {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
      {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {-1, 0, 0}},
      {{0, -1, 0}, {0, -1, 0}, {-1, 0, 0}, {0, -1, 0}, {0, 0, 0}, {0, -1, 0}, {1, -1, 0}, {1, -1, 0}, {0, 0, 0}},
      {{0, -1, 0}, {1, -1, 0}, {1, -1, 0}, {0, 0, 0}, {0, -1, 0}},
      {{0, -1, 0}, {0, -1, 0}, {-1, 0, 0}, {0, -1, 0}, {0, 0, 0}},
      {{1, -1, -1}, {1, -1, -1}, {0, -1, -1}, {0, -1, -1}, {0, -1, -1}, {1, -1, -1}, {1, -1, -1}, {1, -1, -1}, {0, -1, 0}, {1, -1, 0}, {0, 0, -1}, {1, 0, -1}, {0, -1, -1}, {0, -1, 0}, {0, 0, -1}, {-1, 0, 0}, {0, 0, 0}, {0, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 0, 0}},
      {{1, -1, -1}, {1, -1, -1}, {0, -1, -1}, {0, -1, -1}, {0, -1, -1}, {1, -1, -1}, {1, -1, -1}, {1, -1, -1}, {0, -1, 0}, {1, -1, 0}, {0, 0, -1}, {1, 0, -1}, {1, 0, -1}, {0, 0, -1}, {0, 0, 0}},
      {{0, -1, -1}, {0, -1, 0}, {0, -1, -1}, {0, -1, -1}, {0, 0, -1}, {-1, 0, 0}, {0, 0, 0}, {0, 0, -1}, {0, 0, -1}},
      {{0, 0, -1}, {-1, 0, 0}, {0, 0, 0}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 0, 0}, {1, 0, -1}},
      {{0, 0, -1}, {1, 0, -1}, {0, 0, -1}, {0, 0, 0}, {1, 0, -1}},
      {{0, 0, -1}, {-1, 0, 0}, {0, 0, 0}, {0, 0, -1}, {0, 0, -1}},
      {{0, -1, -1}, {0, -1, 0}, {0, -1, -1}, {0, -1, -1}, {0, 0, -1}, {1, -1, -1}, {1, -1, -1}, {0, -1, -1}, {1, -1, -1}, {1, -1, -1}, {1, -1, -1}, {0, -1, 0}, {1, -1, 0}, {0, 0, -1}, {1, 0, -1}},
      {{1, -1, -1}, {1, -1, -1}, {0, -1, -1}, {0, -1, -1}, {0, -1, -1}, {1, -1, -1}, {1, -1, -1}, {1, -1, -1}, {0, -1, 0}, {1, -1, 0}, {0, 0, -1}, {1, 0, -1}},
      {{0, -1, -1}, {0, -1, 0}, {0, -1, -1}, {0, -1, -1}, {0, 0, -1}}
    };
    std::vector<std::vector<int>> lutTriangleIndexOffset;
  
    const int lutVertexTetrahedrons3d[27] = {24,12,12,12,4,8,12,8,4,12,4,8,8,2,6,4,2,2,12,8,4,4,2,2,8,6,2};
    //std::array<int,27> lutVertexTetrahedrons3d = {8,4,4,4,2,2,4,2,2,4,2,2,2,1,1,2,1,1,4,2,2,2,1,1,2,1,1}; //option from quadrangulation to make comparisons

    const std::vector<std::vector<int>> lutTetrahedronDirection3d = {
      {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 2, 0, 1, 1, 5, 2, 3, 3, 4, 4, 5},
      {0, 1, 2, 3, 4, 5, 0, 1, 0, 2, 2, 3},
      {0, 1, 2, 3, 4, 5, 3, 4, 4, 5, 1, 5},
      {0, 1, 2, 3, 4, 5, 0, 2, 3, 4, 2, 3},
      {2, 3, 0, 2},
      {0, 1, 2, 3, 4, 5, 3, 4},
      {0, 1, 2, 3, 4, 5, 4, 5, 0, 1, 1, 5},
      {0, 1, 2, 3, 4, 5, 0, 1},
      {4, 5, 1, 5},
      {0, 1, 2, 3, 4, 5, 1, 5, 0, 1, 0, 2},
      {0, 1, 0, 2},
      {0, 1, 2, 3, 4, 5, 1, 5},
      {0, 1, 2, 3, 4, 5, 0, 2},
      {0, 2},
      {0, 1, 2, 3, 4, 5},
      {1, 5, 0, 1},
      {0, 1},
      {1, 5},
      {0, 1, 2, 3, 4, 5, 4, 5, 3, 4, 2, 3},
      {0, 1, 2, 3, 4, 5, 2, 3},
      {4, 5, 3, 4},
      {3, 4, 2, 3},
      {2, 3},
      {3, 4},
      {0, 1, 2, 3, 4, 5, 4, 5},
      {0, 1, 2, 3, 4, 5},
      {4, 5}
    };

    const std::vector<std::vector<std::array<int,3>>> lutTetraOffset3d= {
      {{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,0,0},{0,0,0},{0,-1,0},{0,-1,0},{-1,-1,0},{-1,-1,0},{0,0,-1},{0,0,-1},{-1,0,-1},{-1,0,-1},{-1,-1,-1},{-1,-1,-1}},
      {{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,0},{0,-1,0},{0,0,0},{0,0,0},{0,0,-1},{0,0,-1}},
      {{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,-1},{-1,0,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,0},{-1,-1,0}},
      {{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{0,0,0},{0,0,0},{-1,0,-1},{-1,0,-1},{0,0,-1},{0,0,-1}},
      {{0,0,-1},{0,0,-1},{0,0,0},{0,0,0}},
      {{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,-1},{-1,0,-1}},
      {{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{-1,-1,-1},{-1,-1,-1},{0,-1,0},{0,-1,0},{-1,-1,0},{-1,-1,0}},
      {{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,0},{0,-1,0}},
      {{-1,-1,-1},{-1,-1,-1},{-1,-1,0},{-1,-1,0}},
      {{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,-1,0},{-1,-1,0},{0,-1,0},{0,-1,0},{0,0,0},{0,0,0}},
      {{0,-1,0},{0,-1,0},{0,0,0},{0,0,0}},
      {{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,-1,0},{-1,-1,0}},
      {{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0}},
      {{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0},{-1,0,0}},
      {{-1,-1,0},{-1,-1,0},{0,-1,0},{0,-1,0}},
      {{0,-1,0},{0,-1,0}},
      {{-1,-1,0},{-1,-1,0}},
      {{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,0,-1},{-1,0,-1},{0,0,-1},{0,0,-1}},
      {{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,0,-1},{0,0,-1}},
      {{-1,-1,-1},{-1,-1,-1},{-1,0,-1},{-1,0,-1}},
      {{-1,0,-1},{-1,0,-1},{0,0,-1},{0,0,-1}},
      {{0,0,-1},{0,0,-1}},
      {{-1,0,-1},{-1,0,-1}},
      {{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{-1,-1,-1},{-1,-1,-1}},
      {{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,-1}},
      {{-1,-1,-1},{-1,-1,-1}},
    };
    std::vector<std::vector<int>> lutTetraIndexOffsets;
  };
}