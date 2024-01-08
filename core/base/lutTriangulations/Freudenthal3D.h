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
      for(std::vector<std::array<int,3>> caseId : lutNeighborOffset){
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
      int ex = extent[0];
      int ey = extent[1];

      // Calculate the coordinates of the given vertexId
      auto xyDim = ex*ey;
      int z = vertexId/xyDim;
      int xy = (vertexId-z*xyDim);
      int y = xy/ex;
      int x = xy-y*ex;
      
      //search for the direction of the edge and tarnsform it to an id
      int lutIndex = lutEdgeDirection3d[getCaseID(vertexId)][localEdgeId]; // generate unique id from coords

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
      return lutVertexTriangles3d[getCaseID(vertexId)];
    };
    
    int getVertexTriangle (
      const SimplexId &vertexId,
      const int &localTriangleId,
      SimplexId &trinagleId
    ) const {
      trinagleId = SimplexId(-1);
      //return 1;

      int dx = extent[0];
      int dy = extent[1];
      int dz = extent[2];

      int plane1 = (dx-1)*2*(dy-1)*dz;
      int plane2 = plane1 + (dx-1)*2*(dy)*(dz-1);
      int plane3 = plane2 + (dx*(dy-1)*(dz-1)*2);
      int plane4 = plane3 + (dx-1)*(dy-1)*(dz-1)*2;
      int plane5 = plane4 + (dx-1)*(dy-1)*(dz-1)*2;


      int caseId = getCaseID(vertexId);

      int lutIndex = lutTriangleDirection3d[caseId][localTriangleId];
      std::array<int,3> vertexOffset = lutTriangleOffset3d[caseId][localTriangleId];

      int x = vertexOffset[0];
      int y = vertexOffset[1];
      int z = vertexOffset[2];


      // Calculate the coordinates of the given vertexId
      auto xyDim = dx*dy;
      int cz = (vertexId/xyDim);
      int xy = (vertexId-cz*xyDim);
      int cy = xy/dx;
      int cx = xy-cy*dx;
      // calculate coordinates from id

      x = x + cx;
      y = y + cy;
      z = z + cz;

      // the cases are expecting the smallest point to be the origin of calculation
      int res = 0;

      switch (lutIndex){
        case 446:
          res = x*2 + 2*y*(dx-1)+2*z*(dx-1)*(dy-1);
          break;
        case 447:
          res = (x*2 + 2*y*(dx-1)+2*z*(dx-1)*(dy-1))-1;
          break;

        case 608:
          res = plane1+x*2 + 2*y*(dx-1)+2*z*(dx-1)*(dy);
          break;
        case 609:
          res = plane1+(x*2 + 2*y*(dx-1)+2*z*(dx-1)*(dy))-1;
          break;

        case 691: 
          res = plane2+x*2 + 2*y*(dx)+2*z*(dx)*(dy-1);
          break;
        case 697:
          res = plane2+(x*2 + 2*y*(dx)+2*z*(dx)*(dy-1))+1;
          break;

        case 659:
          res = plane3+(2*(x-1)+2*y*(dx-1)+2*z*(dx-1)*(dy-1));
          break;
        case 669:
          res = plane3+(2*(x-1)+2*y*(dx-1)+2*z*(dx-1)*(dy-1))+1;
          break;

        case 689:
          res = plane4+(2*(x)+2*y*(dx-1)+2*z*(dx-1)*(dy-1));
          break;
        case 699:
          res = plane4+(2*(x)+2*y*(dx-1)+2*z*(dx-1)*(dy-1))-1;
          break;

        case 717:
          res = plane5+(2*(x-1)+2*y*(dx-1)+2*z*(dx-1)*(dy-1));
          break;
        case 724:
          res = plane5+(2*(x-1)+2*y*(dx-1)+2*z*(dx-1)*(dy-1))+1;
          break;
      
      default:
        return -1;
        break;
      }

      trinagleId = SimplexId(res);
      return 1;
    };


    ttk::SimplexId getVertexStarNumber(const ttk::SimplexId& vertexId) const final {
      return 0;
    };

    int getVertexStar(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId
    ) const final {
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
    
    std::vector<std::vector<int>> lutIndexOffset;
    int d1;
    int d2;
    int d3;
    int d4;
    int d5;
    int d6;

    const int lutNumNeighbor3d[27] = {14,10,10,10,6,8,10,8,6,10,6,8,8,4,7,6,4,4,10,8,6,6,4,4,8,7,4};

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

    const int lutVertexTriangles3d[27] = {36,21,21,21,9,15,21,15,9,21,9,15,15,5,12,9,5,5,21,15,9,9,5,5,15,12,5};

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
  


    const std::vectror<std::vector<std::array<int,3>>> lutTetraOffset3d= {
      {{-1,0,0},{0,0,0},{-1,0,0},{0,0,0},{0,0,0},{0,0,0},{0,-1,-1},{1,-1,-1},{0,-1,-1},{1,-1,-1},{1,-1,-1},{1,-1,-1},{0,0,0},{0,0,0},{0,-1,0},{1,-1,0},{0,-1,0},{0,-1,0},{0,-1,0},{1,0,-1},{0,0,-1},{0,0,-1},{0,-1,-1},{0,-1,-1}},
      {{0,-1,-1},{1,-1,-1},{0,-1,-1},{1,-1,-1},{1,-1,-1},{1,-1,-1},{0,-1,0},{1,-1,0},{0,0,0},{0,0,0},{0,-1,0},{1,0,-1}},
      {{-1,0,0},{0,0,0},{-1,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,-1},{0,0,-1},{0,-1,-1},{0,-1,-1},{0,-1,0},{1,-1,-1}},
      {{-1,0,0},{0,0,0},{-1,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,-1},{0,0,-1},{0,-1,0},{-1,0,0}},
      {{0,-1,0},{1,0,-1},{0,0,0},{0,0,0}},
      {{-1,0,0},{0,0,0},{-1,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,-1},{0,0,-1}},
      {{0,-1,-1},{1,-1,-1},{0,-1,-1},{1,-1,-1},{1,-1,-1},{1,-1,-1},{0,-1,-1},{0,-1,-1},{0,-1,0},{1,-1,0},{0,-1,0},{0,-1,0}},
      {{0,-1,-1},{1,-1,-1},{0,-1,-1},{1,-1,-1},{1,-1,-1},{1,-1,-1},{0,-1,0},{1,-1,0}},
      {{0,-1,-1},{0,-1,-1},{0,-1,0},{0,-1,0}},
      {{-1,0,0},{0,0,0},{-1,0,0},{0,0,0},{0,0,0},{0,0,0},{0,-1,0},{0,-1,0},{0,-1,0},{1,-1,0},{0,0,0},{0,0,0}},
      {{0,-1,0},{1,-1,0},{0,0,0},{0,0,0}},
      {{-1,0,0},{0,0,0},{-1,0,0},{0,0,0},{0,0,0},{0,0,0},{0,-1,0},{0,-1,0}},
      {{-1,0,0},{0,0,0},{-1,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
      {{0,0,0},{0,0,0}},
      {{-1,0,0},{0,0,0},{-1,0,0},{0,0,0},{0,0,0},{0,0,0}},
      {{0,-1,0},{0,-1,0},{0,-1,0},{1,-1,0}},
      {{0,-1,0},{1,-1,0}},
      {{0,0,-1},{0,-1,0}},
      {{0,-1,-1},{1,-1,-1},{0,-1,-1},{1,-1,-1},{1,-1,-1},{1,-1,-1},{0,-1,-1},{0,-1,-1},{0,0,-1},{0,0,-1},{0,-1,0},{1,0,-1}},
      {{0,-1,-1},{1,-1,-1},{0,-1,-1},{1,-1,-1},{1,-1,-1},{1,-1,-1},{0,-1,0},{1,0,-1}},
      {{0,-1,-1},{0,-1,-1},{0,0,-1},{0,0,-1}},
      {{0,0,-1},{0,0,-1},{0,-1,0},{1,0,-1}},
      {{0,0,0},{0,0,0}},
      {{0,0,-1},{0,0,-1}},
      {{0,-1,-1},{1,-1,-1},{0,-1,-1},{1,-1,-1},{1,-1,-1},{1,-1,-1},{0,-1,-1},{0,-1,-1}},
      {{0,-1,-1},{1,-1,-1},{0,-1,-1},{1,-1,-1},{1,-1,-1},{1,-1,-1}},
      {{0,-1,-1},{0,-1,-1}}
    };
  };
}