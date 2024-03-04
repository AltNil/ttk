/// \ingroup base
/// \class ttk::ZombieTriangulation
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date January 2016.
///
/// \brief ZombieTriangulation is a class that provides time and memory
/// efficient traversal methods on triangulations of piecewise linear
/// manifolds represented by regular grids.
/// \sa Triangulation

#pragma once

#include <array>

// base code includes
#include <RegularGridTriangulation.h>

namespace ttk {

  class ZombieTriangulation : public RegularGridTriangulation {

  public:
    ZombieTriangulation();
    ~ZombieTriangulation() override;

    ZombieTriangulation(const ZombieTriangulation &) = default;
    ZombieTriangulation(ZombieTriangulation &&) = default;
    ZombieTriangulation &operator=(const ZombieTriangulation &) = default;
    ZombieTriangulation &operator=(ZombieTriangulation &&) = default;
    ZombieTriangulation(std::array<int,3> extent_) : extent(extent_){
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

    inline const std::array<SimplexId, 3> &getGridDimensions() const override {
      return this->dimensions_;
    }

    int getCellEdgeInternal(const SimplexId &cellId,
                            const int &id,
                            SimplexId &edgeId) const override;

    SimplexId getCellEdgeNumberInternal(const SimplexId &cellId) const override;

    const std::vector<std::vector<SimplexId>> *getCellEdgesInternal() override;

    int TTK_TRIANGULATION_INTERNAL(getCellNeighbor)(
      const SimplexId &cellId,
      const int &localNeighborId,
      SimplexId &neighborId) const override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getCellNeighborNumber)(
      const SimplexId &cellId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getCellNeighbors)() override;

    int getCellTriangleInternal(const SimplexId &cellId,
                                const int &id,
                                SimplexId &triangleId) const override;

    SimplexId getCellTriangleNumberInternal(
      const SimplexId & /*cellId*/) const override {
      // NOTE: the output is always 4 here. let's keep the function in there
      // in case of further generalization to CW-complexes
      return 4;
    }

    const std::vector<std::vector<SimplexId>> *
      getCellTrianglesInternal() override;

    int TTK_TRIANGULATION_INTERNAL(getCellVertex)(
      const SimplexId &cellId,
      const int &localVertexId,
      SimplexId &vertexId) const override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getCellVertexNumber)(
      const SimplexId &cellId) const override;

    int TTK_TRIANGULATION_INTERNAL(getDimensionality)() const override {
      return dimensionality_;
    }

    SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeLinkNumber)(
      const SimplexId &edgeId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeLinks)() override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getEdgeStars)() override;

    const std::vector<std::vector<SimplexId>> *
      getEdgeTrianglesInternal() override;

    const std::vector<std::array<SimplexId, 2>> *
      TTK_TRIANGULATION_INTERNAL(getEdges)() override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getNumberOfCells)() const override {
      return cellNumber_;
    }

    SimplexId getNumberOfEdgesInternal() const override {
      return edgeNumber_;
    }

    SimplexId getNumberOfTrianglesInternal() const override {
      return triangleNumber_;
    }

    SimplexId TTK_TRIANGULATION_INTERNAL(getNumberOfVertices)() const override {
      return vertexNumber_;
    }

    virtual int getTetrahedronEdge(const SimplexId &tetId,
                                   const int &id,
                                   SimplexId &edgeId) const = 0;

    int getTetrahedronEdges(std::vector<std::vector<SimplexId>> &edges) const;

    virtual int getTetrahedronTriangle(const SimplexId &tetId,
                                       const int &id,
                                       SimplexId &triangleId) const = 0;

    int getTetrahedronTriangles(
      std::vector<std::vector<SimplexId>> &triangles) const;

    virtual int getTetrahedronNeighbor(const SimplexId &tetId,
                                       const int &localNeighborId,
                                       SimplexId &neighborId) const = 0;

    virtual SimplexId
      getTetrahedronNeighborNumber(const SimplexId &tetId) const = 0;

    int getTetrahedronNeighbors(std::vector<std::vector<SimplexId>> &neighbors);

    virtual int getTetrahedronVertex(const SimplexId &tetId,
                                     const int &localVertexId,
                                     SimplexId &vertexId) const = 0;

    SimplexId getTriangleEdgeNumberInternal(
      const SimplexId & /*triangleId*/) const override {
      // NOTE: the output is always 3 here. let's keep the function in there
      // in case of further generalization to CW-complexes
      return 3;
    }

    const std::vector<std::vector<SimplexId>> *
      getTriangleEdgesInternal() override;

    int getTriangleEdgesInternal(
      std::vector<std::vector<SimplexId>> &edges) const;

    SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleLinkNumber)(
      const SimplexId &triangleId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleLinks)() override;

    virtual int getTriangleNeighbor(const SimplexId &triangleId,
                                    const int &localNeighborId,
                                    SimplexId &neighborId) const = 0;

    virtual SimplexId
      getTriangleNeighborNumber(const SimplexId &triangleId) const = 0;

    int getTriangleNeighbors(std::vector<std::vector<SimplexId>> &neighbors);

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getTriangleStars)() override;

    const std::vector<std::array<SimplexId, 3>> *
      TTK_TRIANGULATION_INTERNAL(getTriangles)() override;

    SimplexId
      getVertexEdgeNumberInternal(const SimplexId &vertexId) const override {
      // return the neigbor by the caseID of the given vertex
      return lutNumNeighbor3d[getCaseID(vertexId)];
    };;

    const std::vector<std::vector<SimplexId>> *
      getVertexEdgesInternal() override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getVertexLinkNumber)(
      const SimplexId &vertexId) const override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexLinks)() override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexNeighbors)() override;

    const std::vector<std::vector<SimplexId>> *
      TTK_TRIANGULATION_INTERNAL(getVertexStars)() override;

    const std::vector<std::vector<SimplexId>> *
      getVertexTrianglesInternal() override;

    inline bool isEmpty() const override {
      return !vertexNumber_;
    }

    bool TTK_TRIANGULATION_INTERNAL(isTriangleOnBoundary)(
      const SimplexId &triangleId) const override;

    int setInputGrid(const float &xOrigin,
                     const float &yOrigin,
                     const float &zOrigin,
                     const float &xSpacing,
                     const float &ySpacing,
                     const float &zSpacing,
                     const SimplexId &xDim,
                     const SimplexId &yDim,
                     const SimplexId &zDim) override;

    virtual int preconditionVerticesInternal() = 0;
    int preconditionVertexNeighborsInternal() override;
    int preconditionEdgesInternal() override = 0;
    int preconditionTrianglesInternal() override = 0;
    virtual int preconditionTetrahedronsInternal() = 0;

    inline int preconditionCellsInternal() {
      if(dimensionality_ == 3) {
        return this->preconditionTetrahedronsInternal();
      } else if(dimensionality_ == 2 && !hasPreconditionedTriangles_) {
        hasPreconditionedTriangles_ = true;
        return this->preconditionTrianglesInternal();
      }
      return 0;
    }

    inline int preconditionVerticesAndCells() {
      if(!this->hasPreconditionedVerticesAndCells_) {
        this->preconditionVerticesInternal();
        this->preconditionCellsInternal();
        this->hasPreconditionedVerticesAndCells_ = true;
      }
      return 0;
    }

    inline int getCellVTKIDInternal(const int &ttkId,
                                    int &vtkId) const override {
#ifndef TTK_ENABLE_KAMIKAZE
      if(ttkId < 0) {
        return -1;
      }
#endif // TTK_ENABLE_KAMIKAZE
      const SimplexId nSimplexPerCell{this->getDimensionality() == 3 ? 6 : 2};
      vtkId = ttkId / nSimplexPerCell;
      return 0;
    }

#ifdef TTK_ENABLE_MPI

  protected:
    int preconditionDistributedCells() override;

  public:
    void createMetaGrid(const double *const bounds) override;
    int getCellRankInternal(const SimplexId lcid) const override;

  protected:
    bool isVertexOnGlobalBoundaryInternal(const SimplexId lvid) const override;
    bool isEdgeOnGlobalBoundaryInternal(const SimplexId leid) const override;
    bool
      isTriangleOnGlobalBoundaryInternal(const SimplexId ltid) const override;

  private:
    std::array<SimplexId, 3>
      getVertGlobalCoords(const SimplexId lvid) const override;
    std::array<SimplexId, 3>
      getVertLocalCoords(const SimplexId gvid) const override;

#endif // TTK_ENABLE_MPI

  protected:
    enum class VertexPosition : char {
      // a--------b

      LEFT_CORNER_1D, // a
      RIGHT_CORNER_1D, // b
      CENTER_1D,
      // total: 3 1D cases

      // a--------b
      // |        |
      // |        |
      // |        |
      // c--------d

      // 2D corners
      TOP_LEFT_CORNER_2D, // a
      TOP_RIGHT_CORNER_2D, // b
      BOTTOM_LEFT_CORNER_2D, // c
      BOTTOM_RIGHT_CORNER_2D, // d
      // 2D edges
      TOP_EDGE_2D, // ab
      BOTTOM_EDGE_2D, // cd
      LEFT_EDGE_2D, // ac
      RIGHT_EDGE_2D, // bd
      // 2D central strip
      CENTER_2D,
      // total: 9 2D cases

      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--------b  |
      // |  g-----|--h
      // | /      | /
      // |/       |/
      // c--------d

      // 3D corners
      TOP_LEFT_FRONT_CORNER_3D, // a
      TOP_RIGHT_FRONT_CORNER_3D, // b
      BOTTOM_LEFT_FRONT_CORNER_3D, // c
      BOTTOM_RIGHT_FRONT_CORNER_3D, // d
      TOP_LEFT_BACK_CORNER_3D, // e
      TOP_RIGHT_BACK_CORNER_3D, // f
      BOTTOM_LEFT_BACK_CORNER_3D, // g
      BOTTOM_RIGHT_BACK_CORNER_3D, // h
      // 3D edges
      TOP_FRONT_EDGE_3D, // ab
      BOTTOM_FRONT_EDGE_3D, // cd
      LEFT_FRONT_EDGE_3D, // ac
      RIGHT_FRONT_EDGE_3D, // bd
      TOP_BACK_EDGE_3D, // ef
      BOTTOM_BACK_EDGE_3D, // gh
      LEFT_BACK_EDGE_3D, // eg
      RIGHT_BACK_EDGE_3D, // fh
      TOP_LEFT_EDGE_3D, // ae
      TOP_RIGHT_EDGE_3D, // bf
      BOTTOM_LEFT_EDGE_3D, // cg
      BOTTOM_RIGHT_EDGE_3D, // dh
      // 3D faces
      FRONT_FACE_3D, // abcd
      BACK_FACE_3D, // efgh
      TOP_FACE_3D, // abef
      BOTTOM_FACE_3D, // cdgh
      LEFT_FACE_3D, // aceg
      RIGHT_FACE_3D, // bdfh
      // 3D central part
      CENTER_3D,
      // total: 27 3D cases
    };

    // vertex neighbor shifts
    std::array<SimplexId, 14> vertexNeighborABCDEFGH_{};

    std::array<SimplexId, 10> vertexNeighborABCD_{};
    std::array<SimplexId, 10> vertexNeighborEFGH_{};
    std::array<SimplexId, 10> vertexNeighborAEFB_{};
    std::array<SimplexId, 10> vertexNeighborGHDC_{};
    std::array<SimplexId, 10> vertexNeighborAEGC_{};
    std::array<SimplexId, 10> vertexNeighborBFHD_{};

    std::array<SimplexId, 8> vertexNeighborAB_{};
    std::array<SimplexId, 8> vertexNeighborBD_{};
    std::array<SimplexId, 8> vertexNeighborGH_{};
    std::array<SimplexId, 8> vertexNeighborEG_{};
    std::array<SimplexId, 8> vertexNeighborCG_{};
    std::array<SimplexId, 8> vertexNeighborBF_{};

    std::array<SimplexId, 7> vertexNeighborB_{};
    std::array<SimplexId, 7> vertexNeighborG_{};

    std::array<SimplexId, 6> vertexNeighborEF_{};
    std::array<SimplexId, 6> vertexNeighborCD_{};
    std::array<SimplexId, 6> vertexNeighborAC_{};
    std::array<SimplexId, 6> vertexNeighborAE_{};
    std::array<SimplexId, 6> vertexNeighborFH_{};
    std::array<SimplexId, 6> vertexNeighborDH_{};

    std::array<SimplexId, 4> vertexNeighborA_{};
    std::array<SimplexId, 4> vertexNeighborC_{};
    std::array<SimplexId, 4> vertexNeighborD_{};
    std::array<SimplexId, 4> vertexNeighborE_{};
    std::array<SimplexId, 4> vertexNeighborF_{};
    std::array<SimplexId, 4> vertexNeighborH_{};

    std::array<SimplexId, 6> vertexNeighbor2dABCD_{};
    std::array<SimplexId, 4> vertexNeighbor2dAB_{};
    std::array<SimplexId, 4> vertexNeighbor2dCD_{};
    std::array<SimplexId, 4> vertexNeighbor2dAC_{};
    std::array<SimplexId, 4> vertexNeighbor2dBD_{};
    std::array<SimplexId, 3> vertexNeighbor2dB_{};
    std::array<SimplexId, 3> vertexNeighbor2dC_{};
    std::array<SimplexId, 2> vertexNeighbor2dA_{};
    std::array<SimplexId, 2> vertexNeighbor2dD_{};

    enum class EdgePosition : char {
      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--------b  |
      // |  g-----|--h
      // | /      | /
      // |/       |/
      // c--------d

      // length (ab)
      L_xnn_3D,
      L_xn0_3D,
      L_xnN_3D,
      L_x0n_3D,
      L_x00_3D,
      L_x0N_3D,
      L_xNn_3D,
      L_xN0_3D,
      L_xNN_3D,
      // height (ac)
      H_nyn_3D,
      H_ny0_3D,
      H_nyN_3D,
      H_0yn_3D,
      H_0y0_3D,
      H_0yN_3D,
      H_Nyn_3D,
      H_Ny0_3D,
      H_NyN_3D,
      // depth (ae)
      P_nnz_3D,
      P_n0z_3D,
      P_nNz_3D,
      P_0nz_3D,
      P_00z_3D,
      P_0Nz_3D,
      P_Nnz_3D,
      P_N0z_3D,
      P_NNz_3D,
      // diagonal1 (bc)
      D1_xyn_3D,
      D1_xy0_3D,
      D1_xyN_3D,
      // diagonal2 (ag)
      D2_nyz_3D,
      D2_0yz_3D,
      D2_Nyz_3D,
      // diagonal3 (be)
      D3_xnz_3D,
      D3_x0z_3D,
      D3_xNz_3D,
      // diagonal4 (bg)
      D4_3D,

      // length (ab)
      L_xn_2D,
      L_x0_2D,
      L_xN_2D,
      // height (ac)
      H_ny_2D,
      H_0y_2D,
      H_Ny_2D,
      // diagonal1 (bc)
      D1_2D,

      FIRST_EDGE_1D,
      LAST_EDGE_1D,
      CENTER_1D,
    };

    enum class TrianglePosition : char {
      //    e--------f
      //   /|       /|
      //  / |      / |
      // a--------b  |
      // |  g-----|--h
      // | /      | /
      // |/       |/
      // c--------d

      F_3D, // face (abc, bcd)
      C_3D, // side (abe, bef)
      H_3D, // top (acg, aeg)
      D1_3D, // diagonal1 (bdg, beg)
      D2_3D, // diagonal2 (abg, bgh)
      D3_3D, // diagonal3 (bcg, bfg)

      TOP_2D, // abc
      BOTTOM_2D, // bcd
    };

    bool hasPreconditionedVerticesAndCells_{false};

    float origin_[3]; //
    float spacing_[3]; //
    SimplexId nbvoxels_[3]; // nombre de voxels par axe

    // Vertex helper //
    SimplexId vshift_[2]; // VertexShift

    // Edge helper //
    SimplexId esetdims_[7]; // EdgeSetDimensions
    SimplexId esetshift_[7]; // EdgeSetShift
    SimplexId eshift_[14]; // EdgeShift

    // Triangle helper //
    SimplexId tsetdims_[6]; // TriangleSetDimensions
    SimplexId tsetshift_[6]; // TriangleSetShift
    SimplexId tshift_[12]; // TriangleShift

    // Tetrahedron helper //
    SimplexId tetshift_[2]; // TetrahedronShift

    SimplexId cellNumber_; // number of cells
    SimplexId vertexNumber_; // number of vertices
    SimplexId edgeNumber_; // number of edges
    SimplexId triangleNumber_; // number of triangles
    SimplexId tetrahedronNumber_; // number of tetrahedra

    // 2d helpers
    SimplexId Di_;
    SimplexId Dj_;

    // acceleration variables
    bool isAccelerated_;
    SimplexId mod_[2];
    SimplexId div_[2];

    // acceleration functions
    int checkAcceleration();
    bool isPowerOfTwo(unsigned long long int v, unsigned long long int &r);

    //\cond
    // 2D //
    void vertexToPosition2d(const SimplexId vertex,
                            SimplexId p[2]) const override;
    void
      edgeToPosition2d(const SimplexId edge, const int k, SimplexId p[2]) const;
    void triangleToPosition2d(const SimplexId triangle,
                              SimplexId p[2]) const override;

    SimplexId getVertexEdge2dA(const SimplexId p[2], const int id) const;
    SimplexId getVertexEdge2dB(const SimplexId p[2], const int id) const;
    SimplexId getVertexEdge2dC(const SimplexId p[2], const int id) const;
    SimplexId getVertexEdge2dD(const SimplexId p[2], const int id) const;
    SimplexId getVertexEdge2dAB(const SimplexId p[2], const int id) const;
    SimplexId getVertexEdge2dCD(const SimplexId p[2], const int id) const;
    SimplexId getVertexEdge2dAC(const SimplexId p[2], const int id) const;
    SimplexId getVertexEdge2dBD(const SimplexId p[2], const int id) const;
    SimplexId getVertexEdge2dABCD(const SimplexId p[2], const int id) const;

    SimplexId getVertexStar2dA(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2dB(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2dC(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2dD(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2dAB(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2dCD(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2dAC(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2dBD(const SimplexId p[2], const int id) const;
    SimplexId getVertexStar2dABCD(const SimplexId p[2], const int id) const;

    SimplexId getVertexLink2dA(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2dB(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2dC(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2dD(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2dAB(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2dCD(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2dAC(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2dBD(const SimplexId p[2], const int id) const;
    SimplexId getVertexLink2dABCD(const SimplexId p[2], const int id) const;

    SimplexId getEdgeTriangleL_x0(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_xn(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_xN(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_0y(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_ny(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_Ny(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleD1_xy(const SimplexId p[3], const int id) const;

    SimplexId getEdgeLink2dL(const SimplexId p[2], const int id) const;
    SimplexId getEdgeLink2dH(const SimplexId p[2], const int id) const;
    SimplexId getEdgeLink2dD1(const SimplexId p[2], const int id) const;

    SimplexId getEdgeStar2dL(const SimplexId p[2], const int id) const;
    SimplexId getEdgeStar2dH(const SimplexId p[2], const int id) const;

    // 3D //
    void vertexToPosition(const SimplexId vertex,
                          SimplexId p[3]) const override;
    void
      edgeToPosition(const SimplexId edge, const int k, SimplexId p[3]) const;
    void triangleToPosition(const SimplexId triangle,
                            const int k,
                            SimplexId p[3]) const override;
    void tetrahedronToPosition(const SimplexId tetrahedron,
                               SimplexId p[3]) const override;

    SimplexId getVertexEdgeA(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeB(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeC(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeD(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeE(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeF(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeG(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeH(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeAB(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeCD(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeEF(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeGH(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeAC(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeBD(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeEG(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeFH(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeAE(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeBF(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeCG(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeDH(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeABDC(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeEFHG(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeAEGC(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeBFHD(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeAEFB(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeGHDC(const SimplexId p[3], const int id) const;
    SimplexId getVertexEdgeABCDEFGH(const SimplexId p[3], const int id) const;

    SimplexId getVertexTriangleA(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleB(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleC(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleD(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleE(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleF(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleG(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleH(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleAB(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleCD(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleEF(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleGH(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleAC(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleBD(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleEG(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleFH(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleAE(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleBF(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleCG(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleDH(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleABDC(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleEFHG(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleAEGC(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleBFHD(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleAEFB(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleGHDC(const SimplexId p[3], const int id) const;
    SimplexId getVertexTriangleABCDEFGH(const SimplexId p[3],
                                        const int id) const;

    SimplexId getVertexLinkA(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkB(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkC(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkD(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkE(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkF(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkG(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkH(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkAB(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkCD(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkEF(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkGH(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkAC(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkBD(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkEG(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkFH(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkAE(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkBF(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkCG(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkDH(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkABDC(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkEFHG(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkAEGC(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkBFHD(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkAEFB(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkGHDC(const SimplexId p[3], const int id) const;
    SimplexId getVertexLinkABCDEFGH(const SimplexId p[3], const int id) const;

    SimplexId getVertexStarA(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarB(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarC(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarD(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarE(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarF(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarG(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarH(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarAB(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarCD(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarEF(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarGH(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarAC(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarBD(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarEG(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarFH(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarAE(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarBF(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarCG(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarDH(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarABDC(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarEFHG(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarAEGC(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarBFHD(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarAEFB(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarGHDC(const SimplexId p[3], const int id) const;
    SimplexId getVertexStarABCDEFGH(const SimplexId p[3], const int id) const;

    SimplexId getEdgeTriangleL_x00(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_x0n(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_x0N(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_xn0(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_xnn(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_xnN(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_xN0(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_xNn(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleL_xNN(const SimplexId p[3], const int id) const;

    SimplexId getEdgeTriangleH_0y0(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_0yn(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_0yN(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_ny0(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_nyn(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_nyN(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_Ny0(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_Nyn(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleH_NyN(const SimplexId p[3], const int id) const;

    SimplexId getEdgeTriangleP_00z(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleP_0nz(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleP_0Nz(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleP_n0z(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleP_nnz(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleP_nNz(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleP_N0z(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleP_Nnz(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleP_NNz(const SimplexId p[3], const int id) const;

    SimplexId getEdgeTriangleD1_xy0(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleD1_xyn(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleD1_xyN(const SimplexId p[3], const int id) const;

    SimplexId getEdgeTriangleD2_0yz(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleD2_nyz(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleD2_Nyz(const SimplexId p[3], const int id) const;

    SimplexId getEdgeTriangleD3_x0z(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleD3_xnz(const SimplexId p[3], const int id) const;
    SimplexId getEdgeTriangleD3_xNz(const SimplexId p[3], const int id) const;

    SimplexId getEdgeTriangleD4_xyz(const SimplexId p[3], const int id) const;

    SimplexId getEdgeLinkL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkP(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD1(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD2(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD3(const SimplexId p[3], const int id) const;
    SimplexId getEdgeLinkD4(const SimplexId p[3], const int id) const;

    SimplexId getEdgeStarL(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarH(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarP(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarD1(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarD2(const SimplexId p[3], const int id) const;
    SimplexId getEdgeStarD3(const SimplexId p[3], const int id) const;

    SimplexId getTriangleVertexF(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexH(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexC(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexD1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexD2(const SimplexId p[3], const int id) const;
    SimplexId getTriangleVertexD3(const SimplexId p[3], const int id) const;

    SimplexId getTriangleEdgeF_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeF_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeH_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeH_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeC_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeC_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD1_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD1_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD2_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD2_1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD3_0(const SimplexId p[3], const int id) const;
    SimplexId getTriangleEdgeD3_1(const SimplexId p[3], const int id) const;

    SimplexId getTriangleLinkF(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkH(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkC(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkD1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkD2(const SimplexId p[3], const int id) const;
    SimplexId getTriangleLinkD3(const SimplexId p[3], const int id) const;

    SimplexId getTriangleStarF(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarH(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarC(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarD1(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarD2(const SimplexId p[3], const int id) const;
    SimplexId getTriangleStarD3(const SimplexId p[3], const int id) const;

    SimplexId getTetrahedronVertexABCG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBCDG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexABEG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBEFG(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBFGH(const SimplexId p[3],
                                       const int id) const;
    SimplexId getTetrahedronVertexBDGH(const SimplexId p[3],
                                       const int id) const;

    SimplexId getTetrahedronEdgeABCG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBCDG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeABEG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBEFG(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBFGH(const SimplexId p[3], const int id) const;
    SimplexId getTetrahedronEdgeBDGH(const SimplexId p[3], const int id) const;

    SimplexId getTetrahedronTriangleABCG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBCDG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleABEG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBEFG(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBFGH(const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronTriangleBDGH(const SimplexId p[3],
                                         const int id) const;

    SimplexId getTetrahedronNeighborABCG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBCDG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborABEG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBEFG(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBFGH(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    SimplexId getTetrahedronNeighborBDGH(const SimplexId t,
                                         const SimplexId p[3],
                                         const int id) const;
    //\endcond

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

    const std::array<int,3> extent = {0,0,0};
    

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

  template <typename Derived>
  class ZombieTriangulationCRTP : public ZombieTriangulation {
    inline Derived &underlying() {
      return static_cast<Derived &>(*this);
    }
    inline Derived const &underlying() const {
      return static_cast<Derived const &>(*this);
    }

  public:
    inline SimplexId TTK_TRIANGULATION_INTERNAL(getVertexNeighborNumber)(
      const SimplexId &vertexId) const final {

      return lutNumNeighbor3d[getCaseID(vertexId)];

/*#ifndef TTK _ENABLE_KAMIKAZE
      if(vertexId < 0 or vertexId >= vertexNumber_)
        return -1;
#endif // !TTK_ENABLE_KAMIKAZE

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          return 14;
        case VertexPosition::FRONT_FACE_3D:
        case VertexPosition::BACK_FACE_3D:
        case VertexPosition::TOP_FACE_3D:
        case VertexPosition::BOTTOM_FACE_3D:
        case VertexPosition::LEFT_FACE_3D:
        case VertexPosition::RIGHT_FACE_3D:
          return 10;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          return 8;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
          return 7;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
        case VertexPosition::CENTER_2D:
          return 6;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
        case VertexPosition::TOP_EDGE_2D:
        case VertexPosition::BOTTOM_EDGE_2D:
        case VertexPosition::LEFT_EDGE_2D:
        case VertexPosition::RIGHT_EDGE_2D:
          return 4;
        case VertexPosition::TOP_RIGHT_CORNER_2D: // b
        case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
          return 3;
        case VertexPosition::TOP_LEFT_CORNER_2D: // a
        case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
        case VertexPosition::CENTER_1D:
          return 2;
        case VertexPosition::LEFT_CORNER_1D:
        case VertexPosition::RIGHT_CORNER_1D:
          return 1;
      }

      return -1; */
    }

    bool TTK_TRIANGULATION_INTERNAL(isVertexOnBoundary)(
      const SimplexId &vertexId) const override;

    bool TTK_TRIANGULATION_INTERNAL(isEdgeOnBoundary)(
      const SimplexId &edgeId) const override;

    inline int TTK_TRIANGULATION_INTERNAL(getVertexNeighbor)(
      const SimplexId &vertexId,
      const int &localNeighborId,
      SimplexId &neighborId) const final {

      neighborId = SimplexId(vertexId+lutNeighborIndexOffset[getCaseID(vertexId)][localNeighborId]);
      return 1;

/* 
#ifndef TTK_ENABLE_KAMIKAZE
      if(localNeighborId < 0
         or localNeighborId >= getVertexNeighborNumber(vertexId))
        return -1;
#endif // !TTK_ENABLE_KAMIKAZE

      switch(this->underlying().getVertexPosition(vertexId)) {
        case VertexPosition::CENTER_3D:
          neighborId
            = vertexId + this->vertexNeighborABCDEFGH_[localNeighborId];
          break;
        case VertexPosition::FRONT_FACE_3D:
          neighborId = vertexId + this->vertexNeighborABCD_[localNeighborId];
          break;
        case VertexPosition::BACK_FACE_3D:
          neighborId = vertexId + this->vertexNeighborEFGH_[localNeighborId];
          break;
        case VertexPosition::TOP_FACE_3D:
          neighborId = vertexId + this->vertexNeighborAEFB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_FACE_3D:
          neighborId = vertexId + this->vertexNeighborGHDC_[localNeighborId];
          break;
        case VertexPosition::LEFT_FACE_3D:
          neighborId = vertexId + this->vertexNeighborAEGC_[localNeighborId];
          break;
        case VertexPosition::RIGHT_FACE_3D:
          neighborId = vertexId + this->vertexNeighborBFHD_[localNeighborId];
          break;
        case VertexPosition::TOP_FRONT_EDGE_3D: // ab
          neighborId = vertexId + this->vertexNeighborAB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_FRONT_EDGE_3D: // cd
          neighborId = vertexId + this->vertexNeighborCD_[localNeighborId];
          break;
        case VertexPosition::LEFT_FRONT_EDGE_3D: // ac
          neighborId = vertexId + this->vertexNeighborAC_[localNeighborId];
          break;
        case VertexPosition::RIGHT_FRONT_EDGE_3D: // bd
          neighborId = vertexId + this->vertexNeighborBD_[localNeighborId];
          break;
        case VertexPosition::TOP_BACK_EDGE_3D: // ef
          neighborId = vertexId + this->vertexNeighborEF_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_BACK_EDGE_3D: // gh
          neighborId = vertexId + this->vertexNeighborGH_[localNeighborId];
          break;
        case VertexPosition::LEFT_BACK_EDGE_3D: // eg
          neighborId = vertexId + this->vertexNeighborEG_[localNeighborId];
          break;
        case VertexPosition::RIGHT_BACK_EDGE_3D: // fh
          neighborId = vertexId + this->vertexNeighborFH_[localNeighborId];
          break;
        case VertexPosition::TOP_LEFT_EDGE_3D: // ae
          neighborId = vertexId + this->vertexNeighborAE_[localNeighborId];
          break;
        case VertexPosition::TOP_RIGHT_EDGE_3D: // bf
          neighborId = vertexId + this->vertexNeighborBF_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_LEFT_EDGE_3D: // cg
          neighborId = vertexId + this->vertexNeighborCG_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_RIGHT_EDGE_3D: // dh
          neighborId = vertexId + this->vertexNeighborDH_[localNeighborId];
          break;
        case VertexPosition::TOP_LEFT_FRONT_CORNER_3D: // a
          neighborId = vertexId + this->vertexNeighborA_[localNeighborId];
          break;
        case VertexPosition::TOP_RIGHT_FRONT_CORNER_3D: // b
          neighborId = vertexId + this->vertexNeighborB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_LEFT_FRONT_CORNER_3D: // c
          neighborId = vertexId + this->vertexNeighborC_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_RIGHT_FRONT_CORNER_3D: // d
          neighborId = vertexId + this->vertexNeighborD_[localNeighborId];
          break;
        case VertexPosition::TOP_LEFT_BACK_CORNER_3D: // e
          neighborId = vertexId + this->vertexNeighborE_[localNeighborId];
          break;
        case VertexPosition::TOP_RIGHT_BACK_CORNER_3D: // f
          neighborId = vertexId + this->vertexNeighborF_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_LEFT_BACK_CORNER_3D: // g
          neighborId = vertexId + this->vertexNeighborG_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_RIGHT_BACK_CORNER_3D: // h
          neighborId = vertexId + this->vertexNeighborH_[localNeighborId];
          break;
        case VertexPosition::CENTER_2D:
          neighborId = vertexId + this->vertexNeighbor2dABCD_[localNeighborId];
          break;
        case VertexPosition::TOP_EDGE_2D:
          neighborId = vertexId + this->vertexNeighbor2dAB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_EDGE_2D:
          neighborId = vertexId + this->vertexNeighbor2dCD_[localNeighborId];
          break;
        case VertexPosition::LEFT_EDGE_2D:
          neighborId = vertexId + this->vertexNeighbor2dAC_[localNeighborId];
          break;
        case VertexPosition::RIGHT_EDGE_2D:
          neighborId = vertexId + this->vertexNeighbor2dBD_[localNeighborId];
          break;
        case VertexPosition::TOP_LEFT_CORNER_2D: // a
          neighborId = vertexId + this->vertexNeighbor2dA_[localNeighborId];
          break;
        case VertexPosition::TOP_RIGHT_CORNER_2D: // b
          neighborId = vertexId + this->vertexNeighbor2dB_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_LEFT_CORNER_2D: // c
          neighborId = vertexId + this->vertexNeighbor2dC_[localNeighborId];
          break;
        case VertexPosition::BOTTOM_RIGHT_CORNER_2D: // d
          neighborId = vertexId + this->vertexNeighbor2dD_[localNeighborId];
          break;
        case VertexPosition::CENTER_1D:
          neighborId = (localNeighborId == 0 ? vertexId + 1 : vertexId - 1);
          break;
        case VertexPosition::LEFT_CORNER_1D:
          neighborId = vertexId + 1;
          break;
        case VertexPosition::RIGHT_CORNER_1D:
          neighborId = vertexId - 1;
          break;
        default:
          neighborId = -1;
          break;
      }

      return 0; */
    }

    int getVertexEdgeInternal(const SimplexId &vertexId,
                              const int &localEdgeId,
                              SimplexId &edgeId) const override
    {
      
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

    }
                              
                              ;

    SimplexId
      getVertexTriangleNumberInternal(const SimplexId &vertexId) const override{
      return lutVertexTriangles3d[getCaseID(vertexId)];
    };

    int getVertexTriangleInternal(const SimplexId &vertexId,
                                  const int &localTriangleId,
                                  SimplexId &triangleId) const override {
      triangleId = SimplexId(-1);
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

      triangleId = SimplexId(res);
      return 1;
    };

    int TTK_TRIANGULATION_INTERNAL(getVertexLink)(
      const SimplexId &vertexId,
      const int &localLinkId,
      SimplexId &linkId) const override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getVertexStarNumber)(
      const SimplexId &vertexId) const override{
      return lutVertexTetrahedrons3d[getCaseID(vertexId)];
    };

    int TTK_TRIANGULATION_INTERNAL(getVertexStar)(
      const SimplexId &vertexId,
      const int &localStarId,
      SimplexId &starId) const override {
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

    int TTK_TRIANGULATION_INTERNAL(getVertexPoint)(const SimplexId &vertexId,
                                                   float &x,
                                                   float &y,
                                                   float &z) const override;

    int getEdgeVertexInternal(const SimplexId &edgeId,
                              const int &localVertexId,
                              SimplexId &vertexId) const override;

    SimplexId
      getEdgeTriangleNumberInternal(const SimplexId &edgeId) const override;

    int getEdgeTriangleInternal(const SimplexId &edgeId,
                                const int &id,
                                SimplexId &triangleId) const override;

    int
      TTK_TRIANGULATION_INTERNAL(getEdgeLink)(const SimplexId &edgeId,
                                              const int &localLinkId,
                                              SimplexId &linkId) const override;

    int
      TTK_TRIANGULATION_INTERNAL(getEdgeStar)(const SimplexId &edgeId,
                                              const int &localStarId,
                                              SimplexId &starId) const override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getEdgeStarNumber)(
      const SimplexId &edgeId) const override;

    int getTriangleVertexInternal(const SimplexId &triangleId,
                                  const int &localVertexId,
                                  SimplexId &vertexId) const override;

    int getTriangleEdgeInternal(const SimplexId &triangleId,
                                const int &id,
                                SimplexId &edgeId) const override;

    int TTK_TRIANGULATION_INTERNAL(getTriangleLink)(
      const SimplexId &triangleId,
      const int &localLinkId,
      SimplexId &linkId) const override;

    int TTK_TRIANGULATION_INTERNAL(getTriangleStar)(
      const SimplexId &triangleId,
      const int &localStarId,
      SimplexId &starId) const override;

    SimplexId TTK_TRIANGULATION_INTERNAL(getTriangleStarNumber)(
      const SimplexId &triangleId) const override;

    int getTriangleNeighbor(const SimplexId &triangleId,
                            const int &localNeighborId,
                            SimplexId &neighborId) const override;

    SimplexId
      getTriangleNeighborNumber(const SimplexId &triangleId) const override;

    int getTetrahedronVertex(const SimplexId &tetId,
                             const int &localVertexId,
                             SimplexId &vertexId) const override;

    int getTetrahedronEdge(const SimplexId &tetId,
                           const int &id,
                           SimplexId &edgeId) const override;

    int getTetrahedronTriangle(const SimplexId &tetId,
                               const int &id,
                               SimplexId &triangleId) const override;

    SimplexId
      getTetrahedronNeighborNumber(const SimplexId &tetId) const override;

    int getTetrahedronNeighbor(const SimplexId &tetId,
                               const int &localNeighborId,
                               SimplexId &neighborId) const override;
  };
} // namespace ttk

/// @cond

inline void
  ttk::ZombieTriangulation::vertexToPosition2d(const SimplexId vertex,
                                                 SimplexId p[2]) const {
  if(isAccelerated_) {
    p[0] = vertex & mod_[0];
    p[1] = vertex >> div_[0];
  } else {
    p[0] = vertex % vshift_[0];
    p[1] = vertex / vshift_[0];
  }
}

inline void ttk::ZombieTriangulation::edgeToPosition2d(const SimplexId edge,
                                                         const int k,
                                                         SimplexId p[2]) const {
  const int e = (k) ? edge - esetshift_[k - 1] : edge;
  p[0] = e % eshift_[2 * k];
  p[1] = e / eshift_[2 * k];
}

inline void
  ttk::ZombieTriangulation::triangleToPosition2d(const SimplexId triangle,
                                                   SimplexId p[2]) const {
  p[0] = triangle % tshift_[0];
  p[1] = triangle / tshift_[0];
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdge2dA(const SimplexId p[2],
                                               const int id) const {
  // V(a)={b,c}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // ac-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdge2dB(const SimplexId p[2],
                                               const int id) const {
  // V(b)={a,c,d}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] - 1; // ba-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // bd-H
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1; // bc-D1
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdge2dC(const SimplexId p[2],
                                               const int id) const {
  // V(c)={a,b,d}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // ca-H
    case 1:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0]; // cd-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdge2dD(const SimplexId p[2],
                                               const int id) const {
  // V(d)={c,b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] - 1; // dc-L
    case 1:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // db-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdge2dAB(const SimplexId p[2],
                                                const int id) const {
  // V(ab)=V(b)::{a,c,d}+V(a)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] - 1; // ba-L
    case 1:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // bd-H
    case 3:
      return p[0] + p[1] * eshift_[0]; // ab-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdge2dCD(const SimplexId p[2],
                                                const int id) const {
  // V(cd)=V(c)::{a,b,d}+V(d)::{c}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // ca-H
    case 1:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0]; // cd-L
    case 3:
      return p[0] + p[1] * eshift_[0] - 1; // dc-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdge2dAC(const SimplexId p[2],
                                                const int id) const {
  // V(ac)=V(c)::{a,b,d}+V(a)::{c}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // ca-H
    case 1:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0]; // cd-L
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // ac-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdge2dBD(const SimplexId p[2],
                                                const int id) const {
  // V(bd)=V(b)::{c,d}+V(d)::{b,c}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1; // bc-D1
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // bd-H
    case 2:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // db-H
    case 3:
      return p[0] + p[1] * eshift_[0] - 1; // dc-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdge2dABCD(const SimplexId p[2],
                                                  const int id) const {
  // V(abcd)=V(d)::{b,c}+V(c)::{b,d}+V(a)::{c}+V(b)::{c}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]; // db-H
    case 1:
      return p[0] + p[1] * eshift_[0] - 1; // dc-L
    case 2:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]; // cb-D1
    case 3:
      return p[0] + p[1] * eshift_[0]; // cd-L
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]; // ac-H
    case 5:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] - 1; // bc-D1
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStar2dA(const SimplexId p[2],
                                               const int /*id*/) const {
  return p[0] * 2 + p[1] * tshift_[0];
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStar2dB(const SimplexId p[2],
                                               const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStar2dC(const SimplexId p[2],
                                               const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 1:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStar2dD(const SimplexId p[2],
                                               const int /*id*/) const {
  return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStar2dAB(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStar2dCD(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 1:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 2:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStar2dAC(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 1:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStar2dBD(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 2:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStar2dABCD(const SimplexId p[2],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0];
    case 3:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 5:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLink2dA(const SimplexId p[2],
                                               const int /*id*/) const {
  return esetshift_[1] + p[0] + p[1] * eshift_[4]; // D1::bc
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLink2dB(const SimplexId p[2],
                                               const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1; // H::ac
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1; // L::ab
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLink2dC(const SimplexId p[2],
                                               const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1; // H::ac
    case 1:
      return p[0] + (p[1] - 1) * eshift_[0]; // L::ab
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLink2dD(const SimplexId p[2],
                                               const int /*id*/) const {
  return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1; // D1::bc
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLink2dAB(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1; // H::ac
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1; // L::ab
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLink2dCD(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1; // H::ac
    case 1:
      return p[0] + (p[1] - 1) * eshift_[0]; // L::ab
    case 2:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLink2dAC(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1; // H::ac
    case 1:
      return p[0] + (p[1] - 1) * eshift_[0]; // L::ab
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLink2dBD(const SimplexId p[2],
                                                const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1; // H::ac
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1; // L::ab
    case 2:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLink2dABCD(const SimplexId p[2],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] - 1; // H::ac
    case 1:
      return p[0] + (p[1] + 1) * eshift_[0] - 1; // L::ab
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]; // D1::bc
    case 3:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2] + 1; // H::ac
    case 4:
      return p[0] + (p[1] - 1) * eshift_[0]; // L::ab
    case 5:
      return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4] - 1; // D1::bc
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_x0(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_xn(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return p[Di_] * 2 + (p[Dj_] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_xN(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2 + (p[Dj_] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_0y(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Dj_] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_ny(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return (p[Di_] - 1) * 2 + p[Dj_] * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_Ny(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[Di_] - 1) * 2 + p[Dj_] * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD1_xy(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[Di_] * 2 + p[Dj_] * tshift_[0];
    case 1:
      return p[Di_] * 2 + p[Dj_] * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLink2dL(const SimplexId p[2],
                                             const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[Dj_]) {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * vshift_[0];
      case 1:
        return p[0] + (p[1] - 1) * vshift_[0] + 1;
    }
  } else if(p[1] == 0)
    return p[0] + vshift_[0];
  else
    return p[0] + (p[1] - 1) * vshift_[0] + 1;
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLink2dH(const SimplexId p[2],
                                             const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[Di_]) {
    switch(id) {
      case 0:
        return p[0] + p[1] * vshift_[0] + 1;
      case 1:
        return p[0] + (p[1] + 1) * vshift_[0] - 1;
    }
  } else if(p[0] == 0)
    return p[1] * vshift_[0] + 1;
  else
    return p[0] + (p[1] + 1) * vshift_[0] - 1;
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLink2dD1(const SimplexId p[2],
                                              const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0];
    case 1:
      return p[0] + (p[1] + 1) * vshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeStar2dL(const SimplexId p[2],
                                             const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[Dj_]) {
    if(id == 0)
      return p[0] * 2 + p[1] * tshift_[0];
    else
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
  } else if(p[1] == 0)
    return p[0] * 2 + p[1] * tshift_[0];
  else
    return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeStar2dH(const SimplexId p[2],
                                             const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[Di_]) {
    if(id == 0)
      return p[0] * 2 + p[1] * tshift_[0];
    else
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
  } else if(p[0] == 0)
    return p[0] * 2 + p[1] * tshift_[0];
  else
    return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
}

inline void ttk::ZombieTriangulation::vertexToPosition(const SimplexId vertex,
                                                         SimplexId p[3]) const {
  if(isAccelerated_) {
    p[0] = vertex & mod_[0];
    p[1] = (vertex & mod_[1]) >> div_[0];
    p[2] = vertex >> div_[1];
  } else {
    p[0] = vertex % vshift_[0];
    p[1] = (vertex % vshift_[1]) / vshift_[0];
    p[2] = vertex / vshift_[1];
  }
}

inline void ttk::ZombieTriangulation::edgeToPosition(const SimplexId edge,
                                                       const int k,
                                                       SimplexId p[3]) const {
  const int e = (k) ? edge - esetshift_[k - 1] : edge;
  p[0] = e % eshift_[2 * k];
  p[1] = (e % eshift_[2 * k + 1]) / eshift_[2 * k];
  p[2] = e / eshift_[2 * k + 1];
}

inline void ttk::ZombieTriangulation::triangleToPosition(
  const SimplexId triangle, const int k, SimplexId p[3]) const {
  const SimplexId t = (k) ? triangle - tsetshift_[k - 1] : triangle;
  p[0] = t % tshift_[2 * k];
  p[1] = (t % tshift_[2 * k + 1]) / tshift_[2 * k];
  p[2] = t / tshift_[2 * k + 1];
}

inline void
  ttk::ZombieTriangulation::tetrahedronToPosition(const SimplexId tetrahedron,
                                                    SimplexId p[3]) const {
  p[0] = (tetrahedron % tetshift_[0]) / 6;
  p[1] = (tetrahedron % tetshift_[1]) / tetshift_[0];
  p[2] = tetrahedron / tetshift_[1];
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeA(const SimplexId p[3],
                                             const int id) const {
  // V(a)={b,c,e,g}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeB(const SimplexId p[3],
                                             const int id) const {
  // V(b)={a,c,d,e,f,g,h}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeC(const SimplexId p[3],
                                             const int id) const {
  // V(c)={a,b,d,g}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ca-H
    case 1:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // cd-L
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // cg-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeD(const SimplexId p[3],
                                             const int id) const {
  // V(d)={b,c,g,h}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // dc-L
    case 2:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeE(const SimplexId p[3],
                                             const int id) const {
  // V(e)={a,b,f,g}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // ea-P
    case 1:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // eb-D3
    case 2:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ef-L
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // eg-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeF(const SimplexId p[3],
                                             const int id) const {
  // V(f)={b,e,g,h}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // fe-L
    case 2:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // fg-D1
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // fh-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeG(const SimplexId p[3],
                                             const int id) const {
  // V(g)={a,b,c,d,e,f,h}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeH(const SimplexId p[3],
                                             const int id) const {
  // V(h)={b,d,f,g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // hb-D2
    case 1:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // hd-P
    case 2:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // hf-H
    case 3:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeAB(const SimplexId p[3],
                                              const int id) const {
  // V(ab)=V(b)+V(a)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[11] + p[2] * eshift_[12]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeCD(const SimplexId p[3],
                                              const int id) const {
  // V(cd)=V(d)+V(c)::{b,d}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // dc-L
    case 2:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
    case 4:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // cb-D1
    case 5:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // cd-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeEF(const SimplexId p[3],
                                              const int id) const {
  // V(fe)=V(f)+V(e)::{b,f}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // fe-L
    case 2:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // fg-D1
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // fh-H
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // eb-D3
    case 5:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ef-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeGH(const SimplexId p[3],
                                              const int id) const {
  // V(gh)=V(g)+V(h)::{g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeAC(const SimplexId p[3],
                                              const int id) const {
  // V(ac)=V(c)+V(a)::{c,g}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ca-H
    case 1:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // cb-D1
    case 2:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // cd-L
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // cg-P
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 5:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeBD(const SimplexId p[3],
                                              const int id) const {
  // V(bd)=V(b)+V(d)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeEG(const SimplexId p[3],
                                              const int id) const {
  // V(eg)=V(g)+V(e)::{g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // eg-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeFH(const SimplexId p[3],
                                              const int id) const {
  // V(fh)=V(f)+V(h)::{b,f}
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // fe-L
    case 2:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // fg-D1
    case 3:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // fh-H
    case 4:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // hb-D2
    case 5:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // hf-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeAE(const SimplexId p[3],
                                              const int id) const {
  // V(ae)=V(a)+V(e)::{a,b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // ea-P
    case 5:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // eb-D3
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeBF(const SimplexId p[3],
                                              const int id) const {
  // V(bf)=V(b)+V(f)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeCG(const SimplexId p[3],
                                              const int id) const {
  // V(cg)=V(g)+V(c)::{g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // cg-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeDH(const SimplexId p[3],
                                              const int id) const {
  // V(dh)=V(d)+V(h)::{b,d}
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
    case 1:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // dc-L
    case 2:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 3:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
    case 4:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // hb-D2
    case 5:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // hd-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeABDC(const SimplexId p[3],
                                                const int id) const {
  // V(abdc)=V(b)+V(d)::{b}+V(c)::{b}+V(a)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
    case 8:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // cb-D1
    case 9:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeEFHG(const SimplexId p[3],
                                                const int id) const {
  // V(efhg)=V(g)+V(h)::{g}+V(f)::{g,h}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
    case 8:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // fg-D1
    case 9:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // fh-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeAEGC(const SimplexId p[3],
                                                const int id) const {
  // V(aegc)=V(g)+V(a)::{c,g}+V(c)::{g}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 8:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 9:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // cg-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeBFHD(const SimplexId p[3],
                                                const int id) const {
  // V(bfhd)=V(b)+V(f)::{b}+V(h)::{b}+V(d)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
    case 8:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // hb-D2
    case 9:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // db-H
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeAEFB(const SimplexId p[3],
                                                const int id) const {
  // V(aefb)=V(b)+V(a)::{b}+V(e)::{b}+V(f)::{b}
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // ba-L
    case 1:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 2:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // be-D3
    case 4:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // bf-P
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 6:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 8:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // eb-D3
    case 9:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // fb-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeGHDC(const SimplexId p[3],
                                                const int id) const {
  // V(ghdc)=V(g)+V(h)::{g}+V(d)::{g,h}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
    case 8:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 9:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexEdgeABCDEFGH(const SimplexId p[3],
                                                    const int id) const {
  // V(abcdefgh)=V(g)+V(d)::{g,h}+V(h)::{g}+V(b)::{c,d,g,h}
  switch(id) {
    case 0:
      return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
             + (p[2] - 1) * eshift_[9]; // ga-D2
    case 1:
      return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
             + (p[2] - 1) * eshift_[13]; // gb-D4
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + (p[2] - 1) * eshift_[5]; // gc-P
    case 3:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + (p[2] - 1) * eshift_[11]; // gd-D3
    case 4:
      return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
             + p[2] * eshift_[3]; // ge-H
    case 5:
      return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
             + p[2] * eshift_[7]; // gf-D1
    case 6:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // gh-L
    case 7:
      return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11]
             - 1; // dg-D3
    case 8:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // dh-P
    case 9:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1] - 1; // hg-L
    case 10:
      return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7]
             - 1; // bc-D1
    case 11:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // bd-H
    case 12:
      return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13]
             - 1; // bg-D4
    case 13:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // bh-D2
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleA(const SimplexId /*p*/[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return 0;
    case 1:
      return tsetshift_[0];
    case 2:
      return tsetshift_[1];
    case 3:
      return tsetshift_[3];
    case 4:
      return tsetshift_[1] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleB(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2;
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2;
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2;
    case 5:
      return tsetshift_[1] + p[0] * 2 + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2;
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2;
    case 11:
      return (p[0] - 1) * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleC(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return (p[1] - 1) * tshift_[0];
    case 1:
      return (p[1] - 1) * tshift_[0] + 1;
    case 2:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10];
    case 3:
      return tsetshift_[0] + p[1] * tshift_[2];
    case 4:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleD(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6];
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4];
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleE(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[2] - 1) * tshift_[3];
    case 1:
      return tsetshift_[2] + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return tsetshift_[1] + (p[2] - 1) * tshift_[5] + 1;
    case 3:
      return p[2] * tshift_[1];
    case 4:
      return tsetshift_[0] + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleF(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[3] + 1;
    case 1:
      return (p[0] - 1) * 2 + p[2] * tshift_[1];
    case 2:
      return (p[0] - 1) * 2 + p[2] * tshift_[1] + 1;
    case 3:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[11] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + (p[2] - 1) * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleG(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7];
    case 1:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5]
             + 1;
    case 4:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9];
    case 5:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9]
             + 1;
    case 6:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + (p[2] - 1) * tshift_[11]
             + 1;
    case 8:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3];
    case 11:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleH(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 1:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleAB(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2;
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2;
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2;
    case 5:
      return tsetshift_[1] + p[0] * 2 + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2;
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2;
    case 11:
      return (p[0] - 1) * 2;
    case 12:
      return p[0] * 2;
    case 13:
      return tsetshift_[0] + p[0] * 2;
    case 14:
      return tsetshift_[3] + p[0] * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleCD(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6];
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4];
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2] + 1;
    case 5:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 6:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 7:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10];
    case 8:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleEF(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[3] + 1;
    case 1:
      return (p[0] - 1) * 2 + p[2] * tshift_[1];
    case 2:
      return (p[0] - 1) * 2 + p[2] * tshift_[1] + 1;
    case 3:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[11] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + (p[2] - 1) * tshift_[5] + 1;
    case 5:
      return p[0] * 2 + tsetshift_[0] + (p[2] - 1) * tshift_[3];
    case 6:
      return p[0] * 2 + tsetshift_[2] + (p[2] - 1) * tshift_[7] + 1;
    case 7:
      return p[0] * 2 + p[2] * tshift_[1];
    case 8:
      return p[0] * 2 + tsetshift_[0] + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleGH(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 1:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 5:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7];
    case 6:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 7:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9];
    case 8:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 9:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 10:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 11:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 12:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 13:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3];
    case 14:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleAC(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[1] - 1) * tshift_[0];
    case 1:
      return (p[1] - 1) * tshift_[0] + 1;
    case 2:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10];
    case 3:
      return tsetshift_[0] + p[1] * tshift_[2];
    case 4:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4];
    case 5:
      return p[1] * tshift_[0];
    case 6:
      return tsetshift_[1] + p[1] * tshift_[4];
    case 7:
      return tsetshift_[3] + p[1] * tshift_[8];
    case 8:
      return tsetshift_[1] + p[1] * tshift_[4] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleBD(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6];
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4];
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2] + 1;
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10];
    case 7:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6];
    case 8:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8] + 1;
    case 9:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4];
    case 10:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + 1;
    case 11:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10] + 1;
    case 12:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6] + 1;
    case 13:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8];
    case 14:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleEG(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7];
    case 1:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5]
             + 1;
    case 4:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9];
    case 5:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9]
             + 1;
    case 6:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + (p[2] - 1) * tshift_[11]
             + 1;
    case 8:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3];
    case 11:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[2] + p[1] * tshift_[6] + (p[2] - 1) * tshift_[7] + 1;
    case 13:
      return tsetshift_[1] + p[1] * tshift_[4] + (p[2] - 1) * tshift_[5] + 1;
    case 14:
      return p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleFH(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 1:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 7:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 8:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleAE(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[2] - 1) * tshift_[3];
    case 1:
      return tsetshift_[2] + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return tsetshift_[1] + (p[2] - 1) * tshift_[5] + 1;
    case 3:
      return p[2] * tshift_[1];
    case 4:
      return tsetshift_[0] + (p[2] - 1) * tshift_[3] + 1;
    case 5:
      return tsetshift_[0] + p[2] * tshift_[3];
    case 6:
      return tsetshift_[1] + p[2] * tshift_[5];
    case 7:
      return tsetshift_[3] + p[2] * tshift_[9];
    case 8:
      return tsetshift_[1] + p[2] * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleBF(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[3] + 1;
    case 1:
      return (p[0] - 1) * 2 + p[2] * tshift_[1];
    case 2:
      return (p[0] - 1) * 2 + p[2] * tshift_[1] + 1;
    case 3:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[11] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + (p[2] - 1) * tshift_[5] + 1;
    case 5:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11];
    case 6:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[2] * tshift_[7];
    case 7:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[2] * tshift_[9] + 1;
    case 8:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5];
    case 9:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5] + 1;
    case 10:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11] + 1;
    case 11:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3] + 1;
    case 12:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[2] * tshift_[7] + 1;
    case 13:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[2] * tshift_[9];
    case 14:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleCG(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7];
    case 1:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5]
             + 1;
    case 4:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9];
    case 5:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9]
             + 1;
    case 6:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + (p[2] - 1) * tshift_[11]
             + 1;
    case 8:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3];
    case 11:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + p[2] * tshift_[11];
    case 13:
      return tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 14:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleDH(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 1:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 5:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
    case 6:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3];
    case 7:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 8:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleABDC(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10];
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10] + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2] + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6] + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8];
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2];
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0];
    case 12:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 13:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6];
    case 14:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4];
    case 15:
      return p[0] * 2 + (p[1] - 1) * tshift_[0];
    case 16:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 17:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10];
    case 18:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2];
    case 19:
      return p[0] * 2 + p[1] * tshift_[0];
    case 20:
      return p[0] * 2 + tsetshift_[3] + p[1] * tshift_[8];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleEFHG(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7];
    case 1:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9];
    case 5:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 6:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 8:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3];
    case 11:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 13:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 14:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 15:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 16:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 17:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 18:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 19:
      return p[0] * 2 + tsetshift_[2] + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 20:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleAEGC(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7];
    case 1:
      return tsetshift_[2] + (p[1] - 1) * tshift_[6] + (p[2] - 1) * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5];
    case 3:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + (p[2] - 1) * tshift_[5]
             + 1;
    case 4:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9];
    case 5:
      return tsetshift_[3] + (p[1] - 1) * tshift_[8] + (p[2] - 1) * tshift_[9]
             + 1;
    case 6:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + (p[2] - 1) * tshift_[11]
             + 1;
    case 8:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3];
    case 11:
      return tsetshift_[0] + p[1] * tshift_[2] + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[2] + p[1] * tshift_[6] + (p[2] - 1) * tshift_[7] + 1;
    case 13:
      return tsetshift_[1] + p[1] * tshift_[4] + (p[2] - 1) * tshift_[5] + 1;
    case 14:
      return p[1] * tshift_[0] + p[2] * tshift_[1];
    case 15:
      return tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 16:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 17:
      return tsetshift_[3] + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 18:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5] + 1;
    case 19:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + p[2] * tshift_[11];
    case 20:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleBFHD(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9];
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3];
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 12:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 13:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
    case 14:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 15:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 16:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 17:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 18:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 19:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 20:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleAEFB(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[2] * tshift_[7];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[2] * tshift_[9] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5] + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11] + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3] + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[2] * tshift_[7] + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[2] * tshift_[9];
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3];
    case 11:
      return (p[0] - 1) * 2 + p[2] * tshift_[1];
    case 12:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[3] + 1;
    case 13:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[2] - 1) * tshift_[11] + 1;
    case 14:
      return tsetshift_[1] + p[0] * 2 + (p[2] - 1) * tshift_[5] + 1;
    case 15:
      return p[0] * 2 + tsetshift_[0] + (p[2] - 1) * tshift_[3];
    case 16:
      return p[0] * 2 + tsetshift_[2] + (p[2] - 1) * tshift_[7] + 1;
    case 17:
      return p[0] * 2 + p[2] * tshift_[1];
    case 18:
      return p[0] * 2 + tsetshift_[0] + (p[2] - 1) * tshift_[3] + 1;
    case 19:
      return p[0] * 2 + tsetshift_[0] + p[2] * tshift_[3];
    case 20:
      return p[0] * 2 + tsetshift_[3] + p[2] * tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleGHDC(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7];
    case 1:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 3:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9];
    case 5:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 6:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 7:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 8:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 9:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 10:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3];
    case 11:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 12:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 13:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 14:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 15:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
    case 16:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3];
    case 17:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 18:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 19:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11];
    case 20:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexTriangleABCDEFGH(const SimplexId p[3],
                                                        const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1;
    case 4:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 5:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1;
    case 7:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 8:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7] + 1;
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9];
    case 10:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3];
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 12:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7];
    case 13:
      return p[0] * 2 + tsetshift_[2] + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 14:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5];
    case 15:
      return p[0] * 2 + tsetshift_[1] + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 16:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9];
    case 17:
      return p[0] * 2 + tsetshift_[3] + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 18:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11];
    case 19:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 20:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1];
    case 21:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 22:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3];
    case 23:
      return p[0] * 2 + tsetshift_[0] + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 24:
      return (p[0] - 1) * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 25:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
    case 26:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 27:
      return p[0] * 2 + tsetshift_[2] + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 28:
      return p[0] * 2 + tsetshift_[1] + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 29:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 30:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 31:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 32:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 33:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 34:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
    case 35:
      return p[0] * 2 + tsetshift_[4] + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkA(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkB(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 4:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkC(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkD(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkE(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkF(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 1:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkG(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 5:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkH(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkAB(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 7:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkCD(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
    case 2:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkEF(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkGH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 5:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 7:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkAC(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 3:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkBD(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 4:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 7:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkEG(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 2:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 6:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 7:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkFH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 1:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 2:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 3:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkAE(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkBF(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 4:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 7:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkCG(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
    case 2:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 6:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 7:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkDH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
    case 2:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 3:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkABDC(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 7:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 8:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
    case 10:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 11:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkEFHG(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 2:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 4:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 5:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
    case 6:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 7:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 8:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 9:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 10:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 11:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkAEGC(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 4:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 5:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 6:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 7:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 8:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 9:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
    case 10:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 11:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkBFHD(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 4:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 5:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 7:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 8:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 9:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
    case 10:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 11:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkAEFB(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 4:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 5:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 6:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 7:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 8:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 9:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 10:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 11:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkGHDC(const SimplexId p[3],
                                                const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 1:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 2:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 4:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 5:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
    case 6:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 7:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
    case 8:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 9:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
    case 10:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 11:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexLinkABCDEFGH(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1; // D1::beg
    case 2:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5]; // C::acg
    case 3:
      return tsetshift_[1] + (p[0] - 1) * 2 + p[1] * tshift_[4]
             + p[2] * tshift_[5] + 1; // C::aeg
    case 4:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3]; // H::abe
    case 5:
      return tsetshift_[0] + (p[0] - 1) * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1; // H::bef
    case 6:
      return (p[0] - 1) * 2 + p[1] * tshift_[0]
             + (p[2] + 1) * tshift_[1]; // F::abc
    case 7:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1]
             + 1; // F::bcd
    case 8:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9]; // D2::abg
    case 9:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7]; // D1::bdg
    case 10:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11]; // D3::bcg
    case 11:
      return tsetshift_[3] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[8]
             + p[2] * tshift_[9] + 1; // D2::bgh
    case 12:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9]; // D2::abg
    case 13:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 14:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1; // D1::beg
    case 15:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1; // D2::bgh
    case 16:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5]; // C::acg
    case 17:
      return tsetshift_[1] + (p[0] + 1) * 2 + (p[1] - 1) * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1; // C::aeg
    case 18:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3]; // H::abe
    case 19:
      return tsetshift_[0] + p[0] * 2 + (p[1] - 1) * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1; // H::bef
    case 20:
      return p[0] * 2 + (p[1] - 1) * tshift_[0]
             + (p[2] - 1) * tshift_[1]; // F::abc
    case 21:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + (p[2] - 1) * tshift_[1]
             + 1; // F::bcd
    case 22:
      return tsetshift_[4] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1; // D3::bfg
    case 23:
      return tsetshift_[2] + (p[0] - 1) * 2 + (p[1] - 1) * tshift_[6]
             + (p[2] - 1) * tshift_[7]; // D1::bdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarA(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // abcg
    case 1:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2; // abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarB(const SimplexId p[3],
                                             const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + id;
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarC(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]; // abcg
    case 1:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // bcdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarD(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // bcdg
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarE(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // abeg
    case 1:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarF(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // befg
    case 1:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // bfgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarG(const SimplexId p[3],
                                             const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id;
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarH(const SimplexId p[3],
                                             const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // bfgh
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarAB(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 7:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarCD(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
    case 2:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 3:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarEF(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 1:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 2:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 3:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarGH(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarAC(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 1:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 2:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 3:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarBD(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarEG(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 7:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarFH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 2:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 3:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarAE(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 1:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
    case 2:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 3:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarBF(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 7:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarCG(const SimplexId p[3],
                                              const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 7:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarDH(const SimplexId p[3],
                                              const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 1:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 2:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 3:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarABDC(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
    case 8:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 9:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 10:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 11:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarEFHG(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 8:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 9:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 10:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 11:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarAEGC(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 7:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 8:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 9:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
    case 10:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 11:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarBFHD(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 7:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 8:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 9:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 10:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 11:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarAEFB(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  switch(id) {
    case 6:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 7:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
    case 8:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 9:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 10:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 11:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarGHDC(const SimplexId p[3],
                                                const int id) const {
  if(id >= 0 && id <= 5)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1]
           + id; // tet(g)
  switch(id) {
    case 6:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 7:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
    case 8:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 9:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 10:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 11:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getVertexStarABCDEFGH(const SimplexId p[3],
                                                    const int id) const {
  if(id >= 0 && id <= 5)
    return (p[0] - 1) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
           + id; // tet(b)
  if(id >= 6 && id <= 11)
    return p[0] * 6 + (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1] + id
           - 6; // tet(g)
  switch(id) {
    case 12:
      return p[0] * 6 + p[1] * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(a)::abcg
    case 13:
      return p[0] * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // tet(a)::abeg
    case 14:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0]
             + p[2] * tetshift_[1]; // tet(c)::abcg
    case 15:
      return p[0] * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(c)::bcdg
    case 16:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // tet(d)::bcdg
    case 17:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // tet(d)::bdgh
    case 18:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 2; // tet(e)::abeg
    case 19:
      return p[0] * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(e)::befg
    case 20:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // tet(f)::befg
    case 21:
      return (p[0] - 1) * 6 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // tet(f)::bfgh
    case 22:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 4; // tet(h)::bfgh
    case 23:
      return (p[0] - 1) * 6 + (p[1] - 1) * tetshift_[0]
             + (p[2] - 1) * tetshift_[1] + 5; // tet(h)::bdgh
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_x00(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2;
    case 1:
      return tsetshift_[0] + p[0] * 2;
    case 2:
      return tsetshift_[3] + p[0] * 2;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_x0n(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[2] * tshift_[1];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[2] * tshift_[3];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[2] * tshift_[9];
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_x0N(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[2] * tshift_[1];
    case 1:
      return tsetshift_[0] + p[0] * 2 + (p[2] - 1) * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_xn0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8];
    case 3:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_xnn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 2:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 3:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 4:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 5:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_xnN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 3:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_xN0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_xNn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleL_xNN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + (p[1] - 1) * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2]
             + (p[2] - 1) * tshift_[3] + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + (p[1] - 1) * tshift_[8]
             + (p[2] - 1) * tshift_[9] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_0y0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return p[1] * tshift_[0];
    case 1:
      return tsetshift_[1] + p[1] * tshift_[4];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_0yn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[1] * tshift_[4] + (p[2] - 1) * tshift_[5] + 1;
    case 1:
      return tsetshift_[2] + p[1] * tshift_[6] + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[1] * tshift_[0] + p[2] * tshift_[1];
    case 3:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_0yN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[1] * tshift_[4] + (p[2] - 1) * tshift_[5] + 1;
    case 1:
      return tsetshift_[2] + p[1] * tshift_[6] + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_ny0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4];
    case 3:
      return p[0] * 2 + p[1] * tshift_[0];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_nyn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 4:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 5:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_nyN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6]
             + (p[2] - 1) * tshift_[7] + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 3:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_Ny0(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_Nyn(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[2] + (p[0] - 1) * 2 + p[1] * tshift_[6]
             + p[2] * tshift_[7];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 3:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleH_NyN(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4]
             + (p[2] - 1) * tshift_[5] + 1;
    case 1:
      return (p[0] - 1) * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleP_00z(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[1] + p[2] * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleP_0nz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + p[2] * tshift_[11];
    case 2:
      return tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 3:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleP_0Nz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + (p[1] - 1) * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[4] + (p[1] - 1) * tshift_[10] + p[2] * tshift_[11];
    case 2:
      return tsetshift_[0] + p[1] * tshift_[2] + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleP_n0z(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5] + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleP_nnz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 4:
      return tsetshift_[4] + p[0] * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11];
    case 5:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleP_nNz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
    case 2:
      return tsetshift_[4] + p[0] * 2 + (p[1] - 1) * tshift_[10]
             + p[2] * tshift_[11];
    case 3:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleP_N0z(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[2] * tshift_[11] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[2] * tshift_[5] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleP_Nnz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[4] + (p[0] - 1) * 2 + p[1] * tshift_[10]
             + p[2] * tshift_[11] + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 3:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleP_NNz(const SimplexId p[3],
                                                   const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + (p[0] - 1) * 2 + p[1] * tshift_[2]
             + p[2] * tshift_[3] + 1;
    case 1:
      return tsetshift_[1] + p[0] * 2 + (p[1] - 1) * tshift_[4]
             + p[2] * tshift_[5];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD1_xy0(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0];
    case 1:
      return p[0] * 2 + p[1] * tshift_[0] + 1;
    case 2:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD1_xyn(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 1:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD1_xyN(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 1:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 2:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10]
             + (p[2] - 1) * tshift_[11] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD2_0yz(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[1] + p[1] * tshift_[4] + p[2] * tshift_[5] + 1;
    case 2:
      return tsetshift_[3] + p[1] * tshift_[8] + p[2] * tshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD2_nyz(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 3:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD2_Nyz(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 1:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 2:
      return tsetshift_[3] + (p[0] - 1) * 2 + p[1] * tshift_[8]
             + p[2] * tshift_[9] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD3_x0z(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[2] * tshift_[3];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[2] * tshift_[3] + 1;
    case 2:
      return tsetshift_[2] + p[0] * 2 + p[2] * tshift_[7] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD3_xnz(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3]
             + 1;
    case 2:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 3:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD3_xNz(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3]
             + 1;
    case 2:
      return tsetshift_[2] + p[0] * 2 + (p[1] - 1) * tshift_[6]
             + p[2] * tshift_[7];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeTriangleD4_xyz(const SimplexId p[3],
                                                    const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 3:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 4:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 5:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLinkL(const SimplexId p[3],
                                           const int id) const {
  if(p[2] == 0 and p[1] == 0) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + eshift_[4]; // CG
      case 1:
        return esetshift_[0] + p[0] + eshift_[3]; // EG
    }
  } else if(p[2] == 0 and p[1] == nbvoxels_[1])
    return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]; // BG
  else if(p[2] == nbvoxels_[2] and p[1] == 0)
    return esetshift_[5] + p[0] + (p[2] - 1) * eshift_[13]; // BG
  else if(p[2] == nbvoxels_[2] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]
               + (p[2] - 1) * eshift_[5] + 1; // BF
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + (p[2] - 1) * eshift_[3] + 1; // BD
    }
  } else if(p[2] == 0 and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4]; // CG
      case 1:
        return esetshift_[0] + p[0] + p[1] * eshift_[2] + eshift_[3]; // EG
      case 2:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]; // BG
    }
  } else if(p[2] == nbvoxels_[2] and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]
               + (p[2] - 1) * eshift_[5] + 1; // BF
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + (p[2] - 1) * eshift_[3] + 1; // BD
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
        // BG
    }
  } else if(p[2] > 0 and p[2] < nbvoxels_[2] and p[1] == 0) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + p[2] * eshift_[5] + eshift_[4]; // CG
      case 1:
        return esetshift_[0] + p[0] + (p[2] + 1) * eshift_[3]; // EG
      case 2:
        return esetshift_[5] + p[0] + (p[2] - 1) * eshift_[13]; // BG
    }
  } else if(p[2] > 0 and p[2] < nbvoxels_[2] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] - 1) * eshift_[4]
               + (p[2] - 1) * eshift_[5] + 1; // BF
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + (p[2] - 1) * eshift_[3] + 1; // BD
      case 2:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13];
        // BG
    }
  } else {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4]
               + p[2] * eshift_[5];
        // CG
      case 1:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3];
        // EG
      case 2:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13];
      case 3:
        return esetshift_[1] + p[0] + 1 + (p[1] - 1) * eshift_[4]
               + (p[2] - 1) * eshift_[5]; // CG
      case 4:
        return esetshift_[0] + p[0] + 1 + (p[1] - 1) * eshift_[2]
               + (p[2] - 1) * eshift_[3]; // EG
      case 5:
        return esetshift_[5] + p[0] + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLinkH(const SimplexId p[3],
                                           const int id) const {
  if(p[0] == 0 and p[2] == 0)
    return esetshift_[5] + p[1] * eshift_[12];
  else if(p[0] == nbvoxels_[0] and p[2] == 0) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4] - 1;
      case 1:
        return p[0] + (p[1] + 1) * eshift_[0] + eshift_[1] - 1;
    }
  } else if(p[0] == 0 and p[2] == nbvoxels_[2]) {
    switch(id) {
      case 0:
        return p[1] * eshift_[0] + (p[2] - 1) * eshift_[1];
      case 1:
        return esetshift_[1] + p[1] * eshift_[4] + (p[2] - 1) * eshift_[5] + 1;
    }
  } else if(p[0] == nbvoxels_[0] and p[2] == nbvoxels_[2])
    return esetshift_[5] + p[0] - 1 + p[1] * eshift_[12]
           + (p[2] - 1) * eshift_[13];
  else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[2] == 0) {
    switch(id) {
      case 0:
        return p[0] - 1 + (p[1] + 1) * eshift_[0] + eshift_[1];
      case 1:
        return esetshift_[1] + p[0] - 1 + (p[1] + 1) * eshift_[4];
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12];
    }
  } else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[2] == nbvoxels_[2]) {
    switch(id) {
      case 0:
        return esetshift_[5] + p[0] - 1 + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
      case 1:
        return p[0] + p[1] * eshift_[0] + (p[2] - 1) * eshift_[1];
      case 2:
        return esetshift_[1] + p[0] + 1 + p[1] * eshift_[4]
               + (p[2] - 1) * eshift_[5];
    }
  } else if(p[0] == 0 and p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return esetshift_[1] + p[1] * eshift_[4] + (p[2] - 1) * eshift_[5] + 1;
      case 1:
        return esetshift_[5] + p[1] * eshift_[12] + p[2] * eshift_[13];
      case 2:
        return p[1] * eshift_[0] + (p[2] - 1) * eshift_[1];
    }
  } else if(p[0] == nbvoxels_[0] and p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return esetshift_[5] + p[0] - 1 + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
      case 1:
        return esetshift_[1] + p[0] - 1 + (p[1] + 1) * eshift_[4]
               + p[2] * eshift_[5];
      case 2:
        return p[0] - 1 + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1];
    }
  } else {
    switch(id) {
      case 0:
        return p[0] + p[1] * eshift_[0] + (p[2] - 1) * eshift_[1];
      case 1:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
      case 2:
        return esetshift_[1] + p[0] - 1 + (p[1] + 1) * eshift_[4]
               + p[2] * eshift_[5];
      case 3:
        return esetshift_[1] + p[0] + 1 + p[1] * eshift_[4]
               + (p[2] - 1) * eshift_[5];
      case 4:
        return esetshift_[5] + (p[0] - 1) + p[1] * eshift_[12]
               + (p[2] - 1) * eshift_[13];
      case 5:
        return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLinkP(const SimplexId p[3],
                                           const int id) const {
  if(p[0] == 0 and p[1] == 0)
    return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
  else if(p[0] == 0 and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + p[2] * eshift_[3] + 1;
      case 1:
        return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1];
    }
  } else if(p[0] == nbvoxels_[0] and p[1] == 0) {
    switch(id) {
      case 0:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3] - 1;
      case 1:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
    }
  } else if(p[0] == nbvoxels_[0] and p[1] == nbvoxels_[1])
    return esetshift_[5] + p[0] - 1 + (p[1] - 1) * eshift_[12]
           + p[2] * eshift_[13];
  else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[1] == 0) {
    switch(id) {
      case 0:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
      case 1:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3] - 1;
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
    }
  } else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1];
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + p[2] * eshift_[3] + 1;
      case 2:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13] - 1;
    }
  } else if(p[0] == 0 and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1];
      case 1:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + p[2] * eshift_[3] + 1;
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
    }
  } else if(p[0] == nbvoxels_[0] and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13] - 1;
      case 1:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3] - 1;
      case 2:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] + (p[1] - 1) * eshift_[0] + p[2] * eshift_[1];
      case 1:
        return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1] - 1;
      case 2:
        return esetshift_[5] + p[0] + p[1] * eshift_[12] + p[2] * eshift_[13];
      case 3:
        return esetshift_[5] + p[0] + (p[1] - 1) * eshift_[12]
               + p[2] * eshift_[13] - 1;
      case 4:
        return esetshift_[0] + p[0] + p[1] * eshift_[2]
               + (p[2] + 1) * eshift_[3] - 1;
      case 5:
        return esetshift_[0] + p[0] + (p[1] - 1) * eshift_[2]
               + p[2] * eshift_[3] + 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLinkD1(const SimplexId p[3],
                                            const int id) const {
  if(p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return esetshift_[4] + p[0] + p[1] * eshift_[10]
               + (p[2] - 1) * eshift_[11];
      case 1:
        return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
               + p[2] * eshift_[11];
      case 2:
        return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
      case 3:
        return esetshift_[3] + p[0] + p[1] * eshift_[8]
               + (p[2] - 1) * eshift_[9] + 1;
    }
  } else if(p[2] == 0) {
    switch(id) {
      case 0:
        return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
      case 1:
        return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
               + p[2] * eshift_[11];
    }
  } else {
    switch(id) {
      case 0:
        return esetshift_[3] + p[0] + p[1] * eshift_[8]
               + (p[2] - 1) * eshift_[9] + 1;
      case 1:
        return esetshift_[4] + p[0] + p[1] * eshift_[10]
               + (p[2] - 1) * eshift_[11];
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLinkD2(const SimplexId p[3],
                                            const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[0]) {
    switch(id) {
      case 0:
        return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
      case 1:
        return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
               + p[2] * eshift_[11] - 1;
      case 2:
        return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
      case 3:
        return esetshift_[2] + p[0] + p[1] * eshift_[6]
               + (p[2] + 1) * eshift_[7] - 1;
    }
  } else if(p[0] == 0) {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + p[1] * eshift_[6] + p[2] * eshift_[7];
      case 1:
        return esetshift_[4] + p[0] + p[1] * eshift_[10] + p[2] * eshift_[11];
    }
  } else {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + p[1] * eshift_[6]
               + (p[2] + 1) * eshift_[7] - 1;
      case 1:
        return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
               + p[2] * eshift_[11] - 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLinkD3(const SimplexId p[3],
                                            const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
               + p[2] * eshift_[7];
      case 1:
        return esetshift_[2] + p[0] + p[1] * eshift_[6]
               + (p[2] + 1) * eshift_[7];
      case 2:
        return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
      case 3:
        return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
               + p[2] * eshift_[9] + 1;
    }
  } else if(p[1] == 0) {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + p[1] * eshift_[6]
               + (p[2] + 1) * eshift_[7];
      case 1:
        return esetshift_[3] + p[0] + p[1] * eshift_[8] + p[2] * eshift_[9];
    }
  } else {
    switch(id) {
      case 0:
        return esetshift_[2] + p[0] + (p[1] - 1) * eshift_[6]
               + p[2] * eshift_[7];
      case 1:
        return esetshift_[3] + p[0] + (p[1] - 1) * eshift_[8]
               + p[2] * eshift_[9] + 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeLinkD4(const SimplexId p[3],
                                            const int id) const {
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return p[0] + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1];
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 3:
      return esetshift_[1] + p[0] + 1 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5];
    case 4:
      return esetshift_[0] + p[0] + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 5:
      return esetshift_[0] + p[0] + 1 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeStarL(const SimplexId p[3],
                                           const int id) const {
  if(p[2] == 0 and p[1] == 0) {
    switch(id) {
      case 0:
        return p[0] * 6; // ABCG
      case 1:
        return p[0] * 6 + 2; // ABEG
    }
  } else if(p[2] == 0 and p[1] == nbvoxels_[1])
    return (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1; // BCDG
  else if(p[2] == nbvoxels_[2] and p[1] == 0)
    return (p[2] - 1) * tetshift_[1] + p[0] * 6 + 3; // BEFG
  else if(p[2] == nbvoxels_[2] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1] + p[0] * 6
               + 4;
        // BFGH
      case 1:
        return (p[1] - 1) * tetshift_[0] + (p[2] - 1) * tetshift_[1] + p[0] * 6
               + 5;
        // BDGH
    }
  } else if(p[2] == 0 and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 2:
        return (p[1] - 1) * tetshift_[0] + p[0] * 6 + 1; // BCDG
    }
  } else if(p[2] == nbvoxels_[2] and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 1:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 4;
        // BFGH
      case 2:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5;
        // BDGH
    }
  } else if(p[2] > 0 and p[2] < nbvoxels_[2] and p[1] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[0] * 6 + 2; // ABEG
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[0] * 6 + 3; // BEFG
    }
  } else if(p[2] > 0 and p[2] < nbvoxels_[2] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 1:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 4;
        // BFGH
      case 2:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5;
        // BDGH
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 3:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 4:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 4;
        // BFGH
      case 5:
        return (p[2] - 1) * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5;
        // BDGH
    }
  }

  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeStarH(const SimplexId p[3],
                                           const int id) const {
  if(p[0] == 0 and p[2] == 0)
    return p[1] * tetshift_[0]; // ABCG
  else if(p[0] == nbvoxels_[0] and p[2] == 0) {
    switch(id) {
      case 0:
        return p[1] * tetshift_[0] + (p[0] - 1) * 6 + 1; // BCDG
      case 1:
        return p[1] * tetshift_[0] + (p[0] - 1) * 6 + 5; // BDGH
    }
  } else if(p[0] == 0 and p[2] == nbvoxels_[2]) {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + 2; // ABEG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + 3; // BEFG
    }
  } else if(p[0] == nbvoxels_[0] and p[2] == nbvoxels_[2])
    return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
           + 4; // BFGH
  else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[2] == 0) {
    switch(id) {
      case 0:
        return p[1] * tetshift_[0] + (p[0] - 1) * 6 + 1; // BCDG
      case 1:
        return p[1] * tetshift_[0] + (p[0] - 1) * 6 + 5; // BDGH
      case 2:
        return p[1] * tetshift_[0] + p[0] * 6; // ABCG
    }
  } else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[2] == nbvoxels_[2]) {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 2; // ABEG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4;
        // BFGH
    }
  } else if(p[0] == 0 and p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + 2; // ABEG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + 3; // BEFG
      case 2:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0]; // ABCG
    }
  } else if(p[0] == nbvoxels_[0] and p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 1; // BCDG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 5; // BDGH
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4;
        // BFGH
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 2; // ABEG
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 3:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4;
        // BFGH
      case 4:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 1; // BCDG
      case 5:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 5; // BDGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeStarP(const SimplexId p[3],
                                           const int id) const {
  if(p[0] == 0 and p[1] == 0)
    return p[2] * tetshift_[1] + 2; // ABEG
  else if(p[0] == 0 and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0]; // ABCG
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + 1; // BCDG
    }
  } else if(p[0] == nbvoxels_[0] and p[1] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[0] - 1) * 6 + 3; // BEFG
      case 1:
        return p[2] * tetshift_[1] + (p[0] - 1) * 6 + 4; // BFGH
    }
  } else if(p[0] == nbvoxels_[0] and p[1] == nbvoxels_[1])
    return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
           + 5; // BDGH
  else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[1] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[0] - 1) * 6 + 3; // BEFG
      case 1:
        return p[2] * tetshift_[1] + (p[0] - 1) * 6 + 4; // BFGH
      case 2:
        return p[2] * tetshift_[1] + p[0] * 6 + 2; // ABEG
    }
  } else if(p[0] > 0 and p[0] < nbvoxels_[0] and p[1] == nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0]
               + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
               + 5;
        // BDGH
    }
  } else if(p[0] == 0 and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0]; // ABCG
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + 1; // BCDG
      case 2:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + 2; // ABEG
    }
  } else if(p[0] == nbvoxels_[0] and p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 3; // BEFG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4; // BFGH
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
               + 5;
        // BDGH
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + (p[0] - 1) * 6
               + 5;
        // BDGH
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0]
               + p[0] * 6; // ABCG
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 3:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 4:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 3; // BEFG
      case 5:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4; // BFGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeStarD1(const SimplexId p[3],
                                            const int id) const {
  if(p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 1; // BCDG
      case 2:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 3:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 4; // BFGH
    }
  } else if(p[2] == 0) {
    switch(id) {
      case 0:
        return p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[1] * tetshift_[0] + p[0] * 6 + 1; // BCDG
    }
  } else {
    switch(id) {
      case 0:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 3; // BEFG
      case 1:
        return (p[2] - 1) * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6
               + 4; // BFGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeStarD2(const SimplexId p[3],
                                            const int id) const {
  if(p[0] > 0 and p[0] < nbvoxels_[0]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 2:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 5; // BDGH
      case 3:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4; // BFGH
    }
  } else if(p[0] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0]; // ABCG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + 2; // ABEG
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 5; // BDGH
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + (p[0] - 1) * 6
               + 4; // BFGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getEdgeStarD3(const SimplexId p[3],
                                            const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3; // BEFG
      case 2:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 3:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5; // BDGH
    }
  } else if(p[1] == 0) {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 2; // ABEG
      case 1:
        return p[2] * tetshift_[1] + p[1] * tetshift_[0] + p[0] * 6 + 3; // BEFG
    }
  } else {
    switch(id) {
      case 0:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 1; // BCDG
      case 1:
        return p[2] * tetshift_[1] + (p[1] - 1) * tetshift_[0] + p[0] * 6
               + 5; // BDGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleVertexF(const SimplexId p[3],
                                                 const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1;
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
  }
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleVertexH(const SimplexId p[3],
                                                 const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
  }
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleVertexC(const SimplexId p[3],
                                                 const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
    else
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + vshift_[0];
  } else {
    if(id == 0)
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
    else
      return (p[0] / 2) + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + vshift_[0];
  }
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleVertexD1(const SimplexId p[3],
                                                  const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1];
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + vshift_[0];
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]
             + vshift_[0];
  }
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleVertexD2(const SimplexId p[3],
                                                  const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1];
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1];
  }
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleVertexD3(const SimplexId p[3],
                                                  const int id) const {
  if(p[0] % 2) {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1;
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1];
  } else {
    if(id == 0)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
    else if(id == 1)
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0];
    else
      return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1];
  }
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeF_0(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 2:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeF_1(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3]
             + 1;
    case 2:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeH_0(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 2:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeH_1(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1];
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5]
             + 1;
    case 2:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeC_0(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3];
    case 1:
      return esetshift_[1] + p[0] / 2 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5];
    case 2:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeC_1(const SimplexId p[3],
                                                 const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3];
    case 1:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5];
    case 2:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeD1_0(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2] + p[2] * eshift_[3]
             + 1;
    case 1:
      return esetshift_[4] + p[0] / 2 + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeD1_1(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[0] + p[0] / 2 + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3];
    case 1:
      return esetshift_[4] + p[0] / 2 + p[1] * eshift_[10] + p[2] * eshift_[11];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeD2_0(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + p[1] * eshift_[0] + p[2] * eshift_[1];
    case 1:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeD2_1(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return p[0] / 2 + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1];
    case 1:
      return esetshift_[3] + p[0] / 2 + p[1] * eshift_[8] + p[2] * eshift_[9]
             + 1;
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeD3_0(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] / 2 + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5];
    case 1:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6] + p[2] * eshift_[7];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleEdgeD3_1(const SimplexId p[3],
                                                  const int id) const {
  switch(id) {
    case 0:
      return esetshift_[1] + p[0] / 2 + p[1] * eshift_[4] + p[2] * eshift_[5]
             + 1;
    case 1:
      return esetshift_[2] + p[0] / 2 + p[1] * eshift_[6]
             + (p[2] + 1) * eshift_[7];
    case 2:
      return esetshift_[5] + p[0] / 2 + p[1] * eshift_[12] + p[2] * eshift_[13];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleLinkF(const SimplexId p[3],
                                               const int id) const {
  if(p[2] > 0 and p[2] < nbvoxels_[2]) {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] - 1) * vshift_[1] + 1;
    }
  } else if(p[2] == 0)
    return p[0] / 2 + (p[1] + 1) * vshift_[0] + vshift_[1];
  else
    return p[0] / 2 + p[1] * vshift_[0] + (p[2] - 1) * vshift_[1] + 1;

  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleLinkH(const SimplexId p[3],
                                               const int id) const {
  if(p[1] > 0 and p[1] < nbvoxels_[1]) {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] - 1) * vshift_[0] + p[2] * vshift_[1] + 1;
    }
  } else if(p[1] == 0)
    return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1];
  else
    return p[0] / 2 + (p[1] - 1) * vshift_[0] + p[2] * vshift_[1] + 1;

  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleLinkC(const SimplexId p[3],
                                               const int id) const {
  if(p[0] > 1 and p[0] < (dimensions_[0] * 2 - 2)) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] - 1;
    }
  } else if(p[0] < 2)
    return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1] + 1;
  else
    return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] - 1;

  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleLinkD1(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1] + 1;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] + 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleLinkD2(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1] + 1;
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1] + 1;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1];
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleLinkD3(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + (p[2] + 1) * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + (p[2] + 1) * vshift_[1] + 1;
    }
  } else {
    switch(id) {
      case 0:
        return p[0] / 2 + p[1] * vshift_[0] + p[2] * vshift_[1];
      case 1:
        return p[0] / 2 + (p[1] + 1) * vshift_[0] + p[2] * vshift_[1] + 1;
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleStarF(const SimplexId p[3],
                                               const int id) const {
  if(p[0] % 2) {
    if(p[2] > 0 and p[2] < nbvoxels_[2]) {
      switch(id) {
        case 0:
          return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
                 + 1; // BCDG
        case 1:
          return (p[0] - 1) * 3 + p[1] * tetshift_[0]
                 + (p[2] - 1) * tetshift_[1] + 4; // BFGH
      }
    } else if(p[2] == 0)
      return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // BCDG
    else
      return (p[0] - 1) * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 4; // BFGH
  } else {
    if(p[2] > 0 and p[2] < nbvoxels_[2]) {
      switch(id) {
        case 0:
          return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
        case 1:
          return p[0] * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
                 + 3; // BEFG
      }
    } else if(p[2] == 0)
      return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
    else
      return p[0] * 3 + p[1] * tetshift_[0] + (p[2] - 1) * tetshift_[1]
             + 3; // BEFG
  }

  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleStarH(const SimplexId p[3],
                                               const int id) const {
  if(p[0] % 2) {
    if(p[1] > 0 and p[1] < nbvoxels_[1]) {
      switch(id) {
        case 0:
          return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
                 + 3; // BEFG
        case 1:
          return (p[0] - 1) * 3 + (p[1] - 1) * tetshift_[0]
                 + p[2] * tetshift_[1] + 5; // BDGH
      }
    } else if(p[1] == 0)
      return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 3; // BEFG
    else
      return (p[0] - 1) * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // BDGH
  } else {
    if(p[1] > 0 and p[1] < nbvoxels_[1]) {
      switch(id) {
        case 0:
          return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
                 + 2; // ABEG
        case 1:
          return p[0] * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
                 + 1; // BCDG
      }
    } else if(p[1] == 0)
      return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2; // ABEG
    else
      return p[0] * 3 + (p[1] - 1) * tetshift_[0] + p[2] * tetshift_[1]
             + 1; // BCDG
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleStarC(const SimplexId p[3],
                                               const int id) const {
  if(p[0] % 2) {
    if(p[0] > 1 and p[0] < (dimensions_[0] * 2 - 2)) {
      switch(id) {
        case 0:
          return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
                 + 2; // ABEG
        case 1:
          return ((p[0] - 2) / 2) * 6 + p[1] * tetshift_[0]
                 + p[2] * tetshift_[1] + 4; // BFGH
      }
    } else if(p[0] < 2)
      return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 2; // ABEG
    else
      return ((p[0] - 2) / 2) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 4; // BFGH
  } else {
    if(p[0] > 1 and p[0] < (dimensions_[0] * 2 - 2)) {
      switch(id) {
        case 0:
          return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
        case 1:
          return ((p[0] - 1) / 2) * 6 + p[1] * tetshift_[0]
                 + p[2] * tetshift_[1] + 5; // BDGH
      }
    } else if(p[0] < 2)
      return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
    else
      return ((p[0] - 1) / 2) * 6 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
             + 5; // BDGH
  }

  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleStarD1(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 2; // ABEG
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 3; // BEFG
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1; // BCDG
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 5; // BDGH
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleStarD2(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 5; // BDGH
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 4; // BFGH
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 2; // ABEG
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTriangleStarD3(const SimplexId p[3],
                                                const int id) const {
  if(p[0] % 2) {
    switch(id) {
      case 0:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 3; // BEFG
      case 1:
        return (p[0] - 1) * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]
               + 4; // BFGH
    }
  } else {
    switch(id) {
      case 0:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1]; // ABCG
      case 1:
        return p[0] * 3 + p[1] * tetshift_[0] + p[2] * tetshift_[1] + 1; // BCDG
    }
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronVertexABCG(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1]; // a
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]; // c
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronVertexBCDG(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]; // c
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1; // d
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronVertexABEG(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1]; // a
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]; // e
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronVertexBEFG(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1]; // e
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1; // f
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronVertexBFGH(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[1] + 1; // f
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1; // h
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronVertexBDGH(const SimplexId p[3],
                                                       const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + 1; // b
    case 1:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0] + 1; // d
    case 2:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1]; // g
    case 3:
      return p[0] + p[1] * vshift_[0] + p[2] * vshift_[1] + vshift_[0]
             + vshift_[1] + 1; // h
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronEdgeABCG(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6]
             + p[2] * eshift_[7]; // bc-D1
    case 4:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronEdgeBCDG(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6]
             + p[2] * eshift_[7]; // bc-D1
    case 4:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11]; // be-D3
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronEdgeABEG(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + p[2] * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + p[0] + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[3] + p[0] + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + p[2] * eshift_[11]; // be-D3
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronEdgeBEFG(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + p[1] * eshift_[0] + (p[2] + 1) * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + p[0] + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + (p[0] + 1) + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6]
             + (p[2] + 1) * eshift_[7]; // bc-D1
    case 4:
      return esetshift_[4] + p[0] + p[1] * eshift_[10]
             + p[2] * eshift_[11]; // be-D3
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronEdgeBFGH(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2]
             + (p[2] + 1) * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + (p[0] + 1) + p[1] * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[2] + p[0] + p[1] * eshift_[6]
             + (p[2] + 1) * eshift_[7]; // bc-D1
    case 4:
      return esetshift_[3] + (p[0] + 1) + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronEdgeBDGH(const SimplexId p[3],
                                                     const int id) const {
  switch(id) {
    case 0:
      return p[0] + (p[1] + 1) * eshift_[0] + (p[2] + 1) * eshift_[1]; // ab-L
    case 1:
      return esetshift_[0] + (p[0] + 1) + p[1] * eshift_[2]
             + p[2] * eshift_[3]; // ac-H
    case 2:
      return esetshift_[1] + (p[0] + 1) + (p[1] + 1) * eshift_[4]
             + p[2] * eshift_[5]; // ae-P
    case 3:
      return esetshift_[3] + (p[0] + 1) + p[1] * eshift_[8]
             + p[2] * eshift_[9]; // ag-D2
    case 4:
      return esetshift_[4] + p[0] + (p[1] + 1) * eshift_[10]
             + p[2] * eshift_[11]; // be-D3
    case 5:
      return esetshift_[5] + p[0] + p[1] * eshift_[12]
             + p[2] * eshift_[13]; // bg-D4
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronTriangleABCG(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1];
    case 1:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronTriangleBCDG(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return p[0] * 2 + p[1] * tshift_[0] + p[2] * tshift_[1] + 1;
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11];
    case 2:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3];
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronTriangleABEG(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3];
    case 1:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9];
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 1;
    case 3:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronTriangleBEFG(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7]
             + 1;
    case 1:
      return tsetshift_[0] + p[0] * 2 + p[1] * tshift_[2] + p[2] * tshift_[3]
             + 1;
    case 2:
      return p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1];
    case 3:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronTriangleBFGH(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 1:
      return tsetshift_[4] + p[0] * 2 + p[1] * tshift_[10] + p[2] * tshift_[11]
             + 1;
    case 2:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 3;
    case 3:
      return p[0] * 2 + p[1] * tshift_[0] + (p[2] + 1) * tshift_[1] + 1;
  }
  return -1;
}

inline ttk::SimplexId
  ttk::ZombieTriangulation::getTetrahedronTriangleBDGH(const SimplexId p[3],
                                                         const int id) const {
  switch(id) {
    case 0:
      return tsetshift_[1] + p[0] * 2 + p[1] * tshift_[4] + p[2] * tshift_[5]
             + 2;
    case 1:
      return tsetshift_[2] + p[0] * 2 + p[1] * tshift_[6] + p[2] * tshift_[7];
    case 2:
      return tsetshift_[3] + p[0] * 2 + p[1] * tshift_[8] + p[2] * tshift_[9]
             + 1;
    case 3:
      return tsetshift_[0] + p[0] * 2 + (p[1] + 1) * tshift_[2]
             + p[2] * tshift_[3] + 1;
  }
  return -1;
}

inline ttk::SimplexId ttk::ZombieTriangulation::getTetrahedronNeighborABCG(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t + 1;
    case 1:
      return t + 2;
    case 2:
      if(p[0] > 0)
        return t - 1;
      else
        return t - tetshift_[1] + 3;
    case 3:
      return t - tetshift_[1] + 3;
  }
  return -1;
}

inline ttk::SimplexId ttk::ZombieTriangulation::getTetrahedronNeighborBCDG(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 4;
    case 2:
      if(p[2] > 0)
        return t - tetshift_[1] + 3;
      else
        return t + tetshift_[0] + 1;
    case 3:
      return t + tetshift_[0] + 1;
  }
  return -1;
}

inline ttk::SimplexId ttk::ZombieTriangulation::getTetrahedronNeighborABEG(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 2;
    case 1:
      return t + 1;
    case 2:
      if(p[0] > 0)
        return t - 4;
      else
        return t - tetshift_[0] - 1;
    case 3:
      return t - tetshift_[0] - 1;
  }
  return -1;
}

inline ttk::SimplexId ttk::ZombieTriangulation::getTetrahedronNeighborBEFG(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 1;
    case 2:
      if(p[1] > 0)
        return t - tetshift_[0] + 2;
      else
        return t + tetshift_[1] - 3;
    case 3:
      return t + tetshift_[1] - 3;
  }
  return -1;
}

inline ttk::SimplexId ttk::ZombieTriangulation::getTetrahedronNeighborBFGH(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t + 1;
    case 2:
      if(p[0] < nbvoxels_[0] - 1)
        return t + 4;
      else
        return t + tetshift_[1] - 3;
    case 3:
      return t + tetshift_[1] - 3;
  }
  return -1;
}

inline ttk::SimplexId ttk::ZombieTriangulation::getTetrahedronNeighborBDGH(
  const SimplexId t, const SimplexId p[3], const int id) const {
  switch(id) {
    case 0:
      return t - 1;
    case 1:
      return t - 4;
    case 2:
      if(p[0] < nbvoxels_[0] - 1)
        return t + 1;
      else
        return t + tetshift_[0] - 2;
    case 3:
      return t + tetshift_[0] - 2;
  }
  return -1;
}

#include <ImplicitPreconditions.h>

/// @endcond
