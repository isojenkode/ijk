/*!
 *  @file ijkdualtable_vert.tpp
 *  @brief Lookup tables for dual isosurface and interval volume vertices.
 *  - Each table entry corresponds to an isosurface vertex 
 *    in an isosurface configuration, i.e. an entry in ISODUAL_TABLE_BASE.
 *  - Version 0.6.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2025 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 3.0 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IJKDUALTABLE_VERT_TPP_
#define _IJKDUALTABLE_VERT_TPP_

#include "ijkdualtable.tpp"
#include "ijklist.tpp"

#include <bitset>
#include <vector>


namespace IJKDUALTABLE {

  // ***************************************************************
  // DUAL TABLE VERTEX INFORMATION
  // ***************************************************************

  /// @brief Information about each of the isosurface/interval volume vertices
  ///   for each table entry.
  template <typename NTYPE, typename VINFO_TYPE>
  class DUAL_TABLE_VERTEX_INFO:
    public IJK::LIST_OF_LISTS<VINFO_TYPE,NTYPE> {

  public:
    typedef VINFO_TYPE VERTEX_INFO_TYPE;
    typedef NTYPE NUMBER_TYPE;

  public:
    DUAL_TABLE_VERTEX_INFO() {};
    template <typename DUAL_TABLE_TYPE>
    DUAL_TABLE_VERTEX_INFO(const DUAL_TABLE_TYPE & table);

    template <typename TI_TYPE>
    NTYPE NumVertices(const TI_TYPE ientry) const
    { return(this->ListLength(ientry)); };

    template <typename DUAL_TABLE_TYPE>
    void Set(const DUAL_TABLE_TYPE & table);

    /// Return reference to vertex info.
    template <typename TI_TYPE>
    const VINFO_TYPE & VertexInfo
    (const TI_TYPE ientry, const NTYPE jv) const
    { return(this->ElementRefConst(ientry, jv)); }

    /// Return non-constant reference to vertex info.
    template <typename TI_TYPE>
    VINFO_TYPE & VertexInfoNC
    (const TI_TYPE ientry, const NTYPE jv)
    { return(this->ElementRef(ientry, jv)); }
  };


  // ***************************************************************
  // CLASS DUAL_TABLE_VERTEX_CONNECTIVITY
  // ***************************************************************

  template <typename T>
  concept HasVertexConnectivityInfo = requires (T t, int arg_int) {
    { t.IntersectsFacet(arg_int) };
  };

  
  /// @brief Connectivity information about each of the
  ///   isosurface/interval volume vertices for each table entry.
  /// - May also include other vertex information.
  template <typename NTYPE, typename VINFO_TYPE>
  requires HasVertexConnectivityInfo<VINFO_TYPE>
  class DUAL_TABLE_VERTEX_CONNECTIVITY:
    public DUAL_TABLE_VERTEX_INFO<NTYPE,VINFO_TYPE> {

  public:
    
    /// @brief Return true if some isosurface edge incident
    ///   on isov crosses facet kfacet.
    template <typename TI_TYPE, typename ISOV_INDEX_TYPE,
              typename FTYPE>
    bool IntersectsFacet
    (const TI_TYPE ientry, const ISOV_INDEX_TYPE jv,
     const FTYPE kfacet) const
    { return this->VertexInfo(ientry, jv).IntersectsFacet(kfacet); }
    
    /// @brief Return list of isosurface vertices incident
    ///   on edges that cross facet kfacet.
    template <typename TI_TYPE, typename VTYPE>
    void GetVerticesIncidentOnEdgesCrossingFacet
    (const TI_TYPE ientry, const int kfacet,
     std::vector<VTYPE> & vlist) const;
  };

  
  // ***************************************************************
  // CLASS VERTEX_CONNECTIVITY_INFO
  // ***************************************************************

  /*!
   *  @brief Class storing vertex connectivity information.
   */
  template <const int NUM_FACET_BITS>
  class VERTEX_CONNECTIVITY_INFO {

  public:

    /*!
     *  @brief intersect_facet[kf] is true if some incident edge 
     *    intersects facet kf.
     */
    std::bitset<NUM_FACET_BITS> intersects_facet;

    /*!
     *  @brief Return number of bits in intersects_facet.
     *  - Number of facets must be at most NumIntersectFacetBits().
     */
    static int NumIntersectFacetBits()
    { return NUM_FACET_BITS; }

    /*!
     *  @brief Return true if some incident edge intersects facet kf.
     */
    bool IntersectsFacet(const int kf) const
    { return(intersects_facet[kf]); }

    /*!
     *  @brief Set all intersects_facet bits to false.
     */
    void ResetIntersectsFacet()
    { intersects_facet.reset(); }

    /*!
     *  @brief Set intersects_facet[kf] top flag.
     */
    void SetIntersectsFacet(const int kf, const bool flag)
    { intersects_facet[kf] = flag; }

    /*!
     *  @brief Print the bitset intersects_facet.
     *  - Mainly used for debugging.
     */
    template <typename OSTREAM_TYPE>
    void PrintIntersectsFacet
    (OSTREAM_TYPE & out, const int num_facets) const
    {
      for (int kf = 0; kf < num_facets; kf++)
        { out << intersects_facet[kf]; }
    };

    /*!
     *  @brief Print the bitset intersects_facet.
     *  - Mainly used for debugging.
     *  - Version that prints strings before/after bitset.
     */
    template <typename OSTREAM_TYPE>
    void PrintIntersectsFacet
    (OSTREAM_TYPE & out, const int num_facets,
     const char * s0, const char * s1) const
    {
      out << s0;
      PrintIntersectsFacet(out, num_facets);
      out << s1;
    };
  };


  /*!
   *  @brief Class storing vertex connectivity information for isodual table 
   *    based on cube (hypercube).
   */
  template <const int NUM_FACET_BITS>
  class VERTEX_CONNECTIVITY_INFO_CUBE:
    public VERTEX_CONNECTIVITY_INFO<NUM_FACET_BITS> {

  public:

    /*!
     *  @brief Return true if vertex is connected to both cube facets 
     *    orthogonal to orth_dir.
     *  - Facet orth_dir and (orth_dir+dimension) are parallel.
     *  @param orth_dir Orthogononal direction.
     *    @pre orth_dir < dimension.
     *  @param dimension Dimension of cube.
     */
    bool IntersectsBothFacetsOrthogonalTo
    (const int orth_dir, const int dimension) const 
    {
      return (this->IntersectsFacet(orth_dir) &&
	      this->IntersectsFacet(orth_dir+dimension));
    }
  };


  // ***************************************************************
  // CLASS VERTEX_CONNECTIVITY_AND_DEGREE_INFO
  // ***************************************************************

  /*!
   *  @brief Class storing vertex connectivity and vertex degree
   *    information.
   *  - PROBABLY REDUNDANT: isodual_table.NumIncidentIsopoly()
   *    already returns this information.
   *  - Vertex degree is number of facets (not edges) incident
   *    on the isosurface vertex.
   */
  template <const int NUM_FACET_BITS, typename DEGREE_TYPE>
  class VERTEX_CONNECTIVITY_AND_DEGREE_INFO:
    public VERTEX_CONNECTIVITY_INFO<NUM_FACET_BITS> {

  public:

    /// @brief Number of facets (not edges) incident
    //     on the isosurface vertex
    DEGREE_TYPE degree;

    int Degree() const
    { return degree; }
  };

  
  /*!
   *  @brief Class storing vertex connectivity and vertex degree
   *    information based on cube (hypercube).
   *  - Vertex degree is number of facets (not edges) incident
   *    on the isosurface vertex.
   */
  template <const int NUM_FACET_BITS, typename DEGREE_TYPE>
  class VERTEX_CONNECTIVITY_AND_DEGREE_INFO_CUBE:
    public VERTEX_CONNECTIVITY_INFO_CUBE<NUM_FACET_BITS> {

  public:

    /// @brief Number of facets (not edges) incident
    //     on the isosurface vertex
    DEGREE_TYPE degree;

    int Degree() const
    { return degree; }
  };  

  
  // ***************************************************************
  // COMPUTE VERTEX DEGREE
  // ***************************************************************

  /*!
   *  @brief Compute number of isosurface polytopes incident on each vertex.
   *  - PROBABLY REDUNDANT:
   *    isodual_table.CountNumIncidentIsoPolyDualToFacetEdges()
   *    already computes this information.
   *  - The field degree for each table entry is set by this routine.
   *  @pre Every table ientry in dual_table_vertex_info has a field degree.
   */
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_dual_isotable_vertex_degrees
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & dual_table_vertex_info)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (NTYPE j = 0; j < table.Entry(it).NumVertices(); j++) 
        { dual_table_vertex_info.VertexInfoNC(it, j).degree = 0; }
        

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.IsBipolar(it, ie)) {
          const NTYPE isov = table.IncidentIsoVertex(it, ie);
          dual_table_vertex_info.VertexInfoNC(it, isov).degree++;
        }
      }
    }

  }


  /*!
   *  @brief Compute number of interval volume polytopes incident on each vertex.   
   *  @pre Every table ientry in dual_table_vertex_info has 
   *    a field num_incident_poly.
   *  - The field num_incident_poly for each table entry is set 
   *    by this routine.
   */
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_ivoldual_table_num_incident_poly
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & dual_table_vertex_info)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (NTYPE j = 0; j < table.Entry(it).NumVertices(); j++) 
        { dual_table_vertex_info.VertexInfoNC(it, j).num_incident_poly = 0; }
        

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.EdgeHasDualIVolPoly(it, ie)) {
          const NTYPE isov0 = table.LowerIncident(it, ie);
          const NTYPE isov1 = table.UpperIncident(it, ie);
          dual_table_vertex_info.VertexInfoNC(it, isov0).num_incident_poly++;
          dual_table_vertex_info.VertexInfoNC(it, isov1).num_incident_poly++;
        }
      }

      for (NTYPE iv = 0; iv < table.NumPolyVertices(); iv++) {
        if (table.IsInIntervalVolume(it, iv)) {
          const NTYPE isov = table.IncidentIVolVertex(it, iv);
          dual_table_vertex_info.VertexInfoNC(it, isov).num_incident_poly++;
        }
      }
    }

  }


  /*!
   *  @brief Compute number of isosurface polytopes incident on each vertex.
   *  - Note: An isosurface polytope is a facet of an interval volume polytope.
   *  @pre Every table ientry in vinfo has a field num_incident_isopoly.
   *  - The field num_incident_isopoly for each table entry is set 
   *    by this routine.
   */
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_ivoldual_table_num_incident_isosurface_poly
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & vinfo)
  {
    typedef typename DUAL_TABLE_TYPE::IVOL_VERTEX_INDEX_TYPE IVOLV_TYPE;
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    IVOLV_TYPE ivolv;

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (IVOLV_TYPE j = 0; j < table.Entry(it).NumVertices(); j++) 
        { vinfo.VertexInfoNC(it, j).num_incident_isopoly = 0; }

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {

        if (table.DoesLowerIsosurfaceIntersectEdge(it, ie, ivolv))
          { vinfo.VertexInfoNC(it, ivolv).num_incident_isopoly++; }

        if (table.DoesUpperIsosurfaceIntersectEdge(it, ie, ivolv))
          { vinfo.VertexInfoNC(it, ivolv).num_incident_isopoly++; }
      }
    }
  }


  // ***************************************************************
  // COMPUTE VERTEX CONNECTIVITY
  // ***************************************************************

  /*!
   *  @brief Compute connectivity across each cube facet.
   *  - Facet kf is orthogonal to axis (kf%dimension).
   *  - If kf < dimension, facet contains cube vertex 0.
   *  - If kf >= dimension, facet contains cube vertex with largest index.
   *  @pre Mesh polytope is a cube (hypercube.)
   *  @tparam DUAL_TABLE_TYPE Dual table type. 
   *    - Usually derived from ISODUAL_TABLE_BASE.
   *  @tparam DUAL_TABLE_VINFO_TYPE Type of dual table storing
   *      vertex information.
   *    - Usually derived from DUAL_TABLE_VERTEX_INFO.
   *  @pre DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE has a
   *    static member function NumIntersectsFacetBits() that
   *    returns the number of intersects_facet[] bits.
   *    - NumIntersectsFacetBits() returns the number of
   *      intersects_facet[] bits.
   *  @pre NumIntersectsFacetBits() <= num facets = 2*dimension.
   *  @pre DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE has a
   *    member function ResetIntersectsFacet() that sets
   *    all the intersects facet bits to false.
   *  @pre DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE has a
   *    member function SetIntersectsFacet(kf,flag) that sets 
   *    the intersects facet bit for facet kf to boolean flag.
   */
  template <typename DUAL_TABLE_TYPE, 
            typename DUAL_TABLE_VINFO_TYPE>
  void compute_dual_cube_isotable_vertex_connectivity
  (const DUAL_TABLE_TYPE & table,
   DUAL_TABLE_VINFO_TYPE & dual_table_vertex_info)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_TABLE_VINFO_TYPE::VERTEX_INFO_TYPE 
      VERTEX_INFO_TYPE;
    typedef typename DUAL_TABLE_VINFO_TYPE::NUMBER_TYPE NTYPE;

    const int dimension = table.Dimension();
    IJK::CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube(dimension);
    IJK::PROCEDURE_ERROR error("compute_dual_isotable_cube_connectivity");

    // Check size of intersets_facet.
    if (VERTEX_INFO_TYPE::NumIntersectFacetBits() < cube.NumFacets()) {
      error.AddMessage
        ("Programming error. Too few bits in intersects_facet[].");
      error.AddMessage
        ("  intersects_facet[] has ", 
         VERTEX_INFO_TYPE::NumIntersectFacetBits(), " bits.");
      error.AddMessage
        ("  intersects_facet[] should have at least ", 
         cube.NumFacets(), " bits");
      error.AddMessage
        ("  for an isosurface lookup table of dimension ",
         dimension, ".");
      throw error;
    }

    for (TABLE_INDEX it = 0; it < table.NumTableEntries(); it++) {

      // Initialize
      for (NTYPE j = 0; j < table.Entry(it).NumVertices(); j++) {
        VERTEX_INFO_TYPE & vertex_info =
          dual_table_vertex_info.VertexInfoNC(it,j);
        vertex_info.ResetIntersectsFacet();
      }

      for (NTYPE ie = 0; ie < table.NumPolyEdges(); ie++) {
        if (table.IsBipolar(it, ie)) {

          const NTYPE isov_index_in_cube = 
            table.IncidentIsoVertex(it, ie);
          VERTEX_INFO_TYPE & vertex_info =
            dual_table_vertex_info.VertexInfoNC
            (it,isov_index_in_cube);

          const NTYPE iv0 = cube.EdgeEndpoint(ie, 0);
          const NTYPE iv1 = cube.EdgeEndpoint(ie, 1);

          for (int d = 0; d < dimension; d++) {
            const long mask = (1L << d);
            if ((mask & iv0) == (mask & iv1)) {
              // Both vertices lie on the same facet.
              if ((mask & iv0) == 0) 
                { vertex_info.SetIntersectsFacet(d, true); }
              else {
                vertex_info.SetIntersectsFacet(d+dimension, true);
              }
            }
          }
        }
      }
    }

  }


  // ***************************************************************
  // CLASS DUAL_TABLE_VERTEX_INFO MEMBER FUNCTIONS
  // ***************************************************************

  template <typename NTYPE, typename VINFO_TYPE>
  template <typename DUAL_TABLE_TYPE>
  DUAL_TABLE_VERTEX_INFO<NTYPE,VINFO_TYPE>::
  DUAL_TABLE_VERTEX_INFO(const DUAL_TABLE_TYPE & table)
  {
    Set(table);
  }

  template <typename NTYPE, typename VINFO_TYPE>
  template <typename DUAL_TABLE_TYPE>
  void DUAL_TABLE_VERTEX_INFO<NTYPE,VINFO_TYPE>::
  Set(const DUAL_TABLE_TYPE & table)
  {
    typedef typename DUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    std::vector<NTYPE> num_isov(table.NumTableEntries());

    for (TABLE_INDEX i = 0; i < table.NumTableEntries(); i++) 
      { num_isov[i] = table.Entry(i).NumVertices(); }
      
    this->CreateLists(num_isov);
  }


  
  // ***************************************************************
  // CLASS DUAL_TABLE_VERTEX_CONNECTIVITY MEMBER FUNCTIONS
  // ***************************************************************


  template <typename NTYPE, typename VINFO_TYPE>
  requires HasVertexConnectivityInfo<VINFO_TYPE>
  template <typename TI_TYPE, typename VTYPE>
  void DUAL_TABLE_VERTEX_CONNECTIVITY<NTYPE,VINFO_TYPE>::
  GetVerticesIncidentOnEdgesCrossingFacet
  (const TI_TYPE ientry, const int kfacet,
   std::vector<VTYPE> & vlist) const
  {
    vlist.clear();
    for (int jv = 0; jv < this->NumVertices(ientry); jv++) {
      if (IntersectsFacet(ientry, jv, kfacet))
        { vlist.push_back(jv); }
    }
  }

}

#endif
