# ifndef __OCL_MATRIX_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_MATRIX_HPP__




  /**************
   ** includes **
   **************/
 
  // CoDEARE
  # include "../../Matrix.hpp"
  
  // ocl
  # include "../oclSettings.hpp"
  # include "../oclConnection.hpp"
  # include "../oclTraits.hpp"



  /************************************
   ** struct: ocl_matrix_type_traits **
   **   (base struct)                **
   ************************************/
  template <class T>
  struct ocl_matrix_type_traits
  {
  
    /* -- */
  
  };
  
  
  /************************************
   ** struct: ocl_matrix_type_traits **
   **   (derived, type: float)       **
   ************************************/
  template <>
  struct ocl_matrix_type_traits <float>
  {
  
    typedef float elem_type;
    
    static
    std::string
    type_name         ( )
    {
      return std::string ("oclMatrix <float>");
    }
  
  };
  
  
  /************************************
   ** struct: ocl_matrix_type_traits **
   **   (derived, type: double)      **
   ************************************/
  template <>
  struct ocl_matrix_type_traits <double>
  {
  
    typedef double elem_type;
    
    static
    std::string 
    type_name         ( )
    {
      return std::string ("oclMatrix <double>");
    }
  
  };
  
  
  /************************************
   ** struct: ocl_matrix_type_traits **
   **   (derived, type: size_t)      **
   ************************************/
  template <>
  struct ocl_matrix_type_traits <size_t>
  {
  
    typedef size_t elem_type;
    
    static
    std::string
    type_name         ( )
    {
      return std::string ("oclMatrix <size_t>");
    }
  
  };


  /************************************
   ** struct: ocl_matrix_type_traits **
   **   (derived, type: size_t)      **
   ************************************/
  template <>
  struct ocl_matrix_type_traits <bool>
  {
  
    typedef bool elem_type;
    
    static
    std::string
    type_name         ( )
    {
      return std::string ("oclMatrix <bool>");
    }
  
  };
  
    


  /*********************
   ** class oclMatrix **
   ** (declaration)   **
   *********************/
  template <class T>
  class oclMatrix : public Matrix <T>
  {



    public: 


      /** ************************************************
       * @name            Constructors and destructors  **
       *                  Constructors and destructors  **
       ** ************************************************/
      //@{ /*************************************************************************************/

  
      /**
       * @brief           Construct with size 1x1
       */
      inline
      oclMatrix           ()
                        : Matrix <T> (),
                          class_name (ocl_matrix_type_traits <T> :: type_name ()),
                          mp_oclData (oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                        Matrix <T> :: Size ()))
      {

        T t;
        Validate (t);
                
      }
      
      
      /**
       * @brief           Construct with dimension array.
       *
       * @param  dim      All 16 dimensions.
       */
      inline
      oclMatrix           (const size_t * dim)
                        : Matrix <T> (dim),
                          class_name (ocl_matrix_type_traits <T> :: type_name ()),
                          mp_oclData (oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                        Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
        
      }


      /**
       * @brief           Construct with dimension and resolution arrays.
       *
       * @param  dim      All 16 dimensions.
       * @param  res      All 16 resolutions.
       */
      inline
      oclMatrix           ( const size_t * dim,
                            const  float * res )
                        : Matrix <T> (dim, res),
                          class_name (ocl_matrix_type_traits <T> :: type_name ()),
                          mp_oclData (oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                        Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
      
      }
      
      
      /**
       * @brief           Construct 2D (square).
       *
       * @param  n        Dimension of rows and columns.
       */
      inline
      oclMatrix           ( const size_t & n )
                         : Matrix <T> (n),
                           class_name (ocl_matrix_type_traits <T> :: type_name ()),
                           mp_oclData (oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                         Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
        
      }
      
      
      /**
       * @brief           Construct 2D (general).
       *
       * @param  rows     Rows.
       * @param  cols     Columns.
       */
      inline
      oclMatrix           ( const size_t & rows,
                            const size_t & cols )
                         : Matrix <T> (rows, cols),
                           class_name (ocl_matrix_type_traits <T> :: type_name ()),
                           mp_oclData (oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                         Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
      
      }
      
      
      /**
       * @brief           Construct 3D volume.
       *
       * @param  rows     Rows.
       * @param  cols     Columns.
       * @param  slices   3rd dimension.
       */
      inline
      oclMatrix           ( const size_t & rows,
                            const size_t & cols,
                            const size_t & slices )
                         : Matrix <T> (rows, cols, slices),
                           class_name (ocl_matrix_type_traits <T> :: type_name ()),
                           mp_oclData (oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                         Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
        
      }
      
      
      /**
       * @brief           Construct 4D volume. (or higher dimensional)
       *
       * @param ...       Refer to definition in Matrix <T>.
       */
      inline
      oclMatrix           ( const size_t & col,
                            const size_t & lin,
                            const size_t & cha,
                            const size_t & set,
                            const size_t & eco = 1,
                            const size_t & phs = 1,
                            const size_t & rep = 1,
                            const size_t & seg = 1,
                            const size_t & par = 1,
                            const size_t & slc = 1,
                            const size_t & ida = 1,
                            const size_t & idb = 1,
                            const size_t & idc = 1,
                            const size_t & idd = 1,
                            const size_t & ide = 1,
                            const size_t & ave = 1 )
                         : Matrix <T> (col, lin, cha, set, eco, phs, rep, seg, par, slc, ida, idb, idc, idd, ide, ave),
                           class_name (ocl_matrix_type_traits <T> :: type_name ()),
                           mp_oclData (oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                         Matrix <T> :: Size ()))
      {
      
        T t;
        Validate (t);
        
      }

    
      /**
       * @brief           Transform from Matrix <T> to oclMatrix <T>
       *
       * @param mat       Matrix to copy.
       */
      inline
      oclMatrix           (const Matrix <T> & mat)
                        : Matrix <T> (mat),
                          class_name (ocl_matrix_type_traits <T> :: type_name ()),
                          mp_oclData (oclOperations <T> :: make_GPU_Obj (& (Matrix <T> :: _M [0]),
                                                                        Matrix <T> :: Size ()))
      {

        T t;
        Validate (t);
        
      }


      /**
       * @brief           Copy constructor.
       *
       * @param  mat      oclMatrix to copy.
       */
      inline
      oclMatrix           (const oclMatrix <T> & mat)
                        : Matrix <T> ((Matrix<T>) mat),
                          class_name (ocl_matrix_type_traits <T> :: type_name ()),
                          mp_oclData (oclOperations <T> :: make_GPU_Obj (& (       Matrix <T> :: _M [0]),
                                                                               Matrix <T> :: Size (),
                                                                                    * mat.mp_oclData,
                                                                        oclDataObject :: COPY_BUFFER ))
      {

        T t;
        Validate (t);

      }

      
      /**
       * @brief           Virtual destructor.
       */
      virtual
      ~oclMatrix          ()
      {

        // delete member oclDataObject (created by oclOperations)
        delete mp_oclData;
    
      }


      //@} /*************************************************************************************/
      
      
      /** **************************************
       * @name            Elementwise access. **
       *                  Elementwise access. **
       ** **************************************/
      //@{ /*************************************************************************************/
      
      
      /**
       * @brief           Copy of p-th element.
       *
       * @param  p        Requested position.
       *
       * @return          Value at p-th position.
       */
      inline
      T
      operator[]          (const size_t & p)
      const;
      
      
      /**
       * @brief           Reference to p-th element.
       *
       * @param  p        Requested position.
       *
       * @return          Reference to value at p-th position.
       */
      inline
      T &
      operator[]          (const size_t & p);
      
      
      /**
       * @brief           Pointer to data starting from p-th (default: 0) element.
       *
       * @param  p        Requested position.
       *
       * @return          Data pointer.
       */
      inline
      const T *
      Data                (const size_t p = 0)
      const;
      
      
      /**
       * @brief           Get data (lvalue).
       *
       * @return          Data.
       */
      inline
      std::valarray <T> &
      Dat                 ( );
      
      
      /**
       * @brief           Get data (rvalue).
       *
       * @return          Data.
       */
      inline
      std::valarray <T>
      Dat                 ( )
      const;
      
      
      /**
       * @brief           Element at position p (rvalue).
       *
       * @param  p        Position.
       *
       * @return          Value at p-th position.
       */
      inline
      T
      At                  (const size_t & p)
      const;
      
      
      /**
       * @brief           Element at position p (lvalue).
       *
       * @param  p        Position.
       *
       * @return          Reference to value at p-th position.
       */
      inline
      T &
      At                  (const size_t & p);
      
      
      /**
       * @brief           Element in (first) slice.
       *
       * @param  x        Column.
       * @param  y        Line.
       *
       * @return          Requested value. (rvalue)
       */
      inline
      T
      At                  (const size_t & x,
                           const size_t & y)
      const;
      
      
      /**
       * @brief           Element in (first) slice.
       *
       * @param  x        Column.
       * @param  y        Line.
       *
       * @return          Reference to requested value. (lvalue)
       */
      inline
      T &
      At                  (const size_t & x,
                           const size_t & y);
      
      
      /**
       * @brief           Element in volume.
       *
       * @param  x        Column.
       * @param  y        Line.
       * @param  z        Slice.
       *
       * @return          Requested value. (rvalue)
       */
      inline
      T
      At                  (const size_t & x,
                           const size_t & y,
                           const size_t & z)
      const;
      
      
      /**
       * @brief           Element in volume.
       *
       * @param  x        Column.
       * @param  y        Line.
       * @param  z        Slice.
       *
       * @return          Reference to requested value. (lvalue)
       */
      inline
      T &
      At                  (const size_t & x,
                           const size_t & y,
                           const size_t & z);
      
      
      /**
       * @brief           Get value in store.
       *
       * @param  col      Column.
       * @param  lin      Line.
       * @param  cha      Channel.
       * @param  set      Set.
       * @param  eco      Echo.
       * @param  phs      Phase.
       * @param  rep      Repetition.
       * @param  seg      Segment.
       * @param  par      Partition.
       * @param  slc      Slice.
       * @param  ida      Free index A.
       * @param  idb      Free index B.
       * @param  idc      Free index C.
       * @param  idd      Free index D.
       * @param  ide      Free index E.
       * @param  ave      Average.
       *
       * @return          Value at requested position. (rvalue)
       */
      inline
      T
      At                  (const size_t & col,
                           const size_t & lin,
                           const size_t & cha,
                           const size_t & set,
                           const size_t & eco,
                           const size_t & phs = 0,
                           const size_t & rep = 0,
                           const size_t & seg = 0,
                           const size_t & par = 0,
                           const size_t & slc = 0,
                           const size_t & ida = 0,
                           const size_t & idb = 0,
                           const size_t & idc = 0,
                           const size_t & idd = 0,
                           const size_t & ide = 0,
                           const size_t & ave = 0)
      const;
      
      
      /**
       * @brief           Get value in store.
       *
       * @param  col      Column.
       * @param  lin      Line.
       * @param  cha      Channel.
       * @param  set      Set.
       * @param  eco      Echo.
       * @param  phs      Phase.
       * @param  rep      Repetition.
       * @param  seg      Segment.
       * @param  par      Partition.
       * @param  slc      Slice.
       * @param  ida      Free index A.
       * @param  idb      Free index B.
       * @param  idc      Free index C.
       * @param  idd      Free index D.
       * @param  ide      Free index E.
       * @param  ave      Average.
       *
       * @return          Value at requested position. (lvalue)
       */
      inline
      T &
      At                  (const size_t & col,
                           const size_t & lin,
                           const size_t & cha,
                           const size_t & set,
                           const size_t & eco,
                           const size_t & phs = 0,
                           const size_t & rep = 0,
                           const size_t & seg = 0,
                           const size_t & par = 0,
                           const size_t & slc = 0,
                           const size_t & ida = 0,
                           const size_t & idb = 0,
                           const size_t & idc = 0,
                           const size_t & idd = 0,
                           const size_t & ide = 0,
                           const size_t & ave = 0);
      
      
      /**
       * @brief           Element at position p.
       *
       * @param  p        Requested position.
       *
       * @return          Value at requested position. (rvalue)
       */
      inline
      T
      operator()          (const size_t & p)
      const;
      
      
      /**
       * @brief           Element at position p.
       *
       * @param  p        Requested position.
       *
       * @return          Reference to value at requested position. (lvalue)
       */
      inline
      T &
      operator()          (const size_t & p);


      /**
       * @brief           Element in (first) slice.
       *
       * @param  x        Column.
       * @param  y        Line.
       *
       * @return          Requested value. (rvalue)
       */
      inline
      T
      operator()          (const size_t & x,
                           const size_t & y)
      const;
      
      
      /**
       * @brief           Element in (first) slice.
       *
       * @param  x        Column.
       * @param  y        Line.
       *
       * @return          Reference to requested value. (lvalue)
       */
      inline
      T &
      operator()          (const size_t & x,
                           const size_t & y);
      
      
      /**
       * @brief           Element in volume.
       *
       * @param  x        Column.
       * @param  y        Line.
       * @param  z        Slice.
       *
       * @return          Requested value. (rvalue)
       */
      inline
      T
      operator()          (const size_t & x,
                           const size_t & y,
                           const size_t & z)
      const;
      
      
      /**
       * @brief           Element in volume.
       *
       * @param  x        Column.
       * @param  y        Line.
       * @param  z        Slice.
       *
       * @return          Reference to requested value. (lvalue)
       */
      inline
      T &
      operator()          (const size_t & x,
                           const size_t & y,
                           const size_t & z);
      
      
      /**
       * @brief           Get value in store.
       *
       * @param  col      Column.
       * @param  lin      Line.
       * @param  cha      Channel.
       * @param  set      Set.
       * @param  eco      Echo.
       * @param  phs      Phase.
       * @param  rep      Repetition.
       * @param  seg      Segment.
       * @param  par      Partition.
       * @param  slc      Slice.
       * @param  ida      Free index A.
       * @param  idb      Free index B.
       * @param  idc      Free index C.
       * @param  idd      Free index D.
       * @param  ide      Free index E.
       * @param  ave      Average.
       *
       * @return          Value at requested position. (rvalue)
       */
      inline
      T
      operator()          (const size_t & col,
                           const size_t & lin,
                           const size_t & cha,
                           const size_t & set,
                           const size_t & eco,
                           const size_t & phs = 0,
                           const size_t & rep = 0,
                           const size_t & seg = 0,
                           const size_t & par = 0,
                           const size_t & slc = 0,
                           const size_t & ida = 0,
                           const size_t & idb = 0,
                           const size_t & idc = 0,
                           const size_t & idd = 0,
                           const size_t & ide = 0,
                           const size_t & ave = 0)
      const;
      
      
      /**
       * @brief           Get value in store.
       *
       * @param  col      Column.
       * @param  lin      Line.
       * @param  cha      Channel.
       * @param  set      Set.
       * @param  eco      Echo.
       * @param  phs      Phase.
       * @param  rep      Repetition.
       * @param  seg      Segment.
       * @param  par      Partition.
       * @param  slc      Slice.
       * @param  ida      Free index A.
       * @param  idb      Free index B.
       * @param  idc      Free index C.
       * @param  idd      Free index D.
       * @param  ide      Free index E.
       * @param  ave      Average.
       *
       * @return          Value at requested position. (lvalue)
       */
      inline
      T &
      operator()          (const size_t & col,
                           const size_t & lin,
                           const size_t & cha,
                           const size_t & set,
                           const size_t & eco,
                           const size_t & phs = 0,
                           const size_t & rep = 0,
                           const size_t & seg = 0,
                           const size_t & par = 0,
                           const size_t & slc = 0,
                           const size_t & ida = 0,
                           const size_t & idb = 0,
                           const size_t & idc = 0,
                           const size_t & idd = 0,
                           const size_t & ide = 0,
                           const size_t & ave = 0);
      
      
      //@} /*************************************************************************************/
      
      
      /** ************************************
       * @name            Friend operators. **
       ** ************************************/
      //@{ /*************************************************************************************/
      
      
      /**
       * @brief           Elementwise multiplication with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m * s.
       */
      inline
      friend
      oclMatrix <T>
      operator*           (const    double     & s,
                           const oclMatrix <T> & m)
      {
      
        return m * s;
      
      }
      
      
      /**
       * @brief           Elementwise multiplication with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m * s.
       */
      inline
      friend
      oclMatrix <T>
      operator*           (const     float     & s,
                           const oclMatrix <T> & m)
      {
      
        return m * s;
      
      }


      /**
       * @brief           Elementwise multiplication with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m * s.
       */
      inline
      friend
      oclMatrix <T>
      operator*           (const     short     & s,
                           const oclMatrix <T> & m)
      {
      
        return m * s;
      
      }


      /**
       * @brief           Elementwise multiplication with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m * s.
       */
      inline
      friend
      oclMatrix <T>
      operator*           (const      long     & s,
                           const oclMatrix <T> & m)
      {
      
        return m * s;
      
      }


      /**
       * @brief           Elementwise multiplication with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m * s.
       */
      inline
      friend
      oclMatrix <T>
      operator*           (const      cxfl     & s,
                           const oclMatrix <T> & m)
      {
      
        return m * s;
      
      }


      /**
       * @brief           Elementwise multiplication with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m * s.
       */
      inline
      friend
      oclMatrix <T>
      operator*           (const      cxdb     & s,
                           const oclMatrix <T> & m)
      {
      
        return m * s;
      
      }
      
      
      /**
       * @brief           Elementwise addition of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m + s.
       */
      inline
      friend
      oclMatrix <T>
      operator+           (const    double     & s,
                           const oclMatrix <T> & m)
      {
      
        return m + s;
      
      }
      
      
      /**
       * @brief           Elementwise addition of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m + s.
       */
      inline
      friend
      oclMatrix <T>
      operator+           (const     float     & s,
                           const oclMatrix <T> & m)
      {
      
        return m + s;
      
      }


      /**
       * @brief           Elementwise addition of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m + s.
       */
      inline
      friend
      oclMatrix <T>
      operator+           (const     short     & s,
                           const oclMatrix <T> & m)
      {
      
        return m + s;
      
      }


      /**
       * @brief           Elementwise addition of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m + s.
       */
      inline
      friend
      oclMatrix <T>
      operator+           (const      long     & s,
                           const oclMatrix <T> & m)
      {
      
        return m + s;
      
      }


      /**
       * @brief           Elementwise addition of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m + s.
       */
      inline
      friend
      oclMatrix <T>
      operator+           (const      cxfl     & s,
                           const oclMatrix <T> & m)
      {
      
        return m + s;
      
      }


      /**
       * @brief           Elementwise addition of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m + s.
       */
      inline
      friend
      oclMatrix <T>
      operator+           (const      cxdb     & s,
                           const oclMatrix <T> & m)
      {
      
        return m + s;
      
      }
      
      
      /**
       * @brief           Elementwise subtraction of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m - s.
       */
      inline
      friend
      oclMatrix <T>
      operator-           (const    double     & s,
                           const oclMatrix <T> & m)
      {
      
        return m - s;
      
      }
      
      
      /**
       * @brief           Elementwise subtraction of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m - s.
       */
      inline
      friend
      oclMatrix <T>
      operator-           (const     float     & s,
                           const oclMatrix <T> & m)
      {
      
        return m - s;
      
      }


      /**
       * @brief           Elementwise subtraction of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m - s.
       */
      inline
      friend
      oclMatrix <T>
      operator-           (const     short     & s,
                           const oclMatrix <T> & m)
      {
      
        return m - s;
      
      }


      /**
       * @brief           Elementwise subtraction of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m - s.
       */
      inline
      friend
      oclMatrix <T>
      operator-           (const      long     & s,
                           const oclMatrix <T> & m)
      {
      
        return m - s;
      
      }


      /**
       * @brief           Elementwise subtraction of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m - s.
       */
      inline
      friend
      oclMatrix <T>
      operator-           (const      cxfl     & s,
                           const oclMatrix <T> & m)
      {
      
        return m - s;
      
      }


      /**
       * @brief           Elementwise subtraction of scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m - s.
       */
      inline
      friend
      oclMatrix <T>
      operator-           (const      cxdb     & s,
                           const oclMatrix <T> & m)
      {
      
        return m - s;
      
      }


      /**
       * @brief           Elementwise division by scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m / s.
       */
      inline
      friend
      oclMatrix <T>
      operator/           (const    double     & s,
                           const oclMatrix <T> & m)
      {
      
        return m / s;
      
      }
      
      
      /**
       * @brief           Elementwise division by scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m / s.
       */
      inline
      friend
      oclMatrix <T>
      operator/           (const     float     & s,
                           const oclMatrix <T> & m)
      {
      
        return m / s;
      
      }


      /**
       * @brief           Elementwise division by scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m / s.
       */
      inline
      friend
      oclMatrix <T>
      operator/           (const     short     & s,
                           const oclMatrix <T> & m)
      {
      
        /* TODO */
        getData ();
      
        return m / s;
      
      }


      /**
       * @brief           Elementwise division by scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m / s.
       */
      inline
      friend
      oclMatrix <T>
      operator/           (const      long     & s,
                           const oclMatrix <T> & m)
      {
      
        /* TODO */
        getData ();
      
        return m / s;
      
      }
      
      
      /**
       * @brief           Elementwise equality with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m == s
       */
      inline
      friend
      oclMatrix <bool>
      operator==          (const         T     & s,
                           const oclMatrix <T> & m)
      {
      
        return m == s;
      
      }


      /**
       * @brief           Elementwise greater/equal with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m >= s
       */
      inline
      friend
      oclMatrix <bool>
      operator>=          (const         T     & s,
                           const oclMatrix <T> & m)
      {
      
        return m >= s;
      
      }
      
      
      /**
       * @brief           Elementwise less/equal with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m <= s
       */
      inline
      friend
      oclMatrix <bool>
      operator<=          (const         T     & s,
                           const oclMatrix <T> & m)
      {
      
        return m <= s;
      
      }


      /**
       * @brief           Elementwise unequality with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m != s
       */
      inline
      friend
      oclMatrix <bool>
      operator!=          (const         T     & s,
                           const oclMatrix <T> & m)
      {
      
        return m != s;
      
      }


      /**
       * @brief           Elementwise greater than with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m > s
       */
      inline
      friend
      oclMatrix <bool>
      operator>          (const         T     & s,
                          const oclMatrix <T> & m)
      {
      
        return m > s;
      
      }


      /**
       * @brief           Elementwise less than with scalar. (lhs)
       *
       * @param  s        Scalar. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m < s
       */
      inline
      friend
      oclMatrix <bool>
      operator<          (const         T     & s,
                          const oclMatrix <T> & m)
      {
      
        return m < s;
      
      }


      /**
       * @brief           Elementwise AND with boolean matrix. (lhs)
       *
       * @param  mb       Boolean matrix. (lhs)
       * @param  m        Matrix. (rhs)
       *
       * @return          m & s
       */
      inline
      friend
      oclMatrix <T>
      operator&          (const oclMatrix <T> & mb,
                          const oclMatrix <T> &  m)
      {
      
        return m & mb;
      
      }
      
      
      //@} /*************************************************************************************/
      
      

      /** ******************************
       * @name            Assignment. **
       ** ******************************/
      //@{ /*************************************************************************************/
      
      
      /**
       * @brief           Assignment operator.
       *
       * @param  mat      Matrix to be assigned.
       *
       * @return          Reference to this-matrix.
       */
      oclMatrix <T> &
      operator=           (const oclMatrix <T> & mat);
      
      
      /**
       * @brief           Assignment of valarray.
       *
       * @param  v        Vallarray to be assigned.
       *
       * @return          Reference to this-matrix.
       */
      oclMatrix <T> &
      operator=           (const std::valarray <T> & v);
      
      
      /**
       * @brief           Elemenwise assignment of scalar.
       *
       * @param  s        Scalar to be assigned.
       *
       * @return          Reference to this-matrix.
       */
      oclMatrix <T> &
      operator=           (const T & s);
      
      
      //@} /*************************************************************************************/
      
      
      
      /**
       * @brief           Resize to  m x n  2D matrix.
       *                  Preserving data while shrinking.
       *                  Adding zeros when growing.
       *
       * @param  m        No. of rows.
       * @param  n        No. of cols.
       */
      inline
      void
      Resize              (const size_t & m,
                           const size_t & n);
      
      
      /**
       * @brief           Purge data and free RAM.
       */
      inline
      void
      Clear               ( );
      
      
      /**
       * @brief           Cast operator.
       *
       * @return          Casted object if possible.
       */
      template <class S>
      operator oclMatrix <S> ( )
      const;
      
      
      
      /** ****************************************
       * @name            Arithmetic operators. **
       ** ****************************************/
      //@{ /*************************************************************************************/

      
      /**
       * @brief           Elementwise addition of two matrices.
       *
       * @param mat       Matrix additive.
       */
      template <class S>
      inline
      oclMatrix <T>
      operator+           (const oclMatrix <S> & mat) const;

      
      /**
       * @brief           Elementwise subtraction of two matrices.
       *
       * @param mat       Matrix substruent.
       */
      template <class S>
      inline
      oclMatrix <T>
      operator-           (const oclMatrix <S> & mat) const;
      
      
      /**
       * @brief           Elementwise increment (uniform).
       *
       * @param inc       Increment.
       */
      template <class S>
      oclMatrix <T> &
      operator+=          (const S & inc);


      /**
       * @brief           Elementwise increment (non uniform).
       *
       * @param inc_mat   Matrix containing increments.
       */
      template <class S>
      oclMatrix <T> &
      operator+=          (const oclMatrix <S> & inc_mat);


      /**
       * @brief           Elementwise decrement (uniform).
       *
       * @param dec       Decrement.
       */
      template <class S>
      oclMatrix <T> &
      operator-=          (const S & dec);
      
      
      /**
       * @brief           Elementwise decrement (non uniform).
       *
       * @param dec_mat   Matrix containing decrements.
       */
      template <class S>
      oclMatrix <T> &
      operator-=          (const oclMatrix <S> & dec_mat);


      /**
       * @brief           Matrix product. (Example: this * M)
       *
       * @param  fac_mat  Factor.
       */
      oclMatrix <T>
      operator->*         (const oclMatrix <T> & fac_mat)
      const;

      
      //@} /*************************************************************************************/
      
      
      
      /** ****************************************
       * @name            Comparison operators. **
       ** ****************************************/
      //@{ /*************************************************************************************/


      /**
       * @brief           Scalar equality.
       *
       * @param  s        Compared scalar.
       *
       * @return          Boolean oclMatrix containing elements' equality with scalar.
       */
      oclMatrix <bool>
      operator==          (const T & s)
      const;
      
      
      /**
       * @brief           Scalar inequality.
       *
       * @param  s        Compared scalar.
       *
       * @return          Boolean oclMatrix containing elements' inequality with scalar.
       */
      oclMatrix <bool>
      operator!=          (const T & s)
      const;
      
      
      /**
       * @brief           Scalar greater comparison.
       *
       * @param  s        Compared scalar.
       *
       * @return          Boolean oclMatrix containing elements' comparisons with scalar.
       */
      oclMatrix <bool>
      operator>           (const T & s)
      const;


      /**
       * @brief           Scalar greater or equal comparison.
       *
       * @param  s        Compared scalar.
       *
       * @return          Boolean oclMatrix containing elements' comparisons with scalar.
       */
      oclMatrix <bool>
      operator>=          (const T & s)
      const;
      
      
      /**
       * @brief           Scalar less comparison.
       *
       * @param  s        Compared scalar.
       *
       * @return          Boolean oclMatrix containing elements' comparisons with scalar.
       */
      oclMatrix <bool>
      operator<           (const T & s)
      const;
      
      
      /**
       * @brief           Scalar less or equal comparison.
       *
       * @param  s        Compared scalar.
       *
       * @return          Boolean oclMatrix containing elements' comparisons with scalar.
       */
      oclMatrix <bool>
      operator<=           (const T & s)
      const;
      
      
      /**
       * @brief           Elementwise equality.
       *
       * @param  mat      Compared matrix.
       *
       * @return          Boolean oclMatrix containing elements' comparisons.
       */
      oclMatrix <bool>
      operator==          (const oclMatrix <T> & mat)
      const;


      /**
       * @brief           Elementwise inequality.
       *
       * @param  mat      Compared matrix.
       *
       * @return          Boolean oclMatrix containing elements' comparisons.
       */
      oclMatrix <bool>
      operator!=          (const oclMatrix <T> & mat)
      const;
      
      
      /**
       * @brief           Elementwise greater comparison.
       *
       * @param  mat      Compared matrix.
       *
       * @return          Boolean oclMatrix containing elements' comparisons.
       */
      oclMatrix <bool>
      operator>           (const oclMatrix <T> & mat)
      const;
      
      
      /**
       * @brief           Elementwise greater or equal comparison.
       *
       * @param  mat      Compared matrix.
       *
       * @return          Boolean oclMatrix containing elements' comparisons.
       */
      oclMatrix <bool>
      operator>=          (const oclMatrix <T> & mat)
      const;


      /**
       * @brief           Elementwise less comparison.
       *
       * @param  mat      Compared matrix.
       *
       * @return          Boolean oclMatrix containing elements' comparisons.
       */
      oclMatrix <bool>
      operator<           (const oclMatrix <T> & mat)
      const;


      /**
       * @brief           Elementwise less or equal comparison.
       *
       * @param  mat      Compared matrix.
       *
       * @return          Boolean oclMatrix containing elements' comparisons.
       */
      oclMatrix <bool>
      operator<=          (const oclMatrix <T> & mat)
      const;


      //@} /*************************************************************************************/
      
      
      /** ****************************************
       * @name            Boolean operators.    **
       ** ****************************************/
      //@{ /*************************************************************************************/
      
      
      /**
       * @brief           
      
      
      //@} /*************************************************************************************/
            
      
      /**
       * @brief           copy relevant data to CPU
       */
      inline
      void
      getData             ()
      const;



    private:


      /**********************
       ** member variables **
       **********************/
      
      // class name of specialised oclMatrix <T>
      const std::string class_name;
      
      // holds data of matrix for calculations on GPU
      oclDataWrapper <T> * mp_oclData;

      
      
      /**********************
       ** static variables **
       **********************/
      
      // verbosity level for operators
      static
      const VerbosityLevel op_v_level = VERB_LOW;
      
      
      /********************************
       ** member functions (private) **
       ********************************/
       
      // allowed element types for an instance of oclMatrix
      void
      Validate            (float  & t)  const {}
      void
      Validate            (size_t & t)  const {}
      void
      Validate            (double & t)  const {}
      void
      Validate            (bool   & t)  const {}


      /*************************************
       ** friend declarations (all typed) **
       *************************************/
      friend class oclMatrix <bool>;
      friend class oclMatrix <float>;
      friend class oclMatrix <double>;
      friend class oclMatrix <size_t>;

      
  };



  /* include function definitions */
  # include "oclMatrix.cpp"
  



# endif __OCL_MATRIX_HPP__
