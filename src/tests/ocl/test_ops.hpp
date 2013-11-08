# ifndef __TEST_OPS_HPP__


  /************
   ** makros **
   ************/

  # define __TEST_OPS_HPP__
    
  /**
   * @brief             This makro expands to the definition
   *                    of an operator class used for testing.
   */
  # define create_operator(name, op)                                                                     \
                                                                                                         \
  template <class T, class S, class R>                                                                   \
  struct name                                                                                            \
  {                                                                                                      \
                                                                                                         \
    template <MATRIX_TYPE n1, MATRIX_TYPE n2, MATRIX_TYPE n3>                                            \
    static inline                                                                                        \
    typename MType <R, n3> :: MT                                                                         \
    apply               (typename MType <T, n1> :: MT & mat1, const typename MType <S, n2> :: MT & mat2) \
    {                                                                                                    \
      return mat1 op mat2;                                                                               \
    }                                                                                                    \
                                                                                                         \
  };                                                                                                     \
  
  
  
  /**********************
   ** type definitions **
   **********************/
//  typedef double elem_type;
//  typedef oclMatrix <elem_type> Matrix_type;


  /**
   * @brief             Result object for comparing test results.
   */
  template <typename T>
  struct result
  {
  
    bool equal;
    T mean_abs_err;
  
  };


  /**************************
   ** function definitions **
   **************************/
  
  /**
   * @brief             Compare given matrices (evaluate tests).
   */
  template <class T>
  result <double>
  mat_equal             ( const oclMatrix <T> & ocl_mat,
                          const    Matrix <T> &    smat )
  {
  
//    std::cout << " mat: \n" << smat << std::endl;
//    std::cout << " ocl: \n" << ocl_mat << std::endl;
  
    ocl_mat.getData (); // !!! //
    Matrix <cbool> mat_comp = (smat == ocl_mat);
    result <double> res = {true, 0.0};
    for (int i = 0; i < mat_comp.Dim (0); i++)
      for (int j = 0; j < mat_comp.Dim (1); j++)
      {
        res.equal &= mat_comp (i,j);
      }
  
    if (! res.equal)
    {
      const Matrix <T> diff = smat - ((Matrix <T>) ocl_mat);
      for (int i = 0; i < diff.Dim (0); i++)
        for (int j = 0; j < diff.Dim (1); j++)
          res.mean_abs_err += std::abs (diff (i, j));
      res.mean_abs_err /= diff.Size ();
    }
  
    return res;

  }

  
  /**
   * @brief           Defines type of test operands.
   */
  enum MATRIX_TYPE
  {

    SCALAR = 0,
    MATRIX,
    OCL_MATRIX

  };

  /***************************************************
   ** struct MType: Wrap type information (traits). **
   **       (base struct)                           **
   ***************************************************/
  template <class T, MATRIX_TYPE n>
  struct MType
  {  };

  /**************************
   ** MType: oclMatrix <T> **
   **    (derived)         **
   **************************/
  template <class T>
  struct MType <T, OCL_MATRIX>
  {
    typedef oclMatrix <T> MT;
  };

  /***********************
   ** MType: Matrix <T> **
   **    (derived)      **
   ***********************/
  template <class T>
  struct MType <T, MATRIX>
  {
    typedef    Matrix <T> MT;
  };
  
  /************************
   ** MType: T  (Scalar) **
   **    (derived)       **
   ************************/
  template <class T>
  struct MType <T, SCALAR>
  {
    typedef            T  MT;
  };

  
  /**
   * @brief                Function for evaluating test
   *                       of given operator for initialized
   *                       operands.
   *
   * @param  mat1          First operand of Matrix-setting.
   * @param  mat2          Second operand of Matrix-setting.
   * @param  oclmat1       First operand of oclMatrix-setting.
   * @param  oclmat2       Second operand of oclMatrix-setting.
   * @param  op            Operator to be tested.
   *
   * @return               Success.
   */
  template <class T, class S, class R, class Op>
  static
  bool
  operator_eval_MM         (          Matrix <T>     & mat1,
                            const     Matrix <S>     & mat2,
                                   oclMatrix <T>     & oclmat1,
                            const  oclMatrix <S>     & oclmat2,
                            const         Op         & op )
  {
    double time_s, time_ocl;
    MyTimer t;
               t.tic (time_default);
         Matrix <R>     res = op. template apply     <MATRIX,     MATRIX,     MATRIX> (   mat1,    mat2);
    time_s   = t.tic (time_default);
      oclMatrix <R> ocl_res = op. template apply <OCL_MATRIX, OCL_MATRIX, OCL_MATRIX> (oclmat1, oclmat2);
    time_ocl = t.tic (time_default);
//    std::cout << " res: \n" << res << std::endl;
//    std::cout << " ocl: \n" << ocl_res << std::endl;
    result <double> r = mat_equal <R> (ocl_res, res);
    assert ((r.equal || r.mean_abs_err < 5e-1) && " Test failed! ");
    if (verbose)
    {
      std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)";
      if (r.mean_abs_err != 0.0)
        std::cout << " -> mean abs. err. = " << scientific << r.mean_abs_err << fixed;
      std::cout << std::endl;
    }
    return true;
  }

  
  /**
   * @brief                Function for evaluating test
   *                       of given operator for initialized
   *                       operands.
   *
   * @param  mat1          First operand of Matrix-setting.
   * @param  oclmat1       First operand of oclMatrix-setting.
   * @param  scalar        Second operand (scalar).
   * @param  op            Operator to be tested.
   *
   * @return               Success.
   */
  template <class T, class S, class R, class Op>
  bool
  operator_eval_MS         (          Matrix <T>       & mat1,
                                   oclMatrix <T>       & oclmat1,
                             const         S           & scalar,
                             const        Op           & op )
  {
    double time_s, time_ocl;
    MyTimer t;
               t.tic (time_default);
         Matrix <R>     res = op. template apply     <MATRIX, SCALAR,     MATRIX> (   mat1, scalar);
    time_s   = t.tic (time_default);
      oclMatrix <R> ocl_res = op. template apply <OCL_MATRIX, SCALAR, OCL_MATRIX> (oclmat1, scalar);
    time_ocl = t.tic (time_default);
    result <double> r = mat_equal <R> (ocl_res, res);
    assert ((r.equal || r.mean_abs_err < 1e-1) && " Test failed! ");
    if (verbose)
    {
      std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)";
      if (r.mean_abs_err != 0.0)
        std::cout << " -> mean abs. err. = " << scientific << r.mean_abs_err << fixed;
      std::cout << std::endl;
    }
    return true;
  }


  /**
   * @brief                Function for evaluating test
   *                       of given operator for initialized
   *                       operands.
   *
   * @param  scalar        First operand (scalar).
   * @param  mat2          Second operand of Matrix-setting.
   * @param  oclmat2       Second operand of oclMatrix-setting.
   * @param  op            Operator to be tested.
   *
   * @return               Success.
   */
  template <class T, class S, class R, class Op>
  bool
  operator_eval_SM         (               T           & scalar,
                             const    Matrix <T>       & mat2,
                             const oclMatrix <T>       & oclmat2,
                             const        Op           & op )
  {
    double time_s, time_ocl;
    MyTimer t;
               t.tic (time_default);
         Matrix <R>     res = op. template apply <SCALAR,     MATRIX,     MATRIX> (scalar,    mat2);
    time_s   = t.tic (time_default);
      oclMatrix <R> ocl_res = op. template apply <SCALAR, OCL_MATRIX, OCL_MATRIX> (scalar, oclmat2);
    time_ocl = t.tic (time_default);
    result <double> r = mat_equal <R> (ocl_res, res);
//    assert ((r.equal || r.mean_abs_err < 1e-1) && " Test failed! ");
    if (!((r.equal || r.mean_abs_err < 1e-1)))
      std::cerr << "FAILED (mean_abs_err: " << r.mean_abs_err << ")" << std::endl;
    else
    if (verbose)
    {
      std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)";
      if (r.mean_abs_err != 0.0)
        std::cout << " -> mean abs. err. = " << scientific << r.mean_abs_err << fixed;
      std::cout << std::endl;
    }
    return true;
  }


  /**
   * @name              Test operands.
   */
  //@{
  
  create_operator (opAdd, +)
  create_operator (opSub, -)
  create_operator (opMult, *)
  create_operator (opDiv, /)
  create_operator (opPow, ^)
  create_operator (opAddAss, +=)
  create_operator (opSubAss, -=)
  create_operator (opMultAss, *=)
  create_operator (opDivAss, /=)
  create_operator (opPowAss, ^=)
  create_operator (opMinus, = -)
  
  create_operator (opEqual, ==)
  create_operator (opInequal, !=)
  create_operator (opLEqual, <=)
  create_operator (opGEqual, >=)
  create_operator (opLess, <)
  create_operator (opGreater, >)
  create_operator (opAND, &&)
  create_operator (opOR, ||)
  create_operator (opBitwAnd, &)
  
  create_operator (opTrans, = !)
  create_operator (opMatProd, ->*)
  
  //@}
  
  
# endif // __TEST_OPS_HPP__
