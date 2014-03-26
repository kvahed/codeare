/*!

  \mainpage  Developer's guide and API reference

  \section dguide Developer's guide

  This is the developer's guide and API reference to codeare, the
  <em>co</em>mmon <em>d</em>ata <em>e</em>xchange and
  <em>re</em>construction platform. If you are looking for the general
  documentation, tutorials etc, please refer to http://codeare.org

  codeare is most readily devided in three major realms:

  \li Matrix template and operators thereof
  \li Algorithmic strategies
  \li Network facility

  In a nutshell, all three parts together allow one to have parts of a
  live or retrospective reconstruction or runtime feedback run outside
  your vendor's reconstrcution chain. It could be, however also used
  to have an external device, let's say a camera, an HID device and
  vice versa, to interact with a running sequence. It could on the
  other hand be used to engange with programming languages / paradigms
  that may not be possible with your vendor's setup; think of using
  python, MATLAB, maple, Mathematica, ... LIVE!

  The matrix template is essentially the template for an ND-array that
  comes with an extensive and extensible toolbox to allow one to
  conentrate entirely on the implementation of the algorithm in
  question. Most of the commands that act on matrices are very much
  similiar if not identlical to the way one deals with matrices in
  python, MATLAB and octave etc. (f.e. <code>A = fft(A)</code>);
 
  Implementationwise the major difference to the mentioned programming
  environments and languages, lies in the use of references rather
  than copies of matrices for memory efficiency and the fact that the
  core implementation is inherently bound to the C++ language
  specifications. So that <code>[u,s,d] = svd (a)</code> is not
  option. The SVD case is handled here as<br/>

 ::svd(const Matrix<T>& IN, Matrix<S>& s, Matrix<T>& U, Matrix<T>& V, const char& jobz = 'N').<br/>

  \li <a href="ug_matrix.html">Matrix</a>
  \li <a href="ug_modules.html">Modules</a>
  \li <a href="ug_networking.html">Network IO</a>
  \li <a href="ug_ice.html">Integration with ICE</a>

*/

