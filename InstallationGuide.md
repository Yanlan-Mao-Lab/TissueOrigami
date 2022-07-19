Software:
QT Creator
- Install Qt 5.14.2 - https://download.qt.io/archive/qt/5.14/5.14.2/

Generic Packages:
Could be installed from Homebrew or apt
- boost
- gsl
- openblas
- To solve sparse matrices in the SimulationOnTheGo mode of the code, the sparse matrix Pardiso (https://pardiso-project.org) is required [1-3]. Please note that the visualisation mode of the code does not need Pardiso.

Packages for Mesh Generation:
   - The matlab code, uses Fast Alpha Hulls [4] for triangulation of 3D surfaces.
   - To triangulate an outline, the code uses Triangle [5].


References:
1. C. Alappat, G. Hager, O. Schenk, J. Thies, A. Basermann, A. Bischop, H. Fehske, G. Wellein, A Recursive Algebraic Coloring Technique for Hardware-Efficient Symmetric Sparse Matrix-Vector Multiplication, ACM Transactions on Parallel Computing, Vol. 7, No. 3(19), 2020.
2. M. Bollhöfer, O. Schenk , R. Janalik, S. Hamm, K. Gullapalli, State-of-The-Art Sparse Direct Solvers, Parallel Algorithms in Computational Science and Engineering. Modeling and Simulation in Science, Engineering and Technology, 1-32, Birkhäuser, 2020.
3. M. Bollhöfer, A. Eftekhari, S. Scheidegger, and O. Schenk Large-Scale Sparse Inverse Covariance Matrix Estimation, SIAM J. Sci. Comput., 41(1), A380–A401, 2019
4. Dylan Muir (2022). Fast Alpha Hulls (alpha shapes in 3d; parfor enabled) (https://www.mathworks.com/matlabcentral/fileexchange/32725-fast-alpha-hulls-alpha-shapes-in-3d-parfor-enabled), MATLAB Central File Exchange. Retrieved July 19, 2022.
5. Jonathan Richard Shewchuk, ``Triangle:  Engineering a 2D Quality Mesh Generator and Delaunay Triangulator,'' in Applied Computational Geometry: Towards Geometric Engineering (Ming C. Lin and Dinesh Manocha, editors), volume 1148 of Lecture Notes in Computer Science, pages 203-222, Springer-Verlag, Berlin, May 1996.  (From the First ACM Workshop on Applied Computational Geometry.)'
