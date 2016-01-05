<CENTER><H1>Shape Spherical Harmonics (Version1.0)</H1></CENTER>
<CENTER>
<A HREF="#DESCRIPTION">description</A>
<A HREF="#LINKS">links</A>
<A HREF="#NOTES">notes</A>
</CENTER>
<HR>
<A NAME="DESCRIPTION"><B>CODE DESCRIPTION</B></A><br>
<UL>
This distribution is comprised of three repositories that use signal processing on the sphere and the group of rotations to compute affine-invariant shape descriptors, align models, and compute the reflective and rotational symmetries of a shape.
The efficiency of the code is based on the fast signal processing code <A HREF="http://www.cs.dartmouth.edu/~geelong/sphere/">SpharmonicKit</A> and <A HREF="http://www.cs.dartmouth.edu/~geelong/soft/">SOFT</A> which are included in the distribution.
</UL>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
<A HREF="ShapeSPH.zip">Source Code</A><BR>
<A HREF="http://htmlpreview.github.io/?https://github.com/mkazhdan/ShapeSPH/blob/master/descriptors.html">Rotation Invariant Descriptors</A><BR>
<A HREF="http://htmlpreview.github.io/?https://github.com/mkazhdan/ShapeSPH/blob/master/alignment.html">Shape Alignment</A><BR>
<A HREF="http://htmlpreview.github.io/?https://github.com/mkazhdan/ShapeSPH/blob/master/symmetry.html">Symmetry</A><BR>
Papers: <A href="http://www.cs.jhu.edu/~misha/MyPapers/SGP03.pdf">Rotation Invariant Descriptors (SGP 2003)</A>, 
<A href="http://www.cs.jhu.edu/~misha/MyPapers/SIG04b.pdf">Rotation Alignment (SIGGRAPH 2004</A>,
<A href="http://www.cs.jhu.edu/~misha/MyPapers/SIG04.pdf">Scale Normalization (SIGGRAPH 2004)</A>,
<A href="http://www.cs.jhu.edu/~misha/MyPapers/SGP04.pdf">Symmetry Detection (SGP 2004</A><br>


<HR>
<A NAME="NOTES"><B>NOTES</B></A><br>

<UL>
<LI> The executables provided are compiled for a 64-bit OS and requires the <A HREF="
<html>
<head>
<title>Spherical Harmonic Shape Descriptors</title>
</head>
<body>
<CENTER><H1>Rotation Invariant Shape Descriptors</H1></CENTER>
<CENTER>
<A HREF="#LINKS">links</A>
<A HREF="#DESCRIPTION">description</A>
<A HREF="#EXECUTABLE">executable</A>
<A HREF="#NOTES">notes</A>
</CENTER>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
<A href="http://www.cs.jhu.edu/~misha/MyPapers/SGP03.pdf">SGP 2003 Paper</A><br>
<A href="ShapeDescriptor.exe">Windows (x64) Executable</A><br>
<A href="../ShapeSPH.zip">Source Code</A><br>
<A href="license.txt">License</A><br>
<HR>
<A NAME="DESCRIPTION"><B>CODE DESCRIPTION</B></A><br>
<UL>
This distribution (part of the larger <B>ShapeSPH</B> package) provides the code for computing rotation invariant shape descriptors for 3D models.
</UL>

<HR>
<A NAME="EXECUTABLE"><B>EXECUTABLE ARGUMENTS (ShapeDescriptor.exe)</B></A><br>
<UL>
<DL>

<DT><b>--in</b> &#60;<i>source triangle mesh</i>&#62;
<DD> This string is the the name of the mesh whose shape descriptor is to be computed. The file is assumed to be in <A  HREF="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> format.

<DT>[<b>--out</b> &#60;<i>output shape descriptor</i>&#62;]
<DD> This optional string is the the name of the file to which the shape descriptor is written.

<DT>[<b>--res</b> &#60;<i>voxel resolution</i>&#62;]
<DD> 
This optional integer specifies the resolution of the voxel grid into which the mesh is rasterized. The default value for this parameter is 64. 

<DT>[<b>--bw</b> &#60;<i>spherical band-width</i>&#62;]
<DD> 
This optional integer specifies the number of spherical frequencies that are used to define the descriptor. The default value for this parameter is 16. 

<DT>[<b>--radii</b> &#60;<i>radial resolution</i>&#62;]
<DD> 
This optional integer specifies the number of spheres over which the 3D Gaussian-EDT is to be sampled. The default value for this parameter is 32. 

<DT>[<b>--aScale</b> &#60;<i>anisotropic scaling iterations</i>&#62;]
<DD> 
This optional integer specifies the number of scaling iterations that should be performed in order find the canonical anisotropic scale of the model. The default value for this parameter is 0, indicating that only isotropic scale.

<DT>[<b>--threads</b> &#60;<i>number of threads</i>&#62;]
<DD> 
This optional integer specifies the number of threads across which the solver should be parallelized. The default value for this parameter is the number of threads on the machine.

<DT>[<b>--radius</b> &#60;<i>moment radius</i>&#62;]
<DD> 
This optional floating point value specifies the multiple of the second-order moment radius that should be used for defining the canonical isotropic scale of a shape. The default value for this parameter is 2.0. 

<DT>[<b>--fallOff</b> &#60;<i>Gaussian fall off</i>&#62;]
<DD> 
This optional floating point value specifies the radius of the Gaussian used to define Gaussian-EDT. The default value for this parameter is 2.828427.

<DT>[<b>--noCQ</B>]
<DD> If this optional argument is specified, the rotation-invariant representation of frequencies 0 and 2 is represented by the norms of the two components. Otherwise the triplet of values representing the (orthonormal) covariance is used.

<DT>[<b>--binary</B>]
<DD> If this optional argument is specified, the shape descriptor is written out in binary.

<DT>[<b>--double</B>]
<DD> If this optional argument is specified, the computation is performed using double-precision arithmetic. Otherwise, single-precision arithmetic is used.

<DT>[<b>--verbose</B>]
<DD> If this optional argument is specified, the computation is run in "verbose" mode, outputting more information about the state of the computation of the shape descriptor.

</DL>
</UL>

<HR>
<A NAME="NOTES"><B>NOTES</B></A><br>

<UL>
<LI> The executable provided up was compiled for a 64-bit OS and requires the <A HREF="http://www.fftw.org/">FFTW</A> .dlls to run.
<LI> Compilation requires the <A HREF="http://www.fftw.org/">FFTW</A> libraries.
<LI> The <A HREF="http://www.cs.dartmouth.edu/~geelong/sphere/">SpharmonicKit</A> and <A HREF="http://www.cs.dartmouth.edu/~geelong/soft/">SOFT</A> source code is slightly modified from the original to support faster computation on purely real input. 
</UL>
