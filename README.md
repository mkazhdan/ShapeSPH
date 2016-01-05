<CENTER><H1>Shape Spherical Harmonics (Version1.0)</H1></CENTER>
<CENTER>
<A HREF="#DESCRIPTION">description</A>
<A HREF="#LINKS">links</A>
<A HREF="#NOTES">notes</A>
</CENTER>
<HR>
<A NAME="DESCRIPTION"><B>REPOSITORY DESCRIPTION</B></A><br>
<UL>
This distribution is comprised of three repositories that use signal processing on the sphere and the group of rotations to compute affine-invariant shape descriptors, align models, and compute the reflective and rotational symmetries of a shape.
The efficiency of the code is based on the fast signal processing code <A HREF="http://www.cs.dartmouth.edu/~geelong/sphere/">SpharmonicKit</A> and <A HREF="http://www.cs.dartmouth.edu/~geelong/soft/">SOFT</A> which are included in the distribution.
</UL>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
<A HREF="ShapeSPH.zip">Source Code</A><BR>
Code Description:
<UL>
<LI><A HREF="http://htmlpreview.github.io/?https://github.com/mkazhdan/ShapeSPH/blob/master/descriptors.html">Rotation Invariant Descriptors</A>
<LI><A HREF="http://htmlpreview.github.io/?https://github.com/mkazhdan/ShapeSPH/blob/master/alignment.html">Shape Alignment</A>
<LI> <A HREF="http://htmlpreview.github.io/?https://github.com/mkazhdan/ShapeSPH/blob/master/symmetry.html">Symmetry</A><BR>
</UL>
Papers:
<UL>
<LI><A href="http://www.cs.jhu.edu/~misha/MyPapers/SGP03.pdf">Rotation Invariant Descriptors (SGP 2003)</A>
<LI><A href="http://www.cs.jhu.edu/~misha/MyPapers/SIG04b.pdf">Rotation Alignment (SIGGRAPH 2004)</A>
<LI><A href="http://www.cs.jhu.edu/~misha/MyPapers/SIG04.pdf">Scale Normalization (SIGGRAPH 2004)</A>
<LI><A href="http://www.cs.jhu.edu/~misha/MyPapers/SGP04.pdf">Symmetry Detection (SGP 2004)</A><br>
</UL>

<HR>
<A NAME="NOTES"><B>NOTES</B></A><br>

<UL>
<LI> The executable provided up was compiled for a 64-bit OS and requires the <A HREF="http://www.fftw.org/">FFTW</A> .dlls to run.
<LI> Compilation requires the <A HREF="http://www.fftw.org/">FFTW</A> libraries.
<LI> The <A HREF="http://www.cs.dartmouth.edu/~geelong/sphere/">SpharmonicKit</A> and <A HREF="http://www.cs.dartmouth.edu/~geelong/soft/">SOFT</A> source code is slightly modified from the original to support faster computation on purely real input. 
</UL>
