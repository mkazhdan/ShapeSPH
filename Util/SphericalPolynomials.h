#ifndef SPHERICAL_POLYNOMIALS_INCLUDED
#define SPHERICAL_POLYNOMIALS_INCLUDED

enum{ SPHERICAL_POLYNOMIAL_X , SPHERICAL_POLYNOMIAL_Y , SPHERICAL_POLYNOMIAL_Z };
enum{ SPHERICAL_POLYNOMIAL_XX , SPHERICAL_POLYNOMIAL_YY , SPHERICAL_POLYNOMIAL_ZZ , SPHERICAL_POLYNOMIAL_XY , SPHERICAL_POLYNOMIAL_XZ , SPHERICAL_POLYNOMIAL_YZ };

double TrigonometricMonomialIntegral( int A , int B , bool halfCircle );
double SphericalMonomialIntegral( int A ,  int B , int C );

template< class Real > Point< Real , 3 > LinearCoefficients( const FourierKeyS2< Real >& key );
template< class Real > Point< Real , 6 > QuadraticCoefficients( const FourierKeyS2< Real >& key );
template< class Real > SquareMatrix< Real , 3 > QuadraticCoefficientsToMatrix( const Point< Real , 6 >& c );
template< class Real > Point< Real , 6 > QuadraticMatrixToCoefficients( const SquareMatrix< Real , 3 >& m );
template< class Real > Point3D< Real > ConstantAndQuadratic( const FourierKeyS2< Real >& key );

template< class Real > SquareMatrix< Real , 3 > LinearMassMatrix( void );
template< class Real > SquareMatrix< Real , 6 > QuadraticMassMatrix( void );


//////////////////////////////////////
// Spherical Polynomial definitions //
//////////////////////////////////////

inline double TrigonometricMonomialIntegral( int A , int B , bool halfCircle )
{
	if     ( B==0 ) // \cos^A(\theta)
		switch( A )
	{
		case  0: return halfCircle ? M_PI : 2. * M_PI;
		case  1: return 0.;
		default: return double(A-1)/A * TrigonometricMonomialIntegral( A-2 , B , halfCircle );
	}
	else if( A==0 ) // \sin^B(\theta)
		switch( B )
	{
		case  0: return halfCircle ? M_PI : 2. * M_PI;
		case  1: return halfCircle ? 2. : 0;
		default: return double(B-1)/B * TrigonometricMonomialIntegral( A , B-2 , halfCircle );
	}
	else if( B==1 ) // \cos^A(\theta) * \sin(\theta)
		return -1./(A+1) * ( halfCircle ? -2. : 0 );
	else if( A==1 ) // \cos(\theta) * \sin^B(\theta)
		return 0.;
	else            // \cos^A(\theta) * \sin^B(\theta)
		return double(A-1)/(A+B) * TrigonometricMonomialIntegral( A-2 , B , halfCircle );
}
inline double SphericalMonomialIntegral( int A ,  int B , int C )
{
	// \Phi(\theta,\phi) = ( \cos(\theta) * \sin(\phi) , \cos(\phi) , \sin(\theta) * \sin(\phi) ) ; \theta \in [0,2\pi) , \phi \in [0,\pi]
	// \int_{S^2} x^A * y^B * z*C = \int_0^{2\pi} \int_0^\pi x^A * y^B * z^C * \sin(\phi) d\phi d\theta
	//                            = \int_0^{2\pi} \int_0^\pi \cos^A(\theta) * \sin^A(\phi) * \cos^B(\phi) * \sin^C(\theta) * \sin^C(\phi) * \sin(\phi) d\phi d\theta
	//                            = \int_0^{2\pi} \cos^A(\theta) * \sin^C(\theta) d\theta * \int_0^\pi \sin^{A+C+1}(\phi) * \cos^B(\phi) d\phi
	return TrigonometricMonomialIntegral( A , C , false ) * TrigonometricMonomialIntegral( B , A+C+1 , true ) * ( (B%2)==0 ? 1 : -1 );
}
template< class Real >
SquareMatrix< Real , 3 > LinearMassMatrix( void )
{
	int powers[3][3];
	powers[SPHERICAL_POLYNOMIAL_X][0] = 1 , powers[SPHERICAL_POLYNOMIAL_X][1] = 0 , powers[SPHERICAL_POLYNOMIAL_X][2] = 0;
	powers[SPHERICAL_POLYNOMIAL_Y][0] = 0 , powers[SPHERICAL_POLYNOMIAL_Y][1] = 1 , powers[SPHERICAL_POLYNOMIAL_Y][2] = 0;
	powers[SPHERICAL_POLYNOMIAL_Z][0] = 0 , powers[SPHERICAL_POLYNOMIAL_Z][1] = 0 , powers[SPHERICAL_POLYNOMIAL_Z][2] = 1;

	SquareMatrix< Real , 3 > c;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ )
		c(i,j) = Real( SphericalMonomialIntegral( powers[i][0]+powers[j][0] , powers[i][1]+powers[j][1] , powers[i][2]+powers[j][2] ) );

	return c;
}
template< class Real >
SquareMatrix< Real , 6 > QuadraticMassMatrix( void )
{
	int powers[6][3];
	powers[SPHERICAL_POLYNOMIAL_XX][0] = 2 , powers[SPHERICAL_POLYNOMIAL_XX][1] = 0 , powers[SPHERICAL_POLYNOMIAL_XX][2] = 0;
	powers[SPHERICAL_POLYNOMIAL_YY][0] = 0 , powers[SPHERICAL_POLYNOMIAL_YY][1] = 2 , powers[SPHERICAL_POLYNOMIAL_YY][2] = 0;
	powers[SPHERICAL_POLYNOMIAL_ZZ][0] = 0 , powers[SPHERICAL_POLYNOMIAL_ZZ][1] = 0 , powers[SPHERICAL_POLYNOMIAL_ZZ][2] = 2;
	powers[SPHERICAL_POLYNOMIAL_XY][0] = 1 , powers[SPHERICAL_POLYNOMIAL_XY][1] = 1 , powers[SPHERICAL_POLYNOMIAL_XY][2] = 0;
	powers[SPHERICAL_POLYNOMIAL_XZ][0] = 1 , powers[SPHERICAL_POLYNOMIAL_XZ][1] = 0 , powers[SPHERICAL_POLYNOMIAL_XZ][2] = 1;
	powers[SPHERICAL_POLYNOMIAL_YZ][0] = 0 , powers[SPHERICAL_POLYNOMIAL_YZ][1] = 1 , powers[SPHERICAL_POLYNOMIAL_YZ][2] = 1;

	SquareMatrix< Real , 6 > c;
	for( int i=0 ; i<6 ; i++ ) for( int j=0 ; j<6 ; j++ )
		c(i,j) = Real( SphericalMonomialIntegral( powers[i][0]+powers[j][0] , powers[i][1]+powers[j][1] , powers[i][2]+powers[j][2] ) );

	return c;
}
template< class Real >
Point< Real , 3 > LinearCoefficients( const FourierKeyS2< Real >& key )
{
	// Solve for the coefficients of the linear terms in the basis { x , y , z }
	Point< Real , 3 > coefficients;
	// [NOTE] This formulation swaps the y- and z-axes
	// Y[1][0] =  sqrt( 3 / (4*PI) ) * (     y         )
	// Y[1][1] = -sqrt( 3 / (8*PI) ) * ( x +     i * z )
	if( key.bandWidth()>1 )
	{
		coefficients[SPHERICAL_POLYNOMIAL_Y] += Real(      sqrt( 3. / (M_PI*4) ) * key(1,0).r );

		coefficients[SPHERICAL_POLYNOMIAL_X] -= Real( 2. * sqrt( 3. / (M_PI*8) ) * key(1,1).r );
		coefficients[SPHERICAL_POLYNOMIAL_Z] += Real( 2. * sqrt( 3. / (M_PI*8) ) * key(1,1).i );
	}
	return coefficients;
}
template< class Real >
SquareMatrix< Real , 3 > QuadraticCoefficientsToMatrix( const Point< Real , 6 >& c )
{
	SquareMatrix< Real , 3 > m;
	m(0,0) = c[SPHERICAL_POLYNOMIAL_XX] , m(1,1) = c[SPHERICAL_POLYNOMIAL_YY] , m(2,2) = c[SPHERICAL_POLYNOMIAL_ZZ];
	m(1,0) = m(0,1) = c[SPHERICAL_POLYNOMIAL_XY] / 2 , m(2,0) = m(0,2) = c[SPHERICAL_POLYNOMIAL_XZ] / 2 , m(1,2) = m(2,1) = c[SPHERICAL_POLYNOMIAL_YZ] / 2;
	return m;
}
template< class Real >
Point< Real , 6 > QuadraticMatrixToCoefficients( const SquareMatrix< Real , 3 >& m )
{
	Point< Real , 6 > c;

	c[SPHERICAL_POLYNOMIAL_XX] = m(0,0) , c[SPHERICAL_POLYNOMIAL_YY] = m(1,1) , c[SPHERICAL_POLYNOMIAL_ZZ] = m(2,2);
	c[SPHERICAL_POLYNOMIAL_XY] = m(1,0) + m(0,1) , c[SPHERICAL_POLYNOMIAL_XZ] = m(2,0) + m(0,2) , c[SPHERICAL_POLYNOMIAL_YZ] = m(1,2) + m(2,1);
	return c;
}
template< class Real >
Point< Real , 6 > QuadraticCoefficients( const FourierKeyS2< Real >& key )
{
	// Solve for the coefficients of the constant + quadratic terms in the basis { x^2 , y^2 , z^2 , xy , xz , yz }
	Point< Real , 6 > coefficients;
	// [NOTE] This formulation swaps the y- and z-axes
	// Y[0][0] =  sqrt(  1 / ( 4*PI) ) * (       x^2 +    y^2 +     z^2                               )
	// Y[2][0] =  sqrt(  5 / (16*PI) ) * ( -     x^2 + 2* y^2 -     z^2                               )
	// Y[2][1] = -sqrt( 15 / ( 8*PI) ) * (                                     xy           +  i * yz )
	// Y[2][2] =  sqrt( 15 / (32*PI) ) * (       x^2          -     z^2           + 2i * xz           )
	if( key.bandWidth() ) coefficients[SPHERICAL_POLYNOMIAL_XX] = coefficients[SPHERICAL_POLYNOMIAL_YY] = coefficients[SPHERICAL_POLYNOMIAL_ZZ] = Real( sqrt(1./(4.*M_PI) ) ) * key(0,0).r;
	if( key.bandWidth()>2 )
	{
		coefficients[SPHERICAL_POLYNOMIAL_XX] -= Real(      sqrt(  5. / (M_PI*16) ) * key(2,0).r );
		coefficients[SPHERICAL_POLYNOMIAL_YY] += Real( 2. * sqrt(  5. / (M_PI*16) ) * key(2,0).r );
		coefficients[SPHERICAL_POLYNOMIAL_ZZ] -= Real(      sqrt(  5. / (M_PI*16) ) * key(2,0).r );

		coefficients[SPHERICAL_POLYNOMIAL_XY] -= Real( 2. * sqrt( 15. / (M_PI* 8) ) * key(2,1).r );
		coefficients[SPHERICAL_POLYNOMIAL_YZ] += Real( 2. * sqrt( 15. / (M_PI* 8) ) * key(2,1).i );

		coefficients[SPHERICAL_POLYNOMIAL_XX] += Real( 2. * sqrt( 15. / (M_PI*32) ) * key(2,2).r );
		coefficients[SPHERICAL_POLYNOMIAL_ZZ] -= Real( 2. * sqrt( 15. / (M_PI*32) ) * key(2,2).r );
		coefficients[SPHERICAL_POLYNOMIAL_XZ] -= Real( 4. * sqrt( 15. / (M_PI*32) ) * key(2,2).i );
	}
	return coefficients;
}
template< class Real >
Point3D< Real > ConstantAndQuadratic( const FourierKeyS2< Real >& key )
{
	Point< Real , 6 > coeffs = QuadraticCoefficients( key );
	SquareMatrix< Real , 3 > c;
	c(0,0) = coeffs[SPHERICAL_POLYNOMIAL_XX] , c(1,1) = coeffs[SPHERICAL_POLYNOMIAL_YY] , c(2,2) = coeffs[SPHERICAL_POLYNOMIAL_ZZ];
	c(0,1) = c(1,0) = coeffs[SPHERICAL_POLYNOMIAL_XY] , c(0,2) = c(2,0) = coeffs[SPHERICAL_POLYNOMIAL_XZ] , c(1,2) = c(2,1) = coeffs[SPHERICAL_POLYNOMIAL_YZ];

	Real A[3][3] , d[3];
	// [NOTE] To express the coefficient representation as a symmetric-tensor, we need to
	// divide the off-diagonal terms by 2 since they get counted twice.
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) A[i][j] = c(i,j) / ( (i==j)?1:2 );
	eigdc< Real , 3 >( A , d );

	// The values of d now store the coefficients of the factored quadratic representation, aligned to the principal axes in A.
	// However, the basis {x^2,y^2,z^2} is not orthogonal, so we need to transform to an orthonormal basis, of the form { a+b*x^2 , a+b*y^2 , a+b*z^2 }
	// The values {a,b} must satisfy:
	// 1] < a + b*x^2 , a + b*x^2 > = 1
	// 2] < a + b*x^2 , a + b*y^2 > = 0;
	// Or equivalently:
	// 1] a^2 <1,1> + b^2 <x^2,x^2> + 2ab <1,x^2> = 1
	// 2] a^2 <1,1> + b^2 <x^2,y^2> + 2ab <1,x^2> = 0
	// => b^2 <x^2,x^2-y^2> = 1
	// => b = \pm sqrt( 1./ <x^2,x^2-y^2> )
	Real b = Real( sqrt( 1. / ( SphericalMonomialIntegral( 4 , 0 , 0 ) - SphericalMonomialIntegral( 2 , 2 , 0 ) ) ) );
	// Solving: a^2 * <1,1> + 2a * b<1,x^2> + b^2<x^2,y^2> = 0
	Real _a = Real( SphericalMonomialIntegral( 0 , 0 , 0 ) ) , _b = b * Real( SphericalMonomialIntegral( 2 , 0 , 0 ) ) , _c = b * b * Real( SphericalMonomialIntegral( 2 , 2 , 0 ) );
	// We get a = ( -_b +/- sqrt( _b^2 - _a * _c ) ) / _a
	Real disc = _b*_b - _a*_c;
	if( disc<0 )
	{
		fprintf( stderr , "[WARNING] Negative discriminant set to zero\n" );
		disc = 0;
	}
	Real a = ( -_b + Real( sqrt( disc ) ) ) / _a;

	// Finally, to get the coefficients in this basis, we use the fact it is orthonormal:
	Point3D< Real > p;
	for( int i=0 ; i<3 ; i++ )
	{
		int A=(i==0?2:0) , B=(i==1?2:0) , C=(i==2?2:0);
		for( int j=0 ; j<3 ; j++ )
		{
			int _A=(j==0?2:0) , _B=(j==1?2:0) , _C=(j==2?2:0);
			p[i] += Real( a*SphericalMonomialIntegral( _A , _B , _C ) + b*SphericalMonomialIntegral( A+_A , B+_B , C+_C ) )*d[j];
		}
	}
	return p;
}

#endif // SPHERICAL_POLYNOMIALS_INCLUDED