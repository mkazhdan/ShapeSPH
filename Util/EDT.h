#ifndef EDT_INCLUDED
#define EDT_INCLUDED

#include <omp.h>
#include "SignalProcessing/CubeGrid.h"


template< class Real >
void SquaredEDT( const CubeGrid< char >& rasterization , CubeGrid< Real >& edt , int threads=1 );

void SquaredEDT( const CubeGrid< char >& rasterization , CubeGrid< int >& edt , int threads=1 );

template< class Real >
void GaussianEDT( const CubeGrid< char >& rasterization , CubeGrid< Real >& gedt , Real fallOff=Real(sqrt(8.)) , int threads=1 );

///////////////////////////////
// Rasterization definitions //
///////////////////////////////
template< class Real >
void SquaredEDT( const CubeGrid< char >& rasterization , CubeGrid< Real >& edt , int threads )
{
	int res = rasterization.resolution();
	CubeGrid< int > _edt;
	SquaredEDT( rasterization , _edt , threads );
	edt.resize( res );
	const int* _edtPtr = _edt[0];
	Real *edtPtr = edt[0];
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<res*res*res ; i++ ) edtPtr[i] = Real( _edtPtr[i] );
}
template< class Real >
void GaussianEDT( const CubeGrid< char >& rasterization , CubeGrid< Real >& gedt , Real fallOff , int threads )
{
	int res = rasterization.resolution();
	SquaredEDT( rasterization , gedt , threads );
	Real* _gedt = gedt[0];
	fallOff = Real(2.)*fallOff*fallOff;
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<res*res*res; i++ ) _gedt[i] = Real( exp( - _gedt[i] / fallOff) );
}
void SquaredEDT( const CubeGrid< char >& rasterization , CubeGrid< int >& edt , int threads )
{
	threads = std::max< int >( threads , 1 );
	int res = rasterization.resolution();
	if( !res )
	{
		fprintf( stderr , "[WARNING] Cannot compute distance transform of zero resolution rasterization\n" );
		return;
	}
	edt.resize( rasterization.resolution() );

	std::vector< int* > oldBuffer( threads ) , newBuffer( threads );
	for( int i=0 ; i<threads ; i++ ) oldBuffer[i] = new int[res] , newBuffer[i] = new int[res];
	
	// Set the upper bound on the distance values
	{
		int* edtPtr = edt[0];
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<res*res*res ; i++ ) edtPtr[i] = 3 * (res+1) * (res+1);
	}

	// scan along z axis
#pragma omp parallel for num_threads( threads )
	for( int xy=0 ; xy<res*res ; xy++ )
	{
		int x = xy/res , y = xy%res;
		bool first=true;
		int dist = 0;
		int* edtPtr = edt[x] + y*res;
		const char* rasterizationPtr = rasterization[x] + y*res;
		for( int z=0 ; z<res ; z++ )
		{
			if( rasterizationPtr[z] )
			{
				dist = 0;
				first = false;
				edtPtr[z] = 0;
			}
			else if( !first )
			{
				dist++;
				edtPtr[z] = dist*dist;
			}
		}

		// backward scan
		dist = 0;
		first = true;
		for( int z=(res-1) ; z>=0 ; z-- )
		{
			if( rasterizationPtr[z] )
			{
				dist = 0;
				first = false;
				edtPtr[z] = 0;
			}
			else if( !first )
			{
				dist++;
				int square = dist*dist;
				if( square<edtPtr[z] ) edtPtr[z] = square;
			}
		}
	}

	// scan along y axis
#pragma omp parallel for num_threads( threads )
	for( int thread=0 ; thread<threads ; thread++ )
	{
		int *_oldBuffer=oldBuffer[thread] , *_newBuffer=newBuffer[thread];
		for( int xz=(res*res*thread)/threads ; xz<(res*res*(thread+1))/threads ; xz++ )
		{
			int x = xz/res , z = xz%res;
			// forward scan
			int s=0;
			int* edtPtr = edt[x] + z;
			for( int y=0 ; y<res; y++ )
			{
				_oldBuffer[y] = edtPtr[y*res];
				int dist = _oldBuffer[y];
				bool foundCloser=false;
				if( dist )
				{
					for( int t=s ; t<=y ; t++ )
					{
						int new_dist = _oldBuffer[t] + (y - t) * (y - t);
						if( new_dist<=dist ) dist = new_dist , s=t , foundCloser=true;
					}
				}
				if( !foundCloser ) s=y;
				_newBuffer[y] = dist;
			}
			// backward scan
			s=res-1;
			for( int y=res-1 ; y>=0 ; y-- )
			{
				int dist = _newBuffer[y];
				bool foundCloser = false;
				if( dist )
				{
					for( int t=s ; t>y ; t-- )
					{
						int new_dist = _oldBuffer[t] + (y - t) * (y - t);
						if( new_dist<=dist ) dist = new_dist , s=t , foundCloser=true;
					}
					edtPtr[y*res] = dist;
				}
				if( !foundCloser ) s=y;
			}
		}
	}

	// scan along x axis
#pragma omp parallel for num_threads( threads )
	for( int thread=0 ; thread<threads ; thread++ )
	{
		int *_oldBuffer = oldBuffer[thread] , *_newBuffer = newBuffer[thread];
		for( int yz=(res*res*thread)/threads ; yz<(res*res*(thread+1))/threads ; yz++ )
		{
			int y = yz/res , z=yz%res;
			// forward scan
			int s=0;
			int* edtPtr = edt[0] + y*res+z;
			for( int x=0 ; x<res ; x++ )
			{
				int dist = _oldBuffer[x] = edtPtr[x*res*res];

				// If the calculated distance to this point is not zero,
				// start from s and see if you can find something closer.
				bool foundCloser = false;
				if( dist )
				{
					for( int t=s ; t<=x ; t++ )
					{
						// Compute the squared distance that would be obtained if we used the (squared) orthogonal distance to t
						// plus the (squared) parallel distance from t to x
						int new_dist = _oldBuffer[t] + (x - t) * (x - t); // <=> new_dist = _oldBuffer[t] + x*x - 2*t*x + t*t
						// If s has not been updated then: _oldBuffer[t] + (x - t) * (x - t) > _oldBuffer[x] for all t <= x
						// Taking y = x+d (w/ d>0) we get:
						// _oldBuffer[t] + ( y - t ) * ( y - t ) = _oldBuffer[t] + ( x - t + d ) * ( x - t + d )
						//                                       = _oldBuffer[t] + ( x - t ) * ( x - t ) + 2 * d * ( x - t ) + d * d
						//                                       > _oldBuffer[x] + ( y - x ) * ( y - x ) + 2 * d * ( x - t )
						//                                       > _oldBuffer[x] + ( y - x ) * ( y - x )
						// for all t <= x
						// That is, the squared distance through t <= x has to be at least as large as the squared distance through x
						if( new_dist<=dist ) dist = new_dist , s=t , foundCloser=true;
					}
				}
				if( !foundCloser ) s=x;
				_newBuffer[x] = dist;
			}

			// backwards scan
			s = res-1;
			for( int x=res-1 ; x>=0 ; x-- )
			{
				int dist = _newBuffer[x];
				bool foundCloser = false;
				if( dist )
				{
					for( int t=s; t>=x ; t-- )
					{
						int new_dist = _oldBuffer[t] + (x - t) * (x - t);
						if( new_dist<=dist ) dist = new_dist , s=t , foundCloser=true;
					}
					edtPtr[x*res*res] = dist;
				}
				if( !foundCloser ) s=x;
			}
		}
	}

	for( int i=0 ; i<threads ; i++ )
	{
		delete[] oldBuffer[i];
		delete[] newBuffer[i];
	}
}
#endif // EDT_INCLUDED