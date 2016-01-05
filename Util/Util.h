/*
Copyright (c) 2013, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/
#ifndef UTIL_INCLUDED
#define UTIL_INCLUDED

#include <string>
#include <list>
#include <algorithm>

#if WIN32
#include <windows.h>
#include "Shlwapi.h"
#include "CmdLineParser.h"
inline std::string baseName(const char *path)
{
    char tmp[256];
//    strncmp(tmp, path, 255);
	strncpy(tmp, path, 255);
    tmp[255] = 0;
//    PathStripPath(tmp);
	PathStripPath( (LPWSTR)tmp );
    return std::string(tmp);
}

inline std::string directoryName(const char *path)
{
    char tmp[256];
//    strncmp(tmp, path, 255);
	strncpy(tmp, path, 255);
    tmp[255] = 0;
//	PathRemoveFileSpec(tmp);
//	PathRemoveFileSpec( (LPWSTR)tmp);
	DirectoryName( tmp );
    return std::string(tmp);
}

#endif // WIN32


#define LINJIE_COMPILE_COMPLIANCE 0

double Time( void );
void DumpOutput( const char* fileName , bool echoStdout , const char* format , ... ) ;
void DumpOutput( const char* fileName , bool echoStdout , char* str , const char* format , ... );
void DumpOutput( bool echoStdout , char* str , const char* format , ... );
int offset_fprintf( FILE* fp , unsigned int str_size , const char* format , ... );

class ProgressBar
{
	int _bins;
	size_t _total;
	size_t _idx;
	const char* _header;
	double _startTime , _previousTime;
public:
	ProgressBar( int bins , size_t total , const char* header );
	~ProgressBar( void );
	void update( bool output=true );
	void print( void );
};


#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif // M_PI
template< class Real >
void HSV2RGB( const Real* _hsv , Real* rgb )
{
	Real hsv[3];
	hsv[0] = _hsv[0];
	hsv[1] = _hsv[1];
	hsv[2] = _hsv[2];
	int i;
	Real f, p, q, t;
	if( hsv[1] == 0 )
	{
		rgb[0] = rgb[1] = rgb[2] = hsv[2];
		return;
	}
	else if( hsv[1]<0 )
	{
		fprintf( stderr , "[Error] Saturation can't be negative\n" );
		return;
	}
	while( hsv[0]<0 ) hsv[0] += 2. * M_PI;
	hsv[0] /= M_PI / 3.;
	i = (int)floor( hsv[0] );
	f = (Real)(hsv[0] - i);
	p = (Real)(hsv[2] * ( 1 - hsv[1] ));
	q = (Real)(hsv[2] * ( 1 - hsv[1] * f ));
	t = (Real)(hsv[2] * ( 1 - hsv[1] * ( 1 - f ) ));
	switch( i ) {
		case 0:
			rgb[0] = hsv[2] , rgb[1] = t , rgb[2] = p;
//			c=Point3D< Real >(hsv[2],t,p);
			break;
		case 1:
			rgb[0] = q , rgb[1] = hsv[2] , rgb[2] = p;
//			c=Point3D< Real >(q,hsv[2],p);
			break;
		case 2:
			rgb[0] = p , rgb[1] = hsv[2] , rgb[2] = t;
//			c=Point3D< Real >(p,hsv[2],t);
			break;
		case 3:
			rgb[0] = p , rgb[1] = q , rgb[2] = hsv[2];
//			c=Point3D< Real >(p,q,hsv[2]);
			break;
		case 4:
			rgb[0] = t , rgb[1] = p , rgb[2] = hsv[2];
//			c=Point3D< Real >(t,p,hsv[2]);
			break;
		default:
			rgb[0] = hsv[2] , rgb[1] = p , rgb[2] = q;
//			c=Point3D< Real >(hsv[2],p,q);
			break;
	}
}
#endif // UTIL_INCLUDED
