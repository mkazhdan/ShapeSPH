/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <sys/timeb.h>
#ifndef WIN32
#include <sys/time.h>
#endif // WIN32
#include "Util.h"

//#if 1 //LINJIE_COMPILE_COMPLIANCE
double Time( void )
{
#ifdef WIN32
	struct _timeb t;
	_ftime(&t);
	return double(t.time)+double(t.millitm)/1000.0;
#else // WIN32
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+(double)t.tv_usec/1000000;
#endif // WIN32
}

void DumpOutput( const char* fileName , bool echoStdout , const char* format , ... )
{
	if( fileName )
	{
		FILE* fp=fopen( fileName , "a" );
		va_list args;
		va_start( args , format );
		vfprintf( fp , format , args );
		fclose( fp );
		va_end( args );
	}
	if( echoStdout )
	{
		va_list args;
		va_start( args , format );
		vprintf( format , args );
		va_end( args );
	}
}
void DumpOutput( const char* fileName , bool echoStdout , char* str , const char* format , ... )
{
	if( fileName )
	{
		FILE* fp=fopen( fileName , "a" );
		va_list args;
		va_start( args , format );
		vfprintf( fp , format , args );
		fclose( fp );
		va_end( args );
	}
	if( echoStdout )
	{
		va_list args;
		va_start( args , format );
		vprintf( format , args );
		va_end( args );
	}
	if( str )
	{
		va_list args;
		va_start( args , format );
		vsprintf( str , format , args );
		va_end( args );
		if( str[strlen(str)-1]=='\n') str[strlen(str)-1]=0;
	}
}
void DumpOutput( bool echoStdout , char* str , const char* format , ... )
{
	if( echoStdout )
	{
		va_list args;
		va_start( args , format );
		vprintf( format , args );
		va_end( args );
	}
	if( str )
	{
		va_list args;
		va_start( args , format );
		vsprintf( str , format , args );
		va_end( args );
		if( str[strlen(str)-1]=='\n') str[strlen(str)-1]=0;
	}
}
int offset_fprintf( FILE* fp , unsigned int str_size , const char* format , ... )
{
	char temp[1024];
	va_list args;
	va_start( args , format );
	vsprintf( temp , format , args );
	va_end( args );
	for( int i=0 ; i<( (int)str_size )-( (int)strlen(temp) ) ; i++ ) fprintf( fp , " " );
	if( strlen(temp)<=str_size ) return fprintf( fp , "%s" , temp );
	else                         return fprintf( fp , "%s" , temp + strlen(temp) - str_size );
}
/////////////////
// ProgressBar //
/////////////////
ProgressBar::ProgressBar( int bins , size_t total , const char* header )
{
	_startTime = Time();
	_bins = bins;
	_total = total;
	_header = header;
	_idx = 0;
	_previousTime = -1;
}

void ProgressBar::update( bool output )
{
#if 1
#pragma omp atomic
	_idx++;
	if( output ) print( );
#else
	if( output ) print( );
#pragma omp atomic
	_idx++;
#endif
}
void ProgressBar::print( void )
{
#if 1
	int currentBin = int( (_idx*_bins) / _total );
#else
	int currentBin = int( (_idx*_bins) / (_total-1 ) );
#endif
	double currentTime = Time() - _startTime;
	if( int( currentTime*10 ) != int( _previousTime*10 ) )
	{
		printf( "\r[" );
		for( int i=0 ; i<currentBin ; i++ ) printf( "." );
		for( int i=currentBin ; i<_bins ; i++ ) printf( " " );
		printf( "] %s: %.1f (s)\r" , _header , currentTime );
		_previousTime = currentTime;
	}
}
ProgressBar::~ProgressBar( void )
{
	double currentTime = Time() - _startTime;
	printf( "[" );
	for( int i=0 ; i<_bins ; i++ ) printf( "." );
	printf( "] %s: %.1f (s)\n" , _header , currentTime );
}
