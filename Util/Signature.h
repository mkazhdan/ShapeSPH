#ifndef SIGNATURE_INCLUDED
#define SIGNATURE_INCLUDED

#include <vector>
#include <stdio.h>
#include <stdlib.h>

template< class Real >
struct Signature
{
	std::vector< Real > values;
	Signature( void ) {;}
	Signature( int sz ){ resize(sz); }
	Real squareNorm( void ) const;
	Real& operator[] ( int idx ) { return values[idx]; }
	const Real& operator[] ( int idx ) const { return values[idx]; }
	int size( void ) const { return int( values.size() ); }
	void resize( int sz ) { values.resize( sz ); }

	bool read( const char* fileName );
	bool write( const char* fileName , bool binary=false ) const;
};
template< class Real >
Real Signature< Real >::squareNorm( void ) const
{
	Real norm2=0;
	for( int i=0 ; i<values.size() ; i++ ) norm2 += values[i]*values[i];
	return norm2;
}
template< class Real >
bool Signature< Real >::read( const char* fileName )
{
	FILE* fp = fopen( fileName , "rb" );
	if( !fp ) 
	{
		fprintf( stderr , "[ERROR] Failed to open file for reading (%s)\n" , fileName );
		return false;
	}
	char line[256];
	char fileType[100];
	int size;
#if 1
	if( !fgets( line , 255 , fp ) )
	{
		fprintf( stderr , "[ERROR] Could not read initial line (%s)\n" , fileName );
		fclose( fp );
		return false;
	}
	if( sscanf( line , " %d %s " , &size , fileType)!=2 )
	{
		fprintf( stderr , "[ERROR] Failed to read size and file-type (%s)\n" , fileName );
		fclose( fp );
		return false;
	}
#else
	if( fscanf( fp , "%d %s\n" , &size , fileType )!=2 )
	{
		fprintf( stderr , "[ERROR] Failed to read size and file-type (%s)\n" , fileName );
		fclose( fp );
		return false;
	}
#endif
	bool binary;
	if     ( !strcasecmp( fileType , "ASCII"  ) ) binary = false;
	else if( !strcasecmp( fileType , "BINARY" ) ) binary = true;
	else
	{
		fprintf( stderr , "[ERROR] Failed to read file-type (%s): %s\n" , fileName , fileType );
		fclose( fp );
		return false;
	}
	values.resize( size );
	if( binary )
	{
		size_t readSize = fread( &values[0] , sizeof( Real ) , size , fp );
		if( readSize!=size )
		{
			fprintf( stderr , "[ERROR] Failed to read all values (%s): %d != %d\n" , fileName , readSize , size );
			fclose( fp );
			return false;
		}
	}
	else
	{
		for( int i=0 ; i<size ; i++ ) if( fscanf( fp , " %f " , &values[i] )!=1 )
		{
			fprintf( stderr , "[ERROR] Failed to read value (%s): %d / %d\n" , fileName , i+1 , size );
			fclose( fp );
			return false;
		}
	}
	fclose( fp );
	return true;
}

template< class Real >
bool Signature< Real >::write( const char* fileName , bool binary=false ) const
{
	FILE* fp = fopen( fileName , "wb" );
	if( !fp )
	{
		fprintf( stderr , "[ERROR] Failed to open file for writing: %s\n" , fileName );
		return false;
	}
	if( binary ) fprintf( fp , "%d BINARY\n" , values.size() ) , fwrite( &values[0] , sizeof( Real ) , values.size() , fp );
	else
	{
		fprintf( fp , "%d ASCII\n" , values.size() );
		for( int i=0 ; i<values.size() ; i++ ) fprintf( fp , " %f" , values[i] );
		fprintf( fp , "\n" );
	}
	fclose( fp );
	return true;
}
#endif // SIGNATURE_INCLUDED