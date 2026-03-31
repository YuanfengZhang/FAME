// ============================================================================
// zstdstream, C++ iostream classes wrapping the zstd compression library.
// Based on gzstream by Deepak Bandyopadhyay, Lutz Kettner
// Adapted for zstd by FAME project
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : zstdstream.h
// Author(s)     : Based on gzstream, adapted for zstd
//
// Standard streambuf implementation following Nicolai Josuttis, "The
// Standard C++ Library".
// ============================================================================

#ifndef ZSTDSTREAM_H
#define ZSTDSTREAM_H 1

// standard C++ with new header file names and std:: namespace
#include <iostream>
#include <fstream>
#include <cstdio>
#include <zstd.h>

#ifdef ZSTDSTREAM_NAMESPACE
namespace ZSTDSTREAM_NAMESPACE {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement zstdstream. See below for user classes.
// ----------------------------------------------------------------------------

class zstdstreambuf : public std::streambuf {
private:
    // Buffer sizes recommended by zstd library
    static const int bufferSize = 64 * 1024;    // output buffer size (64KB)

    FILE*            file;               // raw file handle for compressed file
    ZSTD_DStream*    dstream;            // zstd decompression context
    ZSTD_inBuffer    input;              // input buffer struct for zstd
    ZSTD_outBuffer   output;             // output buffer struct for zstd
    char*            inBuf;              // raw input buffer for compressed data
    char             buffer[bufferSize]; // output buffer for streambuf
    char             opened;             // open/close state of stream
    int              mode;               // I/O mode

    int flush_buffer();
public:
    zstdstreambuf() : file(nullptr), dstream(nullptr), inBuf(nullptr), opened(0), mode(0) {
        setp( buffer, buffer + (bufferSize-1));
        setg( buffer + 4,     // beginning of putback area
              buffer + 4,     // read position
              buffer + 4);    // end position
        input.src = nullptr;
        input.size = 0;
        input.pos = 0;
        output.dst = buffer;
        output.size = bufferSize;
        output.pos = 0;
    }
    int is_open() { return opened; }
    zstdstreambuf* open( const char* name, int open_mode);
    zstdstreambuf* close();
    ~zstdstreambuf() { close(); }

    virtual int     overflow( int c = EOF);
    virtual int     underflow();
    virtual int     sync();
};

class zstdstreambase : virtual public std::ios {
protected:
    zstdstreambuf buf;
public:
    zstdstreambase() { init(&buf); }
    zstdstreambase( const char* name, int open_mode);
    ~zstdstreambase();
    void open( const char* name, int open_mode);
    void close();
    zstdstreambuf* rdbuf() { return &buf; }
};

// ----------------------------------------------------------------------------
// User classes. Use izstdstream analogously to ifstream.
// They read files based on the zstd streaming API.
// Files are compatible with zstd compression.
// ----------------------------------------------------------------------------

class izstdstream : public zstdstreambase, public std::istream {
public:
    izstdstream() : std::istream( &buf) {}
    izstdstream( const char* name, int open_mode = std::ios::in)
        : zstdstreambase( name, open_mode), std::istream( &buf) {}
    zstdstreambuf* rdbuf() { return zstdstreambase::rdbuf(); }
    void open( const char* name, int open_mode = std::ios::in) {
        zstdstreambase::open( name, open_mode);
    }
};

// Note: Output stream (ozstdstream) is not implemented as FAME only needs input

#ifdef ZSTDSTREAM_NAMESPACE
} // namespace ZSTDSTREAM_NAMESPACE
#endif

#endif // ZSTDSTREAM_H
// ============================================================================
// EOF //